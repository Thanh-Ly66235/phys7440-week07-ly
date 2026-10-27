"""AlphaFold sequence-to-structure MCP server (PHYS 7440, Assignment 7).

Author: Thanh Ly (University of Manitoba)

Two read-only MCP tools that connect a nucleotide sequence to the AlphaFold
Protein Structure Database (AFDB):

``find_alphafold_models_for_dna``
    Translate a coding DNA/RNA sequence (longest ORF, or a chosen frame), look
    the protein up in AFDB by exact sequence (MD5 checksum), and annotate every
    hit with UniProt organism / review status.  Answers: "is there already a
    predicted structure for the protein this DNA encodes?"

``compare_sequence_to_alphafold``
    Fetch the AFDB model for a UniProt accession, summarise its per-residue
    confidence (pLDDT), and - if a query protein or DNA is supplied - align the
    query to the model sequence and report identity, coverage, and every
    substitution/indel with the local pLDDT and the AlphaMissense pathogenicity
    score at that position.  Answers: "how far can I trust the AlphaFold model
    as a structure for *my* sequence, and where do my variants fall?"

Every call returns the same JSON envelope (``ok``, ``tool``, ``data``,
``error``, ``warnings``, ``provenance``).  Failures never raise out of a tool:
they come back as ``ok=false`` with a machine-readable error ``code`` and the
MCP ``isError`` flag set.

Run ``python mcp_server_ly.py --help`` for transports and the offline replay
mode.  Environment variables are documented in README.md.
"""

import argparse
import dataclasses
import hashlib
import json
import logging
import os
import re
import statistics
import sys
import time
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Any, Literal
from urllib.parse import urlsplit

import anyio
import httpx
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.tools.base import Tool
from mcp.types import CallToolResult, TextContent, ToolAnnotations
from pydantic import BaseModel, ConfigDict, Field, ValidationError

SERVER_NAME = "afdb-sequence-mcp"
SERVER_VERSION = "1.0.0"

log = logging.getLogger(SERVER_NAME)
logging.getLogger("httpx").setLevel(logging.WARNING)  # per-request INFO lines would flood the MCP log

# --------------------------------------------------------------------------------------
# Configuration (environment variables; see README.md)
# --------------------------------------------------------------------------------------


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, default))
    except ValueError:
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except ValueError:
        return default


@dataclasses.dataclass(frozen=True)
class Settings:
    afdb_base: str = "https://alphafold.ebi.ac.uk"
    uniprot_base: str = "https://rest.uniprot.org"
    http_timeout: float = 20.0  # seconds per request
    max_retries: int = 2  # extra attempts on timeout / 5xx / 429
    backoff: float = 0.75  # seconds, doubled on every retry
    cache_ttl: float = 3600.0  # seconds; 0 disables the in-memory response cache
    contact: str = ""  # optional e-mail added to the User-Agent (polite API use)
    replay_dir: str = ""  # serve recorded responses instead of the network

    @classmethod
    def from_env(cls) -> "Settings":
        return cls(
            afdb_base=os.environ.get("AFMCP_AFDB_BASE", cls.afdb_base).rstrip("/"),
            uniprot_base=os.environ.get("AFMCP_UNIPROT_BASE", cls.uniprot_base).rstrip("/"),
            http_timeout=_env_float("AFMCP_HTTP_TIMEOUT", cls.http_timeout),
            max_retries=max(0, _env_int("AFMCP_MAX_RETRIES", cls.max_retries)),
            backoff=max(0.0, _env_float("AFMCP_BACKOFF", cls.backoff)),
            cache_ttl=max(0.0, _env_float("AFMCP_CACHE_TTL", cls.cache_ttl)),
            contact=os.environ.get("AFMCP_CONTACT", ""),
            replay_dir=os.environ.get("AFMCP_REPLAY_DIR", ""),
        )


SETTINGS = Settings.from_env()

# Tests (or advanced users) may inject any httpx transport here, e.g. httpx.MockTransport.
TRANSPORT_OVERRIDE: httpx.AsyncBaseTransport | None = None

# --------------------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------------------

ErrorCode = Literal[
    "INVALID_INPUT",  # argument missing / wrong type / out of range / bad sequence
    "NO_ORF_FOUND",  # DNA contains no open reading frame of the requested length
    "NOT_FOUND",  # AFDB has no model for the accession
    "UPSTREAM_TIMEOUT",  # external API did not answer within AFMCP_HTTP_TIMEOUT (after retries)
    "UPSTREAM_UNAVAILABLE",  # connection failure or HTTP 5xx (after retries)
    "RATE_LIMITED",  # HTTP 429 from the external API (after retries)
    "UPSTREAM_BAD_RESPONSE",  # external API answered with something we cannot parse
    "INTERNAL_ERROR",  # bug in this server (still reported, never crashes the session)
]

RETRYABLE: dict[str, bool] = {
    "INVALID_INPUT": False,
    "NO_ORF_FOUND": False,
    "NOT_FOUND": False,
    "UPSTREAM_TIMEOUT": True,
    "UPSTREAM_UNAVAILABLE": True,
    "RATE_LIMITED": True,
    "UPSTREAM_BAD_RESPONSE": False,
    "INTERNAL_ERROR": False,
}


class ToolFailure(Exception):
    """An anticipated failure that becomes a structured error envelope."""

    def __init__(self, code: str, message: str, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details or {}


# --------------------------------------------------------------------------------------
# Input schemas (the JSON Schemas advertised to clients are generated from these models)
# --------------------------------------------------------------------------------------

UNIPROT_ACCESSION_PATTERN = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$"
Frame = Literal["+1", "+2", "+3", "-1", "-2", "-3"]
FRAMES: tuple[Frame, ...] = ("+1", "+2", "+3", "-1", "-2", "-3")


class FindModelsInput(BaseModel):
    """Arguments of find_alphafold_models_for_dna."""

    model_config = ConfigDict(extra="forbid", title="find_alphafold_models_for_dna")

    dna_sequence: str = Field(
        min_length=30,
        max_length=100_000,
        description=(
            "Coding DNA/cDNA/mRNA in the 5'->3' direction (a CDS works best). Plain text or a single FASTA "
            "record; whitespace, digits and lower case are accepted, U is read as T, IUPAC ambiguity codes are "
            "allowed. Genomic DNA with introns is NOT spliced - pass the CDS."
        ),
        examples=["ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG..."],
    )
    reading_frame: Literal["auto", "+1", "+2", "+3", "-1", "-2", "-3"] = Field(
        default="auto",
        description="'auto' = longest ATG-initiated ORF over all six frames; or force one frame (+1..+3 forward, -1..-3 reverse complement).",
    )
    min_protein_length: int = Field(
        default=30, ge=10, le=5000, description="Shortest ORF (amino acids, excluding stop) that counts as a protein."
    )
    taxon_id: int | None = Field(
        default=None,
        ge=1,
        description="Optional NCBI taxonomy id (e.g. 9606 = human) used to rank hits from that organism first.",
    )
    include_complexes: bool = Field(
        default=False, description="Also return AFDB multimer / complex models that contain the protein."
    )
    max_results: int = Field(default=10, ge=1, le=50, description="Maximum number of AFDB models to return.")


class CompareInput(BaseModel):
    """Arguments of compare_sequence_to_alphafold."""

    model_config = ConfigDict(extra="forbid", title="compare_sequence_to_alphafold")

    uniprot_accession: str = Field(
        pattern=UNIPROT_ACCESSION_PATTERN,
        description="UniProtKB accession of the AlphaFold model, optionally with an isoform suffix (P68871, P04637-2).",
        examples=["P68871"],
    )
    query_sequence: str | None = Field(
        default=None,
        max_length=100_000,
        description=(
            "Optional sequence to compare with the model: one-letter protein, or coding DNA if query_type='dna'. "
            "FASTA headers and whitespace are ignored. Omit to get only the model's confidence profile."
        ),
    )
    query_type: Literal["protein", "dna"] = Field(
        default="protein", description="How to read query_sequence. 'dna' translates the longest ORF first."
    )
    plddt_threshold: float = Field(
        default=70.0, ge=0, le=100, description="pLDDT cut-off used to call 'confident' segments (AFDB convention: 70)."
    )
    include_per_residue_plddt: bool = Field(
        default=False, description="Return the full per-residue pLDDT array (long for big proteins)."
    )
    include_alphamissense: bool = Field(
        default=True,
        description="Look up AlphaMissense pathogenicity for each substitution (available for human canonical proteins).",
    )


# --------------------------------------------------------------------------------------
# Output schemas
# --------------------------------------------------------------------------------------


class ErrorInfo(BaseModel):
    code: ErrorCode = Field(description="Machine-readable error class.")
    message: str = Field(description="Human/LLM-readable explanation and what to do next.")
    retryable: bool = Field(description="True if the same call may succeed later (network/upstream problems).")
    details: dict[str, Any] = Field(default_factory=dict, description="Extra context, e.g. per-field problems.")


class Provenance(BaseModel):
    server: str = SERVER_NAME
    server_version: str = SERVER_VERSION
    sources: list[str] = Field(default_factory=list, description="External URLs queried for this answer.")
    retrieved_at: str = Field(description="UTC timestamp of the call (ISO 8601).")
    elapsed_ms: int
    replay_mode: bool = Field(description="True when answers come from recorded fixtures (offline demo mode).")


class QuerySummary(BaseModel):
    input_characters: int
    nucleotide_length: int
    gc_fraction: float
    ambiguous_bases: int
    was_rna: bool


class OrfCandidate(BaseModel):
    reading_frame: Frame
    strand: Literal["+", "-"]
    nt_start: int = Field(description="1-based start on the submitted (forward) sequence.")
    nt_end: int = Field(description="1-based end on the submitted sequence (includes the stop codon if present).")
    protein_length: int
    has_start_codon: bool
    has_stop_codon: bool


class Translation(OrfCandidate):
    protein_sequence: str
    protein_md5: str = Field(description="MD5 of the protein sequence - the key AFDB uses for exact lookups.")
    unknown_residues: int = Field(description="Number of X residues produced by ambiguous codons.")


class ModelMatch(BaseModel):
    model_id: str
    uniprot_accession: str | None
    protein_name: str | None
    gene: str | None
    organism: str | None
    taxon_id: int | None
    reviewed: bool | None = Field(description="True = Swiss-Prot (curated), False = TrEMBL, None = unknown.")
    taxon_match: bool | None
    oligomeric_state: str | None
    mean_plddt: float | None
    sequence_identity: float | None
    coverage: float | None
    model_url: str | None
    model_page_url: str | None
    created: str | None


class FindModelsData(BaseModel):
    query: QuerySummary
    translation: Translation
    alternative_orfs: list[OrfCandidate] = Field(description="Next-longest ORFs (other frames), for sanity checks.")
    match_status: Literal["exact_match", "no_exact_match"]
    total_matches: int = Field(description="Number of AFDB models with an identical sequence (after filtering).")
    matches: list[ModelMatch]
    best_match: ModelMatch | None
    next_step: str


class ModelUrls(BaseModel):
    entry_page: str
    pdb: str | None
    cif: str | None
    plddt_json: str | None
    pae_json: str | None
    pae_image: str | None
    alphamissense_csv: str | None


class ModelInfo(BaseModel):
    entry_id: str
    uniprot_accession: str
    uniprot_id: str | None
    gene: str | None
    protein_name: str | None
    organism: str | None
    taxon_id: int | None
    reviewed: bool | None
    model_version: int | None
    model_created: str | None
    tool_used: str | None
    provider: str | None
    sequence_length: int
    model_start: int
    model_end: int
    afdb_reported_mean_plddt: float | None
    urls: ModelUrls
    other_models_for_accession: list[str]


class Segment(BaseModel):
    start: int
    end: int
    length: int
    mean_plddt: float


class ConfidenceSummary(BaseModel):
    mean_plddt: float
    median_plddt: float
    min_plddt: float
    max_plddt: float
    fraction_very_high: float = Field(description="pLDDT >= 90")
    fraction_confident: float = Field(description="70 <= pLDDT < 90")
    fraction_low: float = Field(description="50 <= pLDDT < 70")
    fraction_very_low: float = Field(description="pLDDT < 50 (often intrinsically disordered)")
    threshold: float
    fraction_above_threshold: float
    confident_segments: list[Segment] = Field(description=">= 8 consecutive residues with pLDDT >= threshold.")
    disordered_segments: list[Segment] = Field(description=">= 8 consecutive residues with pLDDT < 50.")
    interpretation: str


PlddtBand = Literal["very_high", "confident", "low", "very_low"]


class Substitution(BaseModel):
    variant: str = Field(description="Model residue, model position, query residue, e.g. E7V.")
    model_position: int
    query_position: int
    reference_residue: str
    query_residue: str
    blosum62: int | None
    plddt: float | None
    plddt_band: PlddtBand | None
    alphamissense_pathogenicity: float | None
    alphamissense_class: Literal["likely_benign", "ambiguous", "likely_pathogenic"] | None


class Indel(BaseModel):
    kind: Literal["insertion", "deletion"]
    model_start: int = Field(description="Deletion: first missing model residue. Insertion: model residue before it.")
    model_end: int
    query_start: int
    query_end: int
    length: int
    sequence: str
    local_mean_plddt: float | None


class Comparison(BaseModel):
    query_source: Literal["protein", "dna"]
    query_length: int
    translation: Translation | None
    alignment_method: str
    alignment_score: float
    aligned_pairs: int
    identical: int
    identity: float = Field(description="identical / aligned pairs")
    similarity: float = Field(description="pairs with BLOSUM62 > 0 / aligned pairs")
    query_coverage: float
    model_coverage: float
    n_substitutions: int
    n_insertions: int
    n_deletions: int
    substitutions: list[Substitution]
    indels: list[Indel]
    lists_truncated: bool
    mean_plddt_aligned_region: float | None
    mean_plddt_at_substitutions: float | None
    structure_transfer_confidence: Literal["high", "medium", "low"]
    interpretation: str


class CompareData(BaseModel):
    model: ModelInfo
    confidence: ConfidenceSummary
    per_residue_plddt: list[float] | None
    comparison: Comparison | None
    next_step: str


class _Envelope(BaseModel):
    ok: bool = Field(description="True on success; False means 'error' is set and 'data' is null.")
    tool: str
    warnings: list[str] = Field(default_factory=list)
    error: ErrorInfo | None = None
    provenance: Provenance


class FindModelsResponse(_Envelope):
    """Result envelope of find_alphafold_models_for_dna."""

    model_config = ConfigDict(title="find_alphafold_models_for_dna result")
    data: FindModelsData | None = None


class CompareResponse(_Envelope):
    """Result envelope of compare_sequence_to_alphafold."""

    model_config = ConfigDict(title="compare_sequence_to_alphafold result")
    data: CompareData | None = None


# --------------------------------------------------------------------------------------
# HTTP layer: one client per tool call, retries, small TTL cache, offline replay
# --------------------------------------------------------------------------------------


class ReplayTransport(httpx.AsyncBaseTransport):
    """Serve recorded responses described by <dir>/manifest.json (offline demo + tests)."""

    def __init__(self, directory: str | Path):
        self.directory = Path(directory)
        manifest = json.loads((self.directory / "manifest.json").read_text())
        self.routes: list[dict[str, Any]] = manifest["routes"]

    def _match(self, request: httpx.Request) -> dict[str, Any] | None:
        best, best_score = None, -1
        query = dict(request.url.params)
        for route in self.routes:
            if route.get("host") and route["host"] != request.url.host:
                continue
            if "path" in route and route["path"] != request.url.path:
                continue
            if "path_prefix" in route and not request.url.path.startswith(route["path_prefix"]):
                continue
            wanted = route.get("query", {})
            if any(query.get(k) != v for k, v in wanted.items()):
                continue
            score = len(wanted) + (2 if "path" in route else 0)
            if score > best_score:
                best, best_score = route, score
        return best

    async def handle_async_request(self, request: httpx.Request) -> httpx.Response:
        route = self._match(request)
        if route is None:
            raise httpx.ConnectError(f"replay mode: no recorded response for {request.url}", request=request)
        if "file" in route:
            content = (self.directory / route["file"]).read_bytes()
        else:
            content = json.dumps(route.get("body", {})).encode()
        headers = {"content-type": route.get("content_type", "application/json")}
        return httpx.Response(route.get("status", 200), content=content, headers=headers, request=request)


_CACHE: "OrderedDict[str, tuple[float, int, bytes]]" = OrderedDict()
_CACHE_MAX = 256


class Fetcher:
    """Per-call HTTP helper that records provenance and turns transport problems into ToolFailure."""

    def __init__(self, settings: Settings):
        self.settings = settings
        self.sources: list[str] = []
        transport = TRANSPORT_OVERRIDE
        if transport is None and settings.replay_dir:
            manifest = Path(settings.replay_dir) / "manifest.json"
            if not manifest.is_file():
                raise ToolFailure(
                    "UPSTREAM_UNAVAILABLE",
                    f"Replay mode is on but {manifest} does not exist. Copy the fixtures/ folder next to "
                    "mcp_server_ly.py, or start the server without --replay / AFMCP_REPLAY_DIR to use the live APIs.",
                    {"replay_dir": settings.replay_dir},
                )
            transport = ReplayTransport(settings.replay_dir)
        agent = f"{SERVER_NAME}/{SERVER_VERSION} (+PHYS7440 MCP project; python-httpx)"
        if settings.contact:
            agent += f" mailto:{settings.contact}"
        self.client = httpx.AsyncClient(
            transport=transport,
            timeout=httpx.Timeout(settings.http_timeout, connect=min(10.0, settings.http_timeout)),
            headers={"User-Agent": agent},
            follow_redirects=True,
        )
        self.replay = transport is not None and isinstance(transport, ReplayTransport)

    async def __aenter__(self) -> "Fetcher":
        return self

    async def __aexit__(self, *exc: Any) -> None:
        await self.client.aclose()

    async def get(self, url: str, params: dict[str, Any] | None = None, accept: str = "application/json") -> tuple[int, bytes]:
        full = str(httpx.URL(url, params=params or {}))
        self.sources.append(full)
        now = time.monotonic()
        use_cache = self.settings.cache_ttl > 0 and TRANSPORT_OVERRIDE is None and not self.replay
        if use_cache and full in _CACHE:
            expires, status, body = _CACHE[full]
            if expires > now:
                return status, body
            _CACHE.pop(full, None)

        host = urlsplit(url).netloc
        last: ToolFailure | None = None
        for attempt in range(self.settings.max_retries + 1):
            try:
                response = await self.client.get(url, params=params, headers={"Accept": accept})
            except httpx.TimeoutException as exc:
                last = ToolFailure(
                    "UPSTREAM_TIMEOUT",
                    f"{host} did not respond within {self.settings.http_timeout:g} s "
                    f"({attempt + 1} attempt(s)). Try again later or raise AFMCP_HTTP_TIMEOUT.",
                    {"url": full, "exception": type(exc).__name__},
                )
            except httpx.TransportError as exc:
                last = ToolFailure(
                    "UPSTREAM_UNAVAILABLE",
                    f"Could not reach {host}: {exc}. Check the network connection (or use replay mode offline).",
                    {"url": full, "exception": type(exc).__name__},
                )
            else:
                if response.status_code == 429 or response.status_code >= 500:
                    code = "RATE_LIMITED" if response.status_code == 429 else "UPSTREAM_UNAVAILABLE"
                    last = ToolFailure(
                        code,
                        f"{host} answered HTTP {response.status_code}; the service is busy or down. Retry later.",
                        {"url": full, "status": response.status_code},
                    )
                else:
                    if use_cache and response.status_code in (200, 404):
                        _CACHE[full] = (now + self.settings.cache_ttl, response.status_code, response.content)
                        while len(_CACHE) > _CACHE_MAX:
                            _CACHE.popitem(last=False)
                    return response.status_code, response.content
            if attempt < self.settings.max_retries:
                await anyio.sleep(self.settings.backoff * (2**attempt))
        assert last is not None
        raise last

    async def get_json(self, url: str, params: dict[str, Any] | None = None) -> tuple[int, Any]:
        status, body = await self.get(url, params)
        try:
            return status, json.loads(body) if body.strip() else None
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise ToolFailure(
                "UPSTREAM_BAD_RESPONSE",
                f"{urlsplit(url).netloc} returned HTTP {status} with a body that is not valid JSON.",
                {"url": url, "status": status, "parse_error": str(exc), "body_start": body[:200].decode("utf-8", "replace")},
            ) from exc


# --------------------------------------------------------------------------------------
# Sequence utilities
# --------------------------------------------------------------------------------------

DNA_ALPHABET = set("ACGTNRYSWKMBDHV")
AMBIGUOUS_DNA = set("NRYSWKMBDHV")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
_COMPLEMENT = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")


def _strip_fasta(text: str, field: str) -> str:
    lines = text.replace("\r", "\n").split("\n")
    headers = [ln for ln in lines if ln.strip().startswith(">")]
    if len(headers) > 1:
        raise ToolFailure(
            "INVALID_INPUT",
            f"{field} contains {len(headers)} FASTA records; send one sequence per call.",
            {"field": field, "records": len(headers)},
        )
    body = "".join(ln for ln in lines if not ln.strip().startswith((">", ";")))
    return re.sub(r"[\s\d]", "", body).upper()


def _bad_characters(seq: str, alphabet: set[str]) -> list[dict[str, Any]]:
    bad = []
    for i, ch in enumerate(seq):
        if ch not in alphabet:
            bad.append({"position": i + 1, "character": ch})
            if len(bad) >= 10:
                break
    return bad


def clean_dna(raw: str, field: str = "dna_sequence") -> tuple[str, QuerySummary]:
    seq = _strip_fasta(raw, field)
    was_rna = "U" in seq
    seq = seq.replace("U", "T")
    bad = _bad_characters(seq, DNA_ALPHABET)
    if bad:
        hint = ""
        if set(seq) <= PROTEIN_ALPHABET:
            hint = " It looks like a protein sequence - use compare_sequence_to_alphafold with query_type='protein'."
        raise ToolFailure(
            "INVALID_INPUT",
            f"{field} contains characters that are not nucleotides (allowed: A C G T/U and IUPAC codes N R Y S W K M B D H V).{hint}",
            {"field": field, "invalid_characters": bad},
        )
    if len(seq) < 30:
        raise ToolFailure(
            "INVALID_INPUT",
            f"{field} has only {len(seq)} nucleotides after removing headers/whitespace; at least 30 are needed.",
            {"field": field, "nucleotide_length": len(seq)},
        )
    gc = sum(seq.count(b) for b in "GCS") / len(seq)
    summary = QuerySummary(
        input_characters=len(raw),
        nucleotide_length=len(seq),
        gc_fraction=round(gc, 4),
        ambiguous_bases=sum(seq.count(b) for b in AMBIGUOUS_DNA),
        was_rna=was_rna,
    )
    return seq, summary


def clean_protein(raw: str, field: str = "query_sequence") -> str:
    seq = _strip_fasta(raw, field).replace("-", "")
    if seq.endswith("*"):
        seq = seq.rstrip("*")
    if "*" in seq:
        raise ToolFailure(
            "INVALID_INPUT",
            f"{field} contains an internal stop ('*') at position {seq.index('*') + 1}; pass a single translated ORF.",
            {"field": field, "position": seq.index("*") + 1},
        )
    bad = _bad_characters(seq, PROTEIN_ALPHABET)
    if bad:
        raise ToolFailure(
            "INVALID_INPUT",
            f"{field} contains characters that are not one-letter amino-acid codes.",
            {"field": field, "invalid_characters": bad},
        )
    if not seq:
        raise ToolFailure("INVALID_INPUT", f"{field} is empty after removing headers/whitespace.", {"field": field})
    if len(seq) >= 30 and set(seq) <= DNA_ALPHABET:
        raise ToolFailure(
            "INVALID_INPUT",
            f"{field} looks like a nucleotide sequence (only A/C/G/T/N letters). Set query_type='dna' to translate it first.",
            {"field": field, "suggested_query_type": "dna"},
        )
    if len(seq) > 10_000:
        raise ToolFailure(
            "INVALID_INPUT", f"{field} is {len(seq)} residues long; the limit is 10000.", {"field": field, "length": len(seq)}
        )
    return seq


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclasses.dataclass
class _Orf:
    frame: Frame
    protein: str
    nt_start: int
    nt_end: int
    has_start: bool
    has_stop: bool

    def candidate(self) -> OrfCandidate:
        return OrfCandidate(
            reading_frame=self.frame,
            strand=self.frame[0],  # type: ignore[arg-type]
            nt_start=self.nt_start,
            nt_end=self.nt_end,
            protein_length=len(self.protein),
            has_start_codon=self.has_start,
            has_stop_codon=self.has_stop,
        )

    def translation(self) -> Translation:
        return Translation(
            **self.candidate().model_dump(),
            protein_sequence=self.protein,
            protein_md5=hashlib.md5(self.protein.encode()).hexdigest(),
            unknown_residues=self.protein.count("X"),
        )


def find_orfs(seq: str, frames: tuple[Frame, ...], min_aa: int) -> tuple[list[_Orf], int]:
    """Return ORFs >= min_aa (longest first) and the length of the longest stop-free stretch seen.

    ATG-initiated ORFs are preferred.  Only if none is long enough do we fall back to
    5'-open reading frames (no start codon), which is what a CDS fragment looks like.
    """
    n = len(seq)
    rc = reverse_complement(seq)
    with_start: list[_Orf] = []
    open_5p: list[_Orf] = []
    longest = 0
    for frame in frames:
        strand, offset = frame[0], int(frame[1]) - 1
        s = seq if strand == "+" else rc
        usable = (len(s) - offset) // 3 * 3
        if usable < 3:
            continue
        protein = str(Seq(s[offset : offset + usable]).translate())
        start = 0
        for i, aa in enumerate(protein + "*"):
            if aa != "*":
                continue
            segment = protein[start:i]
            has_stop = i < len(protein)
            longest = max(longest, len(segment))
            spans: list[tuple[int, bool, list[_Orf]]] = []
            m = segment.find("M")
            if m != -1:
                spans.append((start + m, True, with_start))
            if start == 0:
                spans.append((0, False, open_5p))
            for aa_start, has_start, bucket in spans:
                orf_protein = protein[aa_start:i]
                if len(orf_protein) < min_aa:
                    continue
                s0 = offset + 3 * aa_start + 1  # 1-based on the strand that was translated
                s1 = offset + 3 * i + (3 if has_stop else 0)
                if strand == "-":
                    s0, s1 = n - s1 + 1, n - s0 + 1
                bucket.append(_Orf(frame, orf_protein, s0, s1, has_start, has_stop))
            start = i + 1
    chosen = with_start or open_5p
    order = {f: k for k, f in enumerate(FRAMES)}
    chosen.sort(key=lambda o: (-len(o.protein), order[o.frame], o.nt_start))
    return chosen, longest


def plddt_band(score: float | None) -> PlddtBand | None:
    if score is None:
        return None
    # AFDB convention (reproduces the fractionPlddt* values AFDB publishes): >= 90 / >= 70 / >= 50 / < 50
    if score >= 90:
        return "very_high"
    if score >= 70:
        return "confident"
    if score >= 50:
        return "low"
    return "very_low"


def _segments(scores: list[float], first_residue: int, predicate, min_len: int = 8) -> list[Segment]:
    out, run_start = [], None
    for i, s in enumerate(scores + [None]):  # type: ignore[list-item]
        inside = s is not None and predicate(s)
        if inside and run_start is None:
            run_start = i
        elif not inside and run_start is not None:
            if i - run_start >= min_len:
                vals = scores[run_start:i]
                out.append(
                    Segment(
                        start=first_residue + run_start,
                        end=first_residue + i - 1,
                        length=i - run_start,
                        mean_plddt=round(statistics.fmean(vals), 2),
                    )
                )
            run_start = None
    return out


def summarise_confidence(scores: list[float], first_residue: int, threshold: float) -> ConfidenceSummary:
    n = len(scores)
    bands = [plddt_band(s) for s in scores]
    frac = lambda band: round(bands.count(band) / n, 3)  # noqa: E731
    above = round(sum(s >= threshold for s in scores) / n, 3)
    mean = statistics.fmean(scores)
    confident = _segments(scores, first_residue, lambda s: s >= threshold)
    disordered = _segments(scores, first_residue, lambda s: s < 50)
    if mean >= 90:
        quality = "very high: backbone and most side chains are expected to be accurate"
    elif mean >= 70:
        quality = "good: the backbone is expected to be correct, details less certain"
    elif mean >= 50:
        quality = "low: treat the fold as a rough hypothesis"
    else:
        quality = "very low: most of the chain is probably disordered or unpredictable"
    text = f"Mean pLDDT {mean:.1f} ({quality}). {above:.0%} of residues have pLDDT >= {threshold:g}"
    if confident:
        text += "; confident segment(s): " + ", ".join(f"{s.start}-{s.end}" for s in confident[:6])
    if disordered:
        text += (
            "; pLDDT < 50 segment(s) "
            + ", ".join(f"{s.start}-{s.end}" for s in disordered[:6])
            + " are likely intrinsically disordered (not a reliable 3D structure)"
        )
    return ConfidenceSummary(
        mean_plddt=round(mean, 2),
        median_plddt=round(statistics.median(scores), 2),
        min_plddt=min(scores),
        max_plddt=max(scores),
        fraction_very_high=frac("very_high"),
        fraction_confident=frac("confident"),
        fraction_low=frac("low"),
        fraction_very_low=frac("very_low"),
        threshold=threshold,
        fraction_above_threshold=above,
        confident_segments=confident,
        disordered_segments=disordered,
        interpretation=text + ".",
    )


_BLOSUM62 = substitution_matrices.load("BLOSUM62")
ALIGNMENT_METHOD = "Biopython PairwiseAligner, BLOSUM62, gap open -10 / extend -0.5, free end gaps (semi-global)"


def _blosum(a: str, b: str) -> int | None:
    try:
        return int(_BLOSUM62[a][b])
    except (KeyError, IndexError):
        return None


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = _BLOSUM62
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.end_gap_score = 0.0  # set last: query fragments are not penalised for overhangs
    return aligner


AM_CLASS = {"LBen": "likely_benign", "Amb": "ambiguous", "LPath": "likely_pathogenic"}
MAX_LISTED = 200


def compare_sequences(
    model_seq: str, query_seq: str, scores: list[float], first_residue: int
) -> tuple[Comparison, list[Substitution]]:
    """Align query to model; returns the comparison (without AlphaMissense yet) and its substitution list."""
    # U/O (selenocysteine/pyrrolysine) are not in BLOSUM62; align them as their closest standard residues.
    q_for_alignment = query_seq.replace("U", "C").replace("O", "K")
    alignment = _aligner().align(model_seq, q_for_alignment)[0]
    t_blocks, q_blocks = alignment.aligned
    score_at = lambda i: scores[i] if 0 <= i < len(scores) else None  # noqa: E731

    substitutions: list[Substitution] = []
    indels: list[Indel] = []
    identical = positives = pairs = 0
    aligned_model_idx: list[int] = []
    prev_t = prev_q = None
    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
        ts, te, qs, qe = int(ts), int(te), int(qs), int(qe)
        if prev_t is not None:
            if ts > prev_t:  # model residues absent from the query
                local = [s for s in (score_at(k) for k in range(prev_t, ts)) if s is not None]
                indels.append(
                    Indel(
                        kind="deletion",
                        model_start=first_residue + prev_t,
                        model_end=first_residue + ts - 1,
                        query_start=prev_q,
                        query_end=prev_q,
                        length=ts - prev_t,
                        sequence=model_seq[prev_t:ts],
                        local_mean_plddt=round(statistics.fmean(local), 2) if local else None,
                    )
                )
            if qs > prev_q:  # extra residues in the query
                local = [s for s in (score_at(prev_t - 1), score_at(ts)) if s is not None]
                indels.append(
                    Indel(
                        kind="insertion",
                        model_start=first_residue + prev_t - 1,
                        model_end=first_residue + prev_t - 1,
                        query_start=prev_q + 1,
                        query_end=qs,
                        length=qs - prev_q,
                        sequence=query_seq[prev_q:qs],
                        local_mean_plddt=round(statistics.fmean(local), 2) if local else None,
                    )
                )
        for k in range(te - ts):
            mi, qi = ts + k, qs + k
            a, b = model_seq[mi], query_seq[qi]
            pairs += 1
            aligned_model_idx.append(mi)
            blosum = _blosum(a, q_for_alignment[qi])
            if a == b:
                identical += 1
                positives += 1
                continue
            if blosum is not None and blosum > 0:
                positives += 1
            s = score_at(mi)
            substitutions.append(
                Substitution(
                    variant=f"{a}{first_residue + mi}{b}",
                    model_position=first_residue + mi,
                    query_position=qi + 1,
                    reference_residue=a,
                    query_residue=b,
                    blosum62=blosum,
                    plddt=s,
                    plddt_band=plddt_band(s),
                    alphamissense_pathogenicity=None,
                    alphamissense_class=None,
                )
            )
        prev_t, prev_q = te, qe

    aligned_q = sum(int(e) - int(s) for s, e in q_blocks)
    aligned_t = sum(int(e) - int(s) for s, e in t_blocks)
    identity = identical / pairs if pairs else 0.0
    q_cov = aligned_q / len(query_seq) if query_seq else 0.0
    aligned_scores = [s for s in (score_at(i) for i in aligned_model_idx) if s is not None]
    sub_scores = [s.plddt for s in substitutions if s.plddt is not None]
    if identity >= 0.95 and q_cov >= 0.9:
        transfer = "high"
    elif identity >= 0.40 and q_cov >= 0.7:
        transfer = "medium"
    else:
        transfer = "low"

    n_ins = sum(1 for d in indels if d.kind == "insertion")
    n_del = len(indels) - n_ins
    comparison = Comparison(
        query_source="protein",
        query_length=len(query_seq),
        translation=None,
        alignment_method=ALIGNMENT_METHOD,
        alignment_score=float(alignment.score),
        aligned_pairs=pairs,
        identical=identical,
        identity=round(identity, 4),
        similarity=round(positives / pairs, 4) if pairs else 0.0,
        query_coverage=round(q_cov, 4),
        model_coverage=round(aligned_t / len(model_seq), 4) if model_seq else 0.0,
        n_substitutions=len(substitutions),
        n_insertions=n_ins,
        n_deletions=n_del,
        substitutions=substitutions[:MAX_LISTED],
        indels=indels[:MAX_LISTED],
        lists_truncated=len(substitutions) > MAX_LISTED or len(indels) > MAX_LISTED,
        mean_plddt_aligned_region=round(statistics.fmean(aligned_scores), 2) if aligned_scores else None,
        mean_plddt_at_substitutions=round(statistics.fmean(sub_scores), 2) if sub_scores else None,
        structure_transfer_confidence=transfer,
        interpretation="",
    )
    return comparison, substitutions[:MAX_LISTED]


def describe_comparison(c: Comparison, entry_id: str) -> str:
    parts = [
        f"Query is {c.identity:.1%} identical to {entry_id} over {c.aligned_pairs} aligned residues "
        f"({c.query_coverage:.0%} of the query, {c.model_coverage:.0%} of the model)"
    ]
    if c.n_substitutions == 0 and c.n_insertions == 0 and c.n_deletions == 0:
        parts.append("no substitutions or indels: the AlphaFold model is a direct structural prediction for this sequence")
    else:
        parts.append(f"{c.n_substitutions} substitution(s), {c.n_insertions} insertion(s), {c.n_deletions} deletion(s)")
        for s in c.substitutions[:5]:
            msg = f"{s.variant} at pLDDT {s.plddt:.1f} ({s.plddt_band})" if s.plddt is not None else s.variant
            if s.alphamissense_pathogenicity is not None:
                msg += f", AlphaMissense {s.alphamissense_pathogenicity:.3f} ({s.alphamissense_class})"
            parts.append(msg)
    verdict = {
        "high": "the model can be used as the structure of the query (point changes are not modelled explicitly)",
        "medium": "homology-level match: the overall fold is likely shared, but local structure may differ",
        "low": "weak match: do not transfer the structure without further evidence (try another accession or de novo prediction)",
    }[c.structure_transfer_confidence]
    parts.append(f"Structure-transfer confidence: {c.structure_transfer_confidence} - {verdict}")
    parts.append(
        "Caveat: single-chain AlphaFold models do not capture changes in oligomerisation, ligand binding or aggregation"
    )
    return ". ".join(parts) + "."


# --------------------------------------------------------------------------------------
# External API wrappers
# --------------------------------------------------------------------------------------


def _afdb_url(settings: Settings, url: str | None, fallback_path: str) -> str:
    """Only follow file URLs that point at AFDB (or the configured mirror); otherwise rebuild them."""
    if url and (url.startswith(settings.afdb_base + "/") or url.startswith("https://alphafold.ebi.ac.uk/")):
        if settings.afdb_base != "https://alphafold.ebi.ac.uk":
            url = url.replace("https://alphafold.ebi.ac.uk", settings.afdb_base, 1)
        return url
    return f"{settings.afdb_base}{fallback_path}"


async def afdb_sequence_summary(f: Fetcher, md5: str) -> list[dict[str, Any]] | None:
    status, payload = await f.get_json(f"{f.settings.afdb_base}/api/sequence/summary", {"id": md5, "type": "md5"})
    if status == 404:
        return None
    if status != 200 or not isinstance(payload, dict) or not isinstance(payload.get("structures"), list):
        raise ToolFailure(
            "UPSTREAM_BAD_RESPONSE",
            f"AlphaFold DB sequence search returned HTTP {status} with an unexpected payload.",
            {"status": status, "payload_type": type(payload).__name__},
        )
    return [s.get("summary", {}) for s in payload["structures"] if isinstance(s, dict)]


async def uniprot_annotations(f: Fetcher, accessions: list[str]) -> dict[str, dict[str, Any]]:
    if not accessions:
        return {}
    status, payload = await f.get_json(
        f"{f.settings.uniprot_base}/uniprotkb/accessions",
        {
            "accessions": ",".join(accessions),
            "fields": "accession,reviewed,organism_name,organism_id,gene_primary,protein_name,length",
            "format": "json",
        },
    )
    if status != 200 or not isinstance(payload, dict):
        raise ToolFailure("UPSTREAM_BAD_RESPONSE", f"UniProt returned HTTP {status}.", {"status": status})
    out: dict[str, dict[str, Any]] = {}
    for rec in payload.get("results", []):
        acc = rec.get("primaryAccession")
        if not acc:
            continue
        desc = rec.get("proteinDescription", {})
        name = (desc.get("recommendedName") or {}).get("fullName", {}).get("value")
        if not name and desc.get("submissionNames"):
            name = desc["submissionNames"][0].get("fullName", {}).get("value")
        genes = rec.get("genes") or [{}]
        org = rec.get("organism", {})
        out[acc] = {
            "protein_name": name,
            "gene": (genes[0].get("geneName") or {}).get("value"),
            "organism": org.get("scientificName"),
            "taxon_id": org.get("taxonId"),
            "reviewed": "reviewed" in rec.get("entryType", "").lower() and "unreviewed" not in rec.get("entryType", "").lower(),
        }
    return out


async def afdb_prediction(f: Fetcher, accession: str) -> tuple[dict[str, Any], list[str]]:
    status, payload = await f.get_json(f"{f.settings.afdb_base}/api/prediction/{accession}")
    if status in (404,) or (status == 200 and not payload):
        raise ToolFailure(
            "NOT_FOUND",
            f"AlphaFold DB has no model for UniProt accession {accession}. Check the accession, or use "
            "find_alphafold_models_for_dna to search by sequence.",
            {"uniprot_accession": accession, "status": status},
        )
    if status == 400:
        message = payload.get("error") if isinstance(payload, dict) else None
        raise ToolFailure(
            "INVALID_INPUT",
            f"AlphaFold DB rejected the identifier {accession}: {message or 'invalid identifier'}",
            {"uniprot_accession": accession, "status": status},
        )
    if status != 200 or not isinstance(payload, list):
        raise ToolFailure(
            "UPSTREAM_BAD_RESPONSE",
            f"AlphaFold DB prediction endpoint returned HTTP {status} with an unexpected payload.",
            {"status": status},
        )
    records = [r for r in payload if isinstance(r, dict)]
    exact = [r for r in records if r.get("uniprotAccession") == accession]

    def rank(r: dict[str, Any]) -> tuple:
        return (
            r.get("providerId") != "GDM",
            not str(r.get("entryId", "")).endswith("-F1"),
            bool(r.get("isComplex")),
            -(r.get("latestVersion") or 0),
        )

    chosen = sorted(exact or records, key=rank)[0]
    others = [str(r.get("entryId")) for r in records if r is not chosen]
    return chosen, others


async def afdb_plddt(f: Fetcher, record: dict[str, Any]) -> list[float]:
    entry, version = record.get("entryId"), record.get("latestVersion")
    url = _afdb_url(f.settings, record.get("plddtDocUrl"), f"/files/{entry}-confidence_v{version}.json")
    status, payload = await f.get_json(url)
    scores = payload.get("confidenceScore") if isinstance(payload, dict) else None
    if status != 200 or not isinstance(scores, list) or not scores:
        raise ToolFailure(
            "UPSTREAM_BAD_RESPONSE",
            f"Could not read per-residue pLDDT for {entry} (HTTP {status}).",
            {"url": url, "status": status},
        )
    try:
        return [float(s) for s in scores]
    except (TypeError, ValueError) as exc:
        raise ToolFailure("UPSTREAM_BAD_RESPONSE", f"Non-numeric pLDDT values for {entry}.", {"url": url}) from exc


async def alphamissense(f: Fetcher, url: str, wanted: set[str]) -> dict[str, tuple[float, str]]:
    status, body = await f.get(url, accept="text/csv")
    if status != 200:
        raise ToolFailure("UPSTREAM_BAD_RESPONSE", f"AlphaMissense table returned HTTP {status}.", {"url": url})
    out: dict[str, tuple[float, str]] = {}
    for line in body.decode("utf-8", "replace").splitlines():
        variant, _, rest = line.partition(",")
        if variant in wanted:
            score, _, klass = rest.partition(",")
            out[variant] = (float(score), AM_CLASS.get(klass.strip(), "ambiguous"))
    return out


# --------------------------------------------------------------------------------------
# Envelope helpers
# --------------------------------------------------------------------------------------


def _validate(model: type[BaseModel], raw: dict[str, Any]) -> Any:
    """Validate tool arguments; every problem becomes one INVALID_INPUT error listing all bad fields."""
    try:
        return model.model_validate(raw)
    except ValidationError as exc:
        raise _validation_failure(exc) from exc


FRIENDLY_PROBLEMS = {
    ("uniprot_accession", "string_pattern_mismatch"): (
        "must be a UniProtKB accession such as P68871, P69905 or P04637-2 (not a gene or entry name like HBB_HUMAN)"
    ),
}


def _validation_failure(exc: ValidationError) -> ToolFailure:
    problems = []
    for err in exc.errors():
        field = ".".join(str(p) for p in err["loc"]) or "(arguments)"
        problem = FRIENDLY_PROBLEMS.get((field, err["type"]), err["msg"])
        problems.append({"field": field, "problem": problem, "type": err["type"]})
    missing = [p["field"] for p in problems if p["type"] == "missing"]
    summary = "; ".join(f"{p['field']}: {p['problem']}" for p in problems[:5])
    message = f"Invalid arguments - {summary}."
    if missing:
        message = f"Missing required argument(s): {', '.join(missing)}. " + message
    return ToolFailure("INVALID_INPUT", message, {"problems": problems})


def _result(
    response_model: type[_Envelope],
    tool: str,
    started: float,
    fetcher: Fetcher | None,
    warnings: list[str],
    data: BaseModel | None = None,
    failure: ToolFailure | None = None,
) -> CallToolResult:
    provenance = Provenance(
        sources=fetcher.sources if fetcher else [],
        retrieved_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        elapsed_ms=int((time.perf_counter() - started) * 1000),
        replay_mode=bool(fetcher.replay) if fetcher else bool(SETTINGS.replay_dir),
    )
    error = None
    if failure is not None:
        error = ErrorInfo(
            code=failure.code,  # type: ignore[arg-type]
            message=failure.message,
            retryable=RETRYABLE.get(failure.code, False),
            details=failure.details,
        )
    envelope = response_model(
        ok=failure is None, tool=tool, data=data, error=error, warnings=warnings, provenance=provenance  # type: ignore[call-arg]
    )
    structured = envelope.model_dump(mode="json")
    text = json.dumps(structured, indent=2)
    return CallToolResult(content=[TextContent(type="text", text=text)], structured_content=structured, is_error=failure is not None)


async def _run_tool(name: str, response_model: type[_Envelope], body) -> CallToolResult:
    """Run a tool body and convert *every* outcome into the structured envelope."""
    started = time.perf_counter()
    warnings: list[str] = []
    fetcher: Fetcher | None = None
    try:
        fetcher = Fetcher(SETTINGS)
        async with fetcher:
            data = await body(fetcher, warnings)
        return _result(response_model, name, started, fetcher, warnings, data=data)
    except ToolFailure as failure:
        log.info("%s failed: %s %s", name, failure.code, failure.message)
        return _result(response_model, name, started, fetcher, warnings, failure=failure)
    except Exception as exc:  # pragma: no cover - last line of defence, keeps the session alive
        log.exception("%s crashed", name)
        failure = ToolFailure("INTERNAL_ERROR", f"Unexpected server error ({type(exc).__name__}): {exc}", {})
        return _result(response_model, name, started, fetcher, warnings, failure=failure)


def _collect(**kwargs: Any) -> dict[str, Any]:
    """Drop arguments the client did not send (or sent as null) so the documented defaults apply."""
    return {k: v for k, v in kwargs.items() if v is not None}

# --------------------------------------------------------------------------------------
# Tool 1: DNA -> ORF -> exact AlphaFold DB models
# --------------------------------------------------------------------------------------

FIND_DESCRIPTION = """\
Translate a coding DNA/RNA sequence and find AlphaFold DB structure models whose protein sequence is IDENTICAL to it.

Use this first when you have a nucleotide sequence (CDS, cDNA or mRNA) and want to know which protein it encodes and \
whether AlphaFold already predicted its structure. The longest ATG-initiated ORF (or the frame you choose) is \
translated, its MD5 is looked up in AlphaFold DB, and each hit is annotated with UniProt organism, gene and \
Swiss-Prot/TrEMBL status (ranked: requested taxon, reviewed, monomer, mean pLDDT).

Returns: translation (frame, coordinates, protein sequence), match_status 'exact_match' or 'no_exact_match', \
ranked matches and best_match. 'no_exact_match' is NOT an error: it usually means the DNA carries a variant - call \
compare_sequence_to_alphafold with the wild-type UniProt accession and the translated protein to locate it.
Errors (ok=false): INVALID_INPUT (non-nucleotide characters, < 30 nt, bad argument), NO_ORF_FOUND, \
UPSTREAM_TIMEOUT / UPSTREAM_UNAVAILABLE / RATE_LIMITED (retryable), UPSTREAM_BAD_RESPONSE.
"""


async def find_alphafold_models_for_dna(
    dna_sequence: Any = None,
    reading_frame: Any = None,
    min_protein_length: Any = None,
    taxon_id: Any = None,
    include_complexes: Any = None,
    max_results: Any = None,
) -> Annotated[CallToolResult, FindModelsResponse]:
    raw = _collect(
        dna_sequence=dna_sequence,
        reading_frame=reading_frame,
        min_protein_length=min_protein_length,
        taxon_id=taxon_id,
        include_complexes=include_complexes,
        max_results=max_results,
    )

    async def body(f: Fetcher, warnings: list[str]) -> FindModelsData:
        args = _validate(FindModelsInput, raw)
        seq, query = clean_dna(args.dna_sequence)
        if query.ambiguous_bases > 0.05 * query.nucleotide_length:
            warnings.append(
                f"{query.ambiguous_bases} ambiguous bases; codons containing them become X and will prevent exact matches."
            )
        frames = FRAMES if args.reading_frame == "auto" else (args.reading_frame,)
        orfs, longest = find_orfs(seq, frames, args.min_protein_length)
        if not orfs:
            raise ToolFailure(
                "NO_ORF_FOUND",
                f"No open reading frame of >= {args.min_protein_length} codons in frame(s) {', '.join(frames)}; "
                f"the longest stop-free stretch is {longest} codons. Check that this is coding sequence (not intronic or "
                "UTR), try another reading_frame, or lower min_protein_length.",
                {"frames_searched": list(frames), "longest_stop_free_codons": longest,
                 "min_protein_length": args.min_protein_length},
            )
        best = orfs[0]
        translation = best.translation()
        alternatives = [o.candidate() for o in orfs[1:] if o.frame != best.frame][:3]
        if not best.has_start:
            warnings.append("ORF has no ATG start codon (treated as a 5'-truncated CDS fragment).")
        if not best.has_stop:
            warnings.append("ORF runs to the end of the sequence without a stop codon (3'-truncated?).")
        if translation.unknown_residues:
            warnings.append(f"Translated protein contains {translation.unknown_residues} X residue(s).")

        summaries = await afdb_sequence_summary(f, translation.protein_md5)
        matches: list[ModelMatch] = []
        if summaries:
            seen: set[str] = set()
            for s in summaries:
                model_id = str(s.get("model_identifier", ""))
                state = s.get("oligomeric_state")
                if not model_id or model_id in seen:
                    continue
                if not args.include_complexes and state not in (None, "MONOMER"):
                    continue
                seen.add(model_id)
                ents = s.get("entities") or [{}]
                acc = next((e.get("identifier") for e in ents if e.get("identifier_category") == "UNIPROT"), None)
                matches.append(
                    ModelMatch(
                        model_id=model_id,
                        uniprot_accession=acc,
                        protein_name=ents[0].get("description"),
                        gene=None,
                        organism=None,
                        taxon_id=None,
                        reviewed=None,
                        taxon_match=None,
                        oligomeric_state=state,
                        mean_plddt=s.get("confidence_avg_local_score"),
                        sequence_identity=s.get("sequence_identity"),
                        coverage=s.get("coverage"),
                        model_url=s.get("model_url"),
                        model_page_url=s.get("model_page_url"),
                        created=s.get("created"),
                    )
                )
            accessions = list(dict.fromkeys(m.uniprot_accession for m in matches if m.uniprot_accession))
            if len(accessions) > 100:
                warnings.append(f"{len(accessions)} UniProt entries share this sequence; only the first 100 were annotated.")
            try:
                notes = await uniprot_annotations(f, accessions[:100])
            except ToolFailure as exc:
                notes = {}
                warnings.append(f"UniProt annotation unavailable ({exc.code}); organisms/genes are missing from matches.")
            for m in matches:
                info = notes.get(m.uniprot_accession or "", {})
                m.protein_name = info.get("protein_name") or m.protein_name
                m.gene = info.get("gene")
                m.organism = info.get("organism")
                m.taxon_id = info.get("taxon_id")
                m.reviewed = info.get("reviewed")
                if args.taxon_id is not None and m.taxon_id is not None:
                    m.taxon_match = m.taxon_id == args.taxon_id
            canonical = lambda m: m.model_id == f"AF-{m.uniprot_accession}-F1"  # noqa: E731
            matches.sort(
                key=lambda m: (
                    not bool(m.taxon_match),
                    not bool(m.reviewed),
                    m.oligomeric_state not in (None, "MONOMER"),
                    not canonical(m),
                    -(m.mean_plddt or 0.0),
                )
            )
            if args.taxon_id is not None and not any(m.taxon_match for m in matches):
                warnings.append(f"No identical model from taxon {args.taxon_id}; showing other organisms.")

        total = len(matches)
        matches = matches[: args.max_results]
        top = matches[0] if matches else None
        if top:
            next_step = (
                f"Call compare_sequence_to_alphafold with uniprot_accession='{top.uniprot_accession}' to get the "
                "pLDDT confidence profile and model files for this protein."
            )
        else:
            next_step = (
                "No identical AlphaFold model. If you know the gene, find its wild-type UniProt accession and call "
                "compare_sequence_to_alphafold(uniprot_accession=..., query_sequence=<translation.protein_sequence>) "
                "to see which residues differ and how confident the model is there."
            )
        return FindModelsData(
            query=query,
            translation=translation,
            alternative_orfs=alternatives,
            match_status="exact_match" if matches else "no_exact_match",
            total_matches=total,
            matches=matches,
            best_match=top,
            next_step=next_step,
        )

    return await _run_tool("find_alphafold_models_for_dna", FindModelsResponse, body)


# --------------------------------------------------------------------------------------
# Tool 2: compare a sequence with the AlphaFold model of a UniProt accession
# --------------------------------------------------------------------------------------

COMPARE_DESCRIPTION = """\
Compare a protein (or coding DNA) with the AlphaFold DB model of a UniProt accession and report how far the \
predicted structure can be trusted for that sequence.

Always returns the model (entry id, version, PDB/mmCIF/PAE URLs) and its confidence profile: mean/median pLDDT, \
fraction of residues per AFDB confidence band, confident segments and likely-disordered segments (pLDDT < 50). \
If query_sequence is given, it is aligned to the model sequence (semi-global, BLOSUM62) and the result lists \
identity, coverage, every substitution with the model's pLDDT at that residue and - for human proteins - the \
AlphaMissense pathogenicity score, plus insertions/deletions and a structure-transfer verdict (high/medium/low).

Typical use: after find_alphafold_models_for_dna, or with the wild-type accession and a mutant sequence \
(e.g. P68871 + sickle-cell HBB) to see where the variant sits in the structure.
Errors (ok=false): INVALID_INPUT (malformed accession, non-amino-acid characters, DNA passed as protein), \
NO_ORF_FOUND (query_type='dna'), NOT_FOUND (no AFDB model), UPSTREAM_TIMEOUT / UPSTREAM_UNAVAILABLE / \
RATE_LIMITED (retryable), UPSTREAM_BAD_RESPONSE.
"""


async def compare_sequence_to_alphafold(
    uniprot_accession: Any = None,
    query_sequence: Any = None,
    query_type: Any = None,
    plddt_threshold: Any = None,
    include_per_residue_plddt: Any = None,
    include_alphamissense: Any = None,
) -> Annotated[CallToolResult, CompareResponse]:
    raw = _collect(
        uniprot_accession=uniprot_accession,
        query_sequence=query_sequence,
        query_type=query_type,
        plddt_threshold=plddt_threshold,
        include_per_residue_plddt=include_per_residue_plddt,
        include_alphamissense=include_alphamissense,
    )
    if isinstance(raw.get("uniprot_accession"), str):  # be forgiving about case/whitespace
        raw["uniprot_accession"] = raw["uniprot_accession"].strip().upper()

    async def body(f: Fetcher, warnings: list[str]) -> CompareData:
        args = _validate(CompareInput, raw)
        # Validate/parse the query before touching the network so bad input fails fast.
        query_protein: str | None = None
        query_translation: Translation | None = None
        if args.query_sequence is not None and args.query_sequence.strip():
            if args.query_type == "dna":
                seq, _ = clean_dna(args.query_sequence, "query_sequence")
                orfs, longest = find_orfs(seq, FRAMES, 10)
                if not orfs:
                    raise ToolFailure(
                        "NO_ORF_FOUND",
                        f"query_sequence (DNA) has no ORF of >= 10 codons (longest stop-free stretch: {longest}).",
                        {"longest_stop_free_codons": longest},
                    )
                query_translation = orfs[0].translation()
                query_protein = query_translation.protein_sequence
            else:
                query_protein = clean_protein(args.query_sequence)

        record, others = await afdb_prediction(f, args.uniprot_accession)
        model_seq = str(record.get("uniprotSequence") or record.get("sequence") or "")
        if not model_seq:
            raise ToolFailure("UPSTREAM_BAD_RESPONSE", "AlphaFold DB record has no sequence.", {"entry": record.get("entryId")})
        scores = await afdb_plddt(f, record)
        if len(scores) != len(model_seq):
            warnings.append(
                f"pLDDT array has {len(scores)} values for a {len(model_seq)}-residue sequence; positions beyond the "
                "shorter of the two have no confidence value."
            )
        first = int(record.get("uniprotStart") or record.get("sequenceStart") or 1)
        entry = str(record.get("entryId"))
        version = record.get("latestVersion")
        am_url = record.get("amAnnotationsUrl")
        model = ModelInfo(
            entry_id=entry,
            uniprot_accession=str(record.get("uniprotAccession") or args.uniprot_accession),
            uniprot_id=record.get("uniprotId"),
            gene=record.get("gene"),
            protein_name=record.get("uniprotDescription"),
            organism=record.get("organismScientificName"),
            taxon_id=record.get("taxId"),
            reviewed=record.get("isReviewed", record.get("isUniProtReviewed")),
            model_version=version,
            model_created=record.get("modelCreatedDate"),
            tool_used=record.get("toolUsed"),
            provider=record.get("providerId"),
            sequence_length=len(model_seq),
            model_start=first,
            model_end=first + len(model_seq) - 1,
            afdb_reported_mean_plddt=record.get("globalMetricValue"),
            urls=ModelUrls(
                entry_page=f"https://alphafold.ebi.ac.uk/entry/{record.get('uniprotAccession') or args.uniprot_accession}",
                pdb=record.get("pdbUrl"),
                cif=record.get("cifUrl"),
                plddt_json=record.get("plddtDocUrl"),
                pae_json=record.get("paeDocUrl"),
                pae_image=record.get("paeImageUrl"),
                alphamissense_csv=am_url,
            ),
            other_models_for_accession=others,
        )
        if others:
            warnings.append(f"{len(others)} other model(s) exist for this accession (isoforms/complexes); using {entry}.")
        confidence = summarise_confidence(scores, first, args.plddt_threshold)

        comparison: Comparison | None = None
        if query_protein is not None:
            comparison, subs = compare_sequences(model_seq, query_protein, scores, first)
            if query_translation is not None:
                comparison.query_source = "dna"
                comparison.translation = query_translation
            if subs and args.include_alphamissense:
                if am_url:
                    try:
                        table = await alphamissense(
                            f, _afdb_url(f.settings, am_url, f"/files/{entry}-aa-substitutions.csv"), {s.variant for s in subs}
                        )
                        for s in comparison.substitutions:
                            if s.variant in table:
                                s.alphamissense_pathogenicity, s.alphamissense_class = table[s.variant]  # type: ignore[assignment]
                    except ToolFailure as exc:
                        warnings.append(f"AlphaMissense scores unavailable ({exc.code}).")
                else:
                    warnings.append("AlphaMissense scores are only published for human canonical UniProt sequences.")
            comparison.interpretation = describe_comparison(comparison, entry)

        if comparison is None:
            next_step = (
                f"Download the model from {model.urls.pdb} (or mmCIF) and colour it by B-factor (= pLDDT). "
                "Pass query_sequence to compare your own protein/DNA with this model."
            )
        elif comparison.structure_transfer_confidence == "low":
            next_step = "The query does not match this model well; search for a closer homolog or predict it de novo (AlphaFold/ColabFold)."
        else:
            next_step = (
                f"Inspect the listed positions in {model.urls.pdb}; for a mutant-specific structure, run AlphaFold/ColabFold "
                "on the query sequence and superimpose it on this model."
            )
        return CompareData(
            model=model,
            confidence=confidence,
            per_residue_plddt=scores if args.include_per_residue_plddt else None,
            comparison=comparison,
            next_step=next_step,
        )

    return await _run_tool("compare_sequence_to_alphafold", CompareResponse, body)


# --------------------------------------------------------------------------------------
# Server assembly
# --------------------------------------------------------------------------------------


def _schema(model: type[BaseModel]) -> dict[str, Any]:
    schema = model.model_json_schema()
    schema.pop("description", None)
    return schema


def _build_tool(fn, input_model: type[BaseModel], title: str, description: str) -> Tool:
    tool = Tool.from_function(
        fn,
        title=title,
        description=description,
        annotations=ToolAnnotations(
            title=title, read_only_hint=True, destructive_hint=False, idempotent_hint=True, open_world_hint=True
        ),
    )
    # The function accepts loosely-typed arguments so that *we* validate them and can always answer with the
    # structured error envelope; clients still see (and LLMs follow) the strict JSON Schema below.
    tool.parameters = _schema(input_model)
    return tool


INSTRUCTIONS = """\
Tools for going from a nucleotide sequence to an AlphaFold-predicted protein structure.
Workflow: (1) find_alphafold_models_for_dna(dna_sequence) translates the ORF and finds AlphaFold DB models with an \
identical sequence; (2) compare_sequence_to_alphafold(uniprot_accession, query_sequence) reports the model's pLDDT \
confidence and, for a variant sequence, where each difference sits and how confident/pathogenic it is.
Every result is a JSON envelope: check 'ok'; on failure read error.code / error.message (retryable errors can be retried).
pLDDT >= 90 very high, 70-90 confident, 50-70 low, < 50 very low (often disordered)."""

TOOLS = [
    _build_tool(
        find_alphafold_models_for_dna,
        FindModelsInput,
        "Find AlphaFold models for a DNA sequence",
        FIND_DESCRIPTION,
    ),
    _build_tool(
        compare_sequence_to_alphafold,
        CompareInput,
        "Compare a sequence with an AlphaFold model",
        COMPARE_DESCRIPTION,
    ),
]

mcp = MCPServer(
    name=SERVER_NAME,
    title="AlphaFold sequence-to-structure tools",
    version=SERVER_VERSION,
    instructions=INSTRUCTIONS,
    tools=TOOLS,
    log_level=os.environ.get("AFMCP_LOG_LEVEL", "INFO").upper(),  # type: ignore[arg-type]
)


def tool_schemas() -> list[dict[str, Any]]:
    """The tool definitions exactly as advertised over MCP (used by --print-schemas and the design doc)."""
    return [
        {
            "name": t.name,
            "title": t.title,
            "inputSchema": t.parameters,
            "outputSchema": t.output_schema,
        }
        for t in TOOLS
    ]


def main(argv: list[str] | None = None) -> None:
    global SETTINGS
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--transport",
        choices=["stdio", "streamable-http"],
        default=os.environ.get("AFMCP_TRANSPORT", "stdio"),
        help="stdio (default; for Claude Desktop/Code, MCP Inspector) or streamable-http",
    )
    parser.add_argument("--host", default="127.0.0.1", help="bind address for streamable-http (default 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000, help="port for streamable-http (default 8000)")
    parser.add_argument(
        "--replay",
        metavar="DIR",
        help="offline mode: answer from recorded API responses in DIR (e.g. ./fixtures) instead of the network",
    )
    parser.add_argument("--print-schemas", action="store_true", help="print the advertised tool schemas as JSON and exit")
    args = parser.parse_args(argv)

    logging.basicConfig(
        stream=sys.stderr,  # stdout is the MCP channel in stdio mode
        level=os.environ.get("AFMCP_LOG_LEVEL", "INFO").upper(),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    if args.replay:
        SETTINGS = dataclasses.replace(SETTINGS, replay_dir=str(Path(args.replay).resolve()))
    if SETTINGS.replay_dir and not (Path(SETTINGS.replay_dir) / "manifest.json").is_file():
        parser.error(f"replay directory {SETTINGS.replay_dir!r} has no manifest.json (copy the fixtures/ folder too)")
    if SETTINGS.replay_dir:
        log.warning("REPLAY MODE: answering from recorded responses in %s (no network calls)", SETTINGS.replay_dir)
    if args.print_schemas:
        print(json.dumps(tool_schemas(), indent=2))
        return
    if args.transport == "stdio":
        mcp.run("stdio")
    else:
        log.info("Streamable HTTP endpoint: http://%s:%d/mcp", args.host, args.port)
        mcp.run("streamable-http", host=args.host, port=args.port)


if __name__ == "__main__":
    main()
