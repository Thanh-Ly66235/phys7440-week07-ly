"""Test suite for mcp_server_ly.py (PHYS 7440, Assignment 7) - Thanh Ly.

Every test talks to the server through a real MCP client (``mcp.Client``), so it exercises tool
listing, argument passing, the advertised JSON Schemas and the structured result exactly as an
LLM client would.  Each result is validated against the tool's advertised ``outputSchema``.

Two groups:

* offline (default, deterministic): HTTP calls are answered from recorded AlphaFold DB / UniProt
  responses in ./fixtures (the server's replay mode), or from a mock transport for failure cases;
* live (marker ``live``): the same scenarios against the real APIs.  They are skipped automatically
  when the APIs are unreachable, or on request with ``-m "not live"`` / ``AFMCP_SKIP_LIVE=1``.

Run:  pytest -v test_mcp_server_ly.py
"""

import dataclasses
import os
import sys
from pathlib import Path

import httpx
import pytest
from jsonschema import Draft202012Validator
from mcp import Client, StdioServerParameters

import mcp_server_ly as server

HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"
EXAMPLES = HERE / "examples"
FIND = "find_alphafold_models_for_dna"
COMPARE = "compare_sequence_to_alphafold"

HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLV"
    "CVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
)
# Random 60-codon ORF (seed 7440) - its protein has no exact AlphaFold DB match.
SYNTHETIC_ORF = (
    "ATGACGGGTCCGGTTTATCACAAAGACGTATTCAATACCCGAAGGCGACTACAGTTCACCTCAAATCCGTCCCCAACTTGTCATTACACGCATACCACCCGTGCGAATC"
    "AGACGGTAATCACCGATTCTCCCCGGACCAAAGCCACAGACCGCGCCCAAATGAAGTTAGCTTCCGATAGGTGA"
)

pytestmark = pytest.mark.anyio


@pytest.fixture
def anyio_backend():
    return "asyncio"


def example(name: str) -> str:
    return (EXAMPLES / name).read_text()


# ----------------------------------------------------------------------------- helpers


async def call(client: Client, tool: str, arguments: dict) -> dict:
    """Call a tool and check the result against the advertised outputSchema; return structured content."""
    listed = {t.name: t for t in (await client.list_tools()).tools}
    result = await client.call_tool(tool, arguments)
    content = result.structured_content
    assert content is not None, "tool must return structuredContent"
    errors = sorted(Draft202012Validator(listed[tool].output_schema).iter_errors(content), key=str)
    assert not errors, f"output does not match outputSchema: {errors[0]}"
    # isError flag, 'ok' and error/data fields must agree
    assert result.is_error is (not content["ok"])
    assert (content["error"] is None) == content["ok"]
    assert (content["data"] is not None) == content["ok"]
    assert result.content and result.content[0].type == "text"  # text fallback for non-structured clients
    return content


def input_errors(tool, arguments: dict) -> list[str]:
    return [e.message for e in Draft202012Validator(tool.input_schema).iter_errors(arguments)]


def assert_error(content: dict, code: str, retryable: bool = False) -> dict:
    assert content["ok"] is False
    assert content["error"]["code"] == code, content["error"]
    assert content["error"]["retryable"] is retryable
    assert content["error"]["message"]
    return content["error"]


@pytest.fixture
def offline(monkeypatch):
    """Answer every HTTP request from ./fixtures; no retries/backoff so failures are fast."""
    settings = dataclasses.replace(server.SETTINGS, replay_dir=str(FIXTURES), max_retries=0, backoff=0.0, cache_ttl=0)
    monkeypatch.setattr(server, "SETTINGS", settings)
    monkeypatch.setattr(server, "TRANSPORT_OVERRIDE", None)
    return settings


@pytest.fixture
def mock_http(monkeypatch):
    """Install an httpx.MockTransport; returns a list that records every request made."""
    seen: list[httpx.Request] = []

    def install(handler, retries: int = 0):
        def wrapped(request: httpx.Request):
            seen.append(request)
            return handler(request)

        monkeypatch.setattr(server, "TRANSPORT_OVERRIDE", httpx.MockTransport(wrapped))
        monkeypatch.setattr(
            server, "SETTINGS", dataclasses.replace(server.SETTINGS, replay_dir="", max_retries=retries, backoff=0.0, cache_ttl=0)
        )
        return seen

    return install


# ----------------------------------------------------------------------------- discovery


async def test_server_advertises_two_tools_with_schemas(offline):
    async with Client(server.mcp) as client:
        tools = {t.name: t for t in (await client.list_tools()).tools}
    assert set(tools) == {FIND, COMPARE}
    for tool in tools.values():
        assert tool.input_schema["type"] == "object" and tool.input_schema["additionalProperties"] is False
        assert tool.output_schema["type"] == "object"
        assert {"ok", "tool", "data", "error", "warnings", "provenance"} <= set(tool.output_schema["properties"])
        assert tool.annotations.read_only_hint is True
        Draft202012Validator.check_schema(tool.input_schema)
        Draft202012Validator.check_schema(tool.output_schema)
    assert tools[FIND].input_schema["required"] == ["dna_sequence"]
    assert tools[COMPARE].input_schema["required"] == ["uniprot_accession"]


# ----------------------------------------------------------------------------- tool 1: valid inputs


async def test_find_hbb_cds_returns_human_hemoglobin_first(offline):
    """VALID: wild-type HBB CDS -> exact AlphaFold match, human Swiss-Prot entry ranked first."""
    args = {"dna_sequence": example("hbb_cds.fasta"), "taxon_id": 9606}
    async with Client(server.mcp) as client:
        tool = {t.name: t for t in (await client.list_tools()).tools}[FIND]
        assert input_errors(tool, args) == []
        out = await call(client, FIND, args)
    data = out["data"]
    t = data["translation"]
    assert (t["reading_frame"], t["nt_start"], t["nt_end"]) == ("+1", 1, 444)
    assert t["protein_sequence"] == HBB_PROTEIN and t["protein_length"] == 147
    assert t["has_start_codon"] and t["has_stop_codon"]
    assert t["protein_md5"] == "209d686939d0b8d1ea089368150140e2"
    assert data["match_status"] == "exact_match"
    best = data["best_match"]
    assert best["model_id"] == "AF-P68871-F1" and best["uniprot_accession"] == "P68871"
    assert best["organism"] == "Homo sapiens" and best["reviewed"] is True and best["taxon_match"] is True
    # complexes are filtered out by default
    assert all(m["oligomeric_state"] == "MONOMER" for m in data["matches"])
    assert "P68871" in data["next_step"]
    assert out["provenance"]["replay_mode"] is True
    assert any("/api/sequence/summary" in u for u in out["provenance"]["sources"])


async def test_find_sickle_allele_matches_variant_entry_and_complexes_optional(offline):
    """VALID: the sickle-cell CDS (GAG>GTG) translates to a different protein with its own AFDB model."""
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": example("hbb_sickle_E7V_cds.fasta")})
        both = await call(client, FIND, {"dna_sequence": example("hbb_cds.fasta"), "include_complexes": True})
    assert out["data"]["translation"]["protein_sequence"][:9] == "MVHLTPVEK"
    assert out["data"]["best_match"]["uniprot_accession"] == "A0A481SHK9"
    assert out["data"]["best_match"]["reviewed"] is False  # TrEMBL entry
    states = {m["oligomeric_state"] for m in both["data"]["matches"]}
    assert "HOMODIMER" in states


async def test_find_accepts_rna_reverse_strand_lowercase(offline):
    """VALID: reverse-complemented, lower-case RNA with a FASTA header is still translated correctly."""
    dna = "".join(l for l in example("hbb_cds.fasta").splitlines() if not l.startswith(">"))
    rna = ">reverse strand\n" + server.reverse_complement(dna).replace("T", "U").lower()
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": rna, "max_results": 1})
    data = out["data"]
    assert data["query"]["was_rna"] is True
    assert data["translation"]["reading_frame"] == "-1" and data["translation"]["strand"] == "-"
    assert (data["translation"]["nt_start"], data["translation"]["nt_end"]) == (1, 444)
    assert data["translation"]["protein_sequence"] == HBB_PROTEIN
    assert len(data["matches"]) == 1 and data["total_matches"] >= 1


async def test_find_without_identical_model_is_success_with_hint(offline):
    """VALID: a novel ORF is not an error - match_status says so and next_step explains what to do."""
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": SYNTHETIC_ORF})
    data = out["data"]
    assert out["ok"] is True
    assert data["match_status"] == "no_exact_match" and data["matches"] == [] and data["best_match"] is None
    assert data["translation"]["protein_length"] == 60
    assert "compare_sequence_to_alphafold" in data["next_step"]


async def test_find_tp53_wild_type_vs_hotspot_mutant(offline):
    """VALID: wild-type TP53 CDS finds the human p53 model; the R175H allele has no identical model."""
    async with Client(server.mcp) as client:
        wt = await call(client, FIND, {"dna_sequence": example("tp53_cds.fasta"), "taxon_id": 9606})
        mut = await call(client, FIND, {"dna_sequence": example("tp53_R175H_cds.fasta"), "taxon_id": 9606})
    ranked = [(m["uniprot_accession"], m["reviewed"], m["taxon_match"]) for m in wt["data"]["matches"]]
    assert ranked[:2] == [("P04637", True, True), ("K7PPA8", False, True)]  # human Swiss-Prot, then human TrEMBL
    assert wt["data"]["translation"]["protein_length"] == 393
    assert mut["data"]["match_status"] == "no_exact_match"
    assert mut["data"]["translation"]["protein_sequence"][174] == "H"


# ----------------------------------------------------------------------------- tool 1: invalid inputs


async def test_find_rejects_non_nucleotide_characters(offline):
    """ERROR: a protein sequence passed as DNA -> INVALID_INPUT naming the offending characters."""
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": HBB_PROTEIN})
    err = assert_error(out, "INVALID_INPUT")
    assert err["details"]["field"] == "dna_sequence"
    # M, V and H are legal IUPAC nucleotide codes, so the first illegal letter is L at position 4
    assert err["details"]["invalid_characters"][0] == {"position": 4, "character": "L"}
    assert "protein" in err["message"]
    assert out["provenance"]["sources"] == []  # rejected before any network call


async def test_find_missing_and_out_of_range_arguments(offline):
    """ERROR: missing required parameter, and several schema violations reported together."""
    async with Client(server.mcp) as client:
        tool = {t.name: t for t in (await client.list_tools()).tools}[FIND]
        missing = await call(client, FIND, {"reading_frame": "+1"})
        bad_args = {"dna_sequence": "ACGT" * 20, "reading_frame": "+4", "max_results": 500}
        assert len(input_errors(tool, bad_args)) == 2  # the advertised schema already forbids these
        bad = await call(client, FIND, bad_args)
    err = assert_error(missing, "INVALID_INPUT")
    assert "Missing required argument(s): dna_sequence" in err["message"]
    err = assert_error(bad, "INVALID_INPUT")
    assert {p["field"] for p in err["details"]["problems"]} == {"reading_frame", "max_results"}


async def test_find_reports_no_orf(offline):
    """ERROR: a sequence with stop codons in all six frames -> NO_ORF_FOUND (not a crash)."""
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": "TTAG" * 25})
    err = assert_error(out, "NO_ORF_FOUND")
    assert err["details"]["longest_stop_free_codons"] < 30
    assert len(err["details"]["frames_searched"]) == 6


# ----------------------------------------------------------------------------- tool 2: valid inputs


async def test_compare_confidence_profile_matches_afdb_summary(offline):
    """VALID: p53 model without a query -> confidence bands identical to AFDB's published fractions."""
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P04637", "include_per_residue_plddt": True})
    data = out["data"]
    assert data["model"]["entry_id"] == "AF-P04637-F1"  # canonical model chosen over isoform models
    assert "AF-P04637-2-F1" in data["model"]["other_models_for_accession"]
    c = data["confidence"]
    assert (c["fraction_very_low"], c["fraction_low"], c["fraction_confident"], c["fraction_very_high"]) == (
        0.298, 0.104, 0.071, 0.527)
    assert abs(c["mean_plddt"] - data["model"]["afdb_reported_mean_plddt"]) < 0.05
    assert any(s["start"] <= 175 <= s["end"] for s in c["confident_segments"])  # DNA-binding domain
    assert any(s["start"] == 1 for s in c["disordered_segments"])  # disordered N-terminal TAD
    assert len(data["per_residue_plddt"]) == 393
    assert data["comparison"] is None


async def test_compare_sickle_protein_to_wild_type_model(offline):
    """VALID: sickle-cell beta-globin vs AF-P68871-F1 -> one substitution E7V with pLDDT + AlphaMissense."""
    sickle_protein = HBB_PROTEIN.replace("MVHLTPEEK", "MVHLTPVEK", 1)
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P68871", "query_sequence": sickle_protein})
    comp = out["data"]["comparison"]
    assert comp["aligned_pairs"] == 147 and comp["identical"] == 146
    assert comp["query_coverage"] == 1.0 and comp["structure_transfer_confidence"] == "high"
    [sub] = comp["substitutions"]
    assert sub["variant"] == "E7V" and sub["plddt"] == 96.56 and sub["plddt_band"] == "very_high"
    assert sub["alphamissense_pathogenicity"] == 0.2229 and sub["alphamissense_class"] == "likely_benign"
    assert "E7V" in comp["interpretation"]


async def test_compare_dna_query_tp53_r175h(offline):
    """VALID: query given as DNA (TP53 c.524G>A) is translated and the hotspot R175H is located."""
    async with Client(server.mcp) as client:
        out = await call(
            client, COMPARE, {"uniprot_accession": "p04637 ", "query_sequence": example("tp53_R175H_cds.fasta"), "query_type": "dna"}
        )
    comp = out["data"]["comparison"]
    assert comp["query_source"] == "dna" and comp["translation"]["protein_length"] == 393
    [sub] = comp["substitutions"]
    assert sub["variant"] == "R175H" and sub["plddt_band"] == "very_high"
    assert sub["alphamissense_class"] == "likely_pathogenic" and sub["alphamissense_pathogenicity"] > 0.9


async def test_compare_fragment_with_deletion(offline):
    """VALID: a partial query with an internal 3-residue deletion -> coverage and indel are reported."""
    fragment = HBB_PROTEIN[20:60] + HBB_PROTEIN[63:120]  # residues 21-60 and 64-120 (61-63 deleted)
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P68871", "query_sequence": fragment})
    comp = out["data"]["comparison"]
    assert comp["identity"] == 1.0 and comp["query_coverage"] == 1.0
    assert comp["model_coverage"] == round(97 / 147, 4)
    assert comp["n_deletions"] == 1 and comp["n_substitutions"] == 0
    [indel] = comp["indels"]
    assert (indel["kind"], indel["model_start"], indel["model_end"], indel["sequence"]) == ("deletion", 61, 63, HBB_PROTEIN[60:63])


# ----------------------------------------------------------------------------- tool 2: invalid inputs / failures


async def test_compare_rejects_malformed_accession_without_network(offline):
    """ERROR: invalid accession format -> INVALID_INPUT before any HTTP request."""
    async with Client(server.mcp) as client:
        tool = {t.name: t for t in (await client.list_tools()).tools}[COMPARE]
        assert input_errors(tool, {"uniprot_accession": "HBB_HUMAN!"})
        out = await call(client, COMPARE, {"uniprot_accession": "HBB_HUMAN!"})
    err = assert_error(out, "INVALID_INPUT")
    assert err["details"]["problems"][0]["field"] == "uniprot_accession"
    assert "UniProtKB accession such as P68871" in err["message"]
    assert out["provenance"]["sources"] == []


async def test_compare_unknown_accession_is_not_found(offline):
    """ERROR: well-formed accession without an AlphaFold model -> NOT_FOUND (HTTP 404 upstream)."""
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "O00000"})
    err = assert_error(out, "NOT_FOUND")
    assert err["details"]["status"] == 404


async def test_compare_dna_sent_as_protein(offline):
    """ERROR: nucleotide text with query_type='protein' -> INVALID_INPUT suggesting query_type='dna'."""
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P68871", "query_sequence": example("hbb_cds.fasta")})
    err = assert_error(out, "INVALID_INPUT")
    assert err["details"]["suggested_query_type"] == "dna"


async def test_compare_upstream_timeout_is_structured_and_retried(mock_http):
    """ERROR: the API never answers -> UPSTREAM_TIMEOUT (retryable) after 1 + max_retries attempts."""

    def handler(request):
        raise httpx.ReadTimeout("simulated timeout", request=request)

    seen = mock_http(handler, retries=2)
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P68871"})
    err = assert_error(out, "UPSTREAM_TIMEOUT", retryable=True)
    assert len(seen) == 3
    assert err["details"]["url"].endswith("/api/prediction/P68871")


async def test_find_upstream_server_error_and_bad_json(mock_http):
    """ERROR: HTTP 503 -> UPSTREAM_UNAVAILABLE; HTML instead of JSON -> UPSTREAM_BAD_RESPONSE."""
    mock_http(lambda request: httpx.Response(503, text="Service Unavailable"))
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": example("hbb_cds.fasta")})
    assert_error(out, "UPSTREAM_UNAVAILABLE", retryable=True)

    mock_http(lambda request: httpx.Response(200, text="<html>maintenance</html>"))
    async with Client(server.mcp) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P68871"})
    assert_error(out, "UPSTREAM_BAD_RESPONSE")


async def test_uniprot_outage_degrades_to_warning(mock_http):
    """Partial failure: AFDB works but UniProt is down -> success with a warning instead of an error."""
    replay = server.ReplayTransport(FIXTURES)

    async def handler(request):
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(500)
        return await replay.handle_async_request(request)

    mock_http(handler)
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": example("hbb_cds.fasta")})
    assert out["ok"] is True and out["data"]["match_status"] == "exact_match"
    assert any("UniProt annotation unavailable" in w for w in out["warnings"])


# ----------------------------------------------------------------------------- real stdio transport


async def test_stdio_subprocess_end_to_end():
    """The script itself starts as an MCP stdio server (replay mode) and answers a tool call."""
    params = StdioServerParameters(
        command=sys.executable,
        args=[str(HERE / "mcp_server_ly.py"), "--replay", str(FIXTURES)],
        env={**os.environ, "AFMCP_LOG_LEVEL": "WARNING"},
        cwd=str(HERE),
    )
    async with Client(params) as client:
        out = await call(client, COMPARE, {"uniprot_accession": "P68871", "query_sequence": HBB_PROTEIN})
    assert out["ok"] and out["data"]["comparison"]["identity"] == 1.0


# ----------------------------------------------------------------------------- live API tests


@pytest.fixture
def live(monkeypatch):
    if os.environ.get("AFMCP_SKIP_LIVE") == "1":
        pytest.skip("AFMCP_SKIP_LIVE=1")
    try:
        httpx.get("https://alphafold.ebi.ac.uk/api/prediction/P68871", timeout=10).raise_for_status()
        httpx.get("https://rest.uniprot.org/uniprotkb/P68871.json", timeout=10).raise_for_status()
    except httpx.HTTPError as exc:
        pytest.skip(f"AlphaFold DB / UniProt not reachable: {exc}")
    monkeypatch.setattr(server, "TRANSPORT_OVERRIDE", None)
    monkeypatch.setattr(server, "SETTINGS", dataclasses.replace(server.SETTINGS, replay_dir=""))


@pytest.mark.live
async def test_live_find_hbb(live):
    async with Client(server.mcp) as client:
        out = await call(client, FIND, {"dna_sequence": example("hbb_cds.fasta"), "taxon_id": 9606})
    assert out["ok"], out["error"]
    assert out["provenance"]["replay_mode"] is False
    assert out["data"]["best_match"]["uniprot_accession"] == "P68871"
    assert out["data"]["best_match"]["organism"] == "Homo sapiens"


@pytest.mark.live
async def test_live_compare_sickle_dna(live):
    async with Client(server.mcp) as client:
        out = await call(
            client, COMPARE,
            {"uniprot_accession": "P68871", "query_sequence": example("hbb_sickle_E7V_cds.fasta"), "query_type": "dna"},
        )
    assert out["ok"], out["error"]
    [sub] = out["data"]["comparison"]["substitutions"]
    assert sub["variant"] == "E7V" and sub["plddt"] > 90
    assert 0 <= sub["alphamissense_pathogenicity"] <= 1
    assert out["data"]["confidence"]["mean_plddt"] > 90


@pytest.mark.live
async def test_live_errors(live):
    async with Client(server.mcp) as client:
        missing = await call(client, COMPARE, {"uniprot_accession": "O00000"})
        bad = await call(client, FIND, {"dna_sequence": "NOT DNA AT ALL, JUST WORDS AND MORE WORDS"})
    assert_error(missing, "NOT_FOUND")
    assert_error(bad, "INVALID_INPUT")
