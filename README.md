# AlphaFold Sequence-to-Structure MCP Server

**Thanh Ly**, PHYS 7440 GenAI in Physics Research, Assignment 7 (Project MCP Server Prototype)

`mcp_server_ly.py` is a [Model Context Protocol](https://modelcontextprotocol.io/specification) server. It lets an LLM client go from a **coding DNA sequence** to an **AlphaFold-predicted protein structure**, and it reports how far that prediction can be trusted at the residues where the sequence differs from the reference. It wraps two public data sources:

* the **AlphaFold Protein Structure Database (AFDB)** REST API: models, per-residue pLDDT and AlphaMissense tables;
* the **UniProt REST API**: organism, gene and Swiss-Prot/TrEMBL status of the matching entries.

Protein folding is a statistical-physics problem: the native structure sits at the free-energy minimum, and AlphaFold's pLDDT is a learned estimate of how well that minimum is pinned down for each residue. My project question is: *given a gene variant, does the AlphaFold structure of the wild-type protein tell me anything reliable about the variant, and where does the variant sit in the confidence landscape?*

| Tool | Question it answers | External API calls |
|---|---|---|
| `find_alphafold_models_for_dna` | "Which protein does this DNA encode, and does AlphaFold already have a structure for exactly that sequence?" | AFDB `GET /api/sequence/summary?id=<md5>&type=md5`, UniProt `GET /uniprotkb/accessions` |
| `compare_sequence_to_alphafold` | "How reliable is the AlphaFold model of protein X, and where do the differences in my sequence fall (pLDDT, AlphaMissense)?" | AFDB `GET /api/prediction/{accession}`, `GET /files/AF-…-confidence_v6.json`, `GET /files/AF-…-aa-substitutions.csv` |

Neither API needs an account or a key.

---

## 1. Folder contents

| File | Purpose |
|---|---|
| `mcp_server_ly.py` | **The MCP server** (Python `mcp` SDK v2, `MCPServer`). Defines both tools, their schemas and the error handling. |
| `test_mcp_server_ly.py` | **pytest suite**: 20 offline tests (valid and invalid inputs for both tools, schema validation, stdio transport) plus 3 live-API tests. |
| `design_ly.md` / `design_ly.pdf` | **One-page design document**: the endpoints, their inputSchemas, the output format and the data sources. |
| `README.md` | This file. |
| `requirements.txt` | Pinned dependency ranges. |
| `pytest.ini` | Test configuration (defines the `live` marker). |
| `demo_client_ly.py` | Command-line MCP client for the class demo. It runs the whole DNA-to-structure pipeline over stdio. |
| `examples/` | Real RefSeq CDS FASTA files (HBB wild type, sickle-cell HBB E7V, TP53 wild type, TP53 R175H) and a Claude Desktop config snippet. |
| `example_outputs/` | Tool results (arguments, `isError`, `structuredContent`) produced by `demo_client_ly.py --replay fixtures --save`, i.e. real tool calls on API data recorded on 26 Oct 2026. |
| `fixtures/` | Recorded AFDB/UniProt responses (captured 26 Oct 2026). Used by the offline tests and by `--replay` mode. |

---

## 2. Installation (fresh Python environment)

You need Python **3.10 or newer** (tested with 3.13) and internet access to reach `alphafold.ebi.ac.uk` and `rest.uniprot.org`. Prebuilt wheels exist for all dependencies on Linux, macOS and Windows, so nothing has to be compiled.

```bash
cd ~/Downloads/"Assignment 7"          # the folder that contains this README
python3 -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
```

> **Copying the project to another machine** (for example a JupyterHub server): copy the *whole* folder, including `examples/`, `fixtures/` and `pytest.ini`. Jupyter's upload button only takes files, so upload `Assignment 7.zip` instead and extract it in a terminal with `python -m zipfile -e "Assignment 7.zip" .` (this works even where `unzip` is not installed). If anything is missing, `pytest` stops straight away and lists the missing files.

Check the installation:

```bash
pytest -v                              # 20 offline tests pass; 3 live tests run if the APIs are reachable
python mcp_server_ly.py --print-schemas | head -40    # the tool schemas exactly as advertised
python demo_client_ly.py               # full pipeline against the live APIs
```

> The server targets **`mcp` 2.x**. In that version `FastMCP` was renamed to `mcp.server.mcpserver.MCPServer`. A plain `pip install mcp` installs 2.x; `requirements.txt` pins it exactly (`mcp==2.3.0`), the version all tests were run with, so a later SDK release cannot change the behaviour.

---

## 3. Running the server

| Command | Use |
|---|---|
| `python mcp_server_ly.py` | **stdio** transport (default). An MCP client starts this process itself. |
| `python mcp_server_ly.py --transport streamable-http --port 8000` | HTTP server at `http://127.0.0.1:8000/mcp`, for remote or web clients. |
| `python mcp_server_ly.py --replay fixtures` | **Offline replay mode.** Answers come from recorded API responses. This is a backup for the in-class demo if Wi-Fi fails. It only knows the recorded sequences (HBB, HBB E7V, TP53, TP53 R175H). Any other request returns `UPSTREAM_UNAVAILABLE`. |
| `python mcp_server_ly.py --print-schemas` | Prints the advertised `inputSchema`/`outputSchema` of both tools and exits. |

Logs go to **stderr**, because stdout carries the MCP protocol in stdio mode.

### Environment variables (all optional)

| Variable | Default | Meaning |
|---|---|---|
| `AFMCP_HTTP_TIMEOUT` | `20` | Timeout in seconds for each HTTP request to AFDB/UniProt. |
| `AFMCP_MAX_RETRIES` | `2` | Extra attempts after a timeout, connection error, HTTP 5xx or HTTP 429. |
| `AFMCP_BACKOFF` | `0.75` | First retry delay in seconds; it doubles on each retry. |
| `AFMCP_CACHE_TTL` | `3600` | In-memory cache lifetime in seconds for successful and 404 responses. `0` disables the cache. |
| `AFMCP_CONTACT` | *(empty)* | E-mail address appended to the User-Agent. This is optional API etiquette. |
| `AFMCP_AFDB_BASE` | `https://alphafold.ebi.ac.uk` | AFDB base URL, for a mirror. |
| `AFMCP_UNIPROT_BASE` | `https://rest.uniprot.org` | UniProt REST base URL. |
| `AFMCP_REPLAY_DIR` | *(empty)* | Same as `--replay DIR`. |
| `AFMCP_LOG_LEVEL` | `INFO` | `DEBUG`, `INFO`, `WARNING` or `ERROR`. |
| `AFMCP_TRANSPORT` | `stdio` | Default for `--transport`. |
| `AFMCP_SKIP_LIVE` | *(unset)* | When set to `1`, the test suite skips the live-API tests. |

No API keys or secrets are required.

---

## 4. Connecting an LLM client (class demo, Meeting 2)

In the commands below, replace the folder path if your copy lives somewhere else, and keep the quotes because the folder name contains a space.

**Claude Desktop**: open *Settings → Developer → Edit Config* and add the block from `examples/claude_desktop_config.example.json`:

```json
{
  "mcpServers": {
    "afdb-sequence": {
      "command": "/home/thanhly/Downloads/Assignment 7/.venv/bin/python",
      "args": ["/home/thanhly/Downloads/Assignment 7/mcp_server_ly.py"]
    }
  }
}
```

Restart Claude Desktop. The two tools then appear under the connectors/tools menu.

**Claude Code**:

```bash
claude mcp add afdb-sequence -- "$PWD/.venv/bin/python" "$PWD/mcp_server_ly.py"
```

**MCP Inspector** (a web UI for testing calls by hand; needs Node.js):

```bash
npx @modelcontextprotocol/inspector "$PWD/.venv/bin/python" "$PWD/mcp_server_ly.py"
```

**Example prompt for the live demo** (paste `examples/hbb_sickle_E7V_cds.fasta`):

> Here is a coding sequence from a patient sample. Which protein does it encode, is there an AlphaFold structure for it, and how does it differ from the wild-type human protein? Tell me how confident the structure is at the changed residue.

The model should call `find_alphafold_models_for_dna`, which returns an exact match to the TrEMBL sickle-cell entry A0A481SHK9. It should then call `compare_sequence_to_alphafold` with `P68871`, which reports E7V at pLDDT 96.6 with AlphaMissense 0.22 (likely benign).

Verified before submission: Claude Code (headless `claude -p`, stdio transport, replay mode) discovered both tools, called `find_alphafold_models_for_dna` and then `compare_sequence_to_alphafold`, and answered with E7V, pLDDT 96.56 and AlphaMissense 0.223.

**Without an LLM** (backup):

```bash
python demo_client_ly.py                                    # live APIs
python demo_client_ly.py --replay fixtures --errors         # offline, also shows the error handling
python demo_client_ly.py --dna examples/tp53_R175H_cds.fasta --reference P04637
```

---

## 5. Tool reference

Both tools are read-only, idempotent and open-world (`ToolAnnotations`). Each returns `structuredContent` that matches its advertised `outputSchema`, plus the same JSON as a text block for clients that ignore structured output.

### Common result envelope

```jsonc
{
  "ok": true,                      // false => "error" is set, "data" is null, MCP isError = true
  "tool": "compare_sequence_to_alphafold",
  "data": { ... },                 // tool-specific, see below
  "error": null,                   // {code, message, retryable, details}
  "warnings": [],                  // non-fatal issues (e.g. UniProt down, ORF without stop codon)
  "provenance": {"server": "afdb-sequence-mcp", "server_version": "1.0.0",
                 "sources": ["https://alphafold.ebi.ac.uk/api/prediction/P68871", "..."],
                 "retrieved_at": "2026-10-26T18:03:07+00:00", "elapsed_ms": 412, "replay_mode": false}   // timing values are illustrative
}
```

### 5.1 `find_alphafold_models_for_dna`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `dna_sequence` | string, 30–100 000 chars | **required** | CDS/cDNA/mRNA as plain text or one FASTA record. Case, whitespace and digits are ignored, U is read as T, IUPAC codes are allowed. |
| `reading_frame` | `"auto"`, `"+1"`…`"-3"` | `"auto"` | `auto` picks the longest ATG-initiated ORF over all 6 frames. A 5′-truncated ORF is used only if no ATG ORF is long enough. |
| `min_protein_length` | integer 10–5000 | 30 | Minimum ORF length in amino acids. |
| `taxon_id` | integer or null | null | Hits from this NCBI taxon are ranked first (9606 = human). |
| `include_complexes` | boolean | false | Also list AFDB homo- and hetero-oligomer models. |
| `max_results` | integer 1–50 | 10 | |

Ranking of matches: requested taxon, then Swiss-Prot before TrEMBL, then monomer, then canonical `AF-<acc>-F1` model, then higher mean pLDDT.

**Example call**: wild-type HBB CDS (`examples/hbb_cds.fasta`, NM_000518.5), `{"dna_sequence": "<FASTA text>", "taxon_id": 9606, "max_results": 3}`.
**Expected output** (abridged; produced in replay mode from API data recorded on 26 Oct 2026, full file `example_outputs/find_hbb_cds.json` shows the same call without `max_results`):

```json
{
  "ok": true,
  "data": {
    "query": {"nucleotide_length": 444, "gc_fraction": 0.5608, "ambiguous_bases": 0, "was_rna": false},
    "translation": {"reading_frame": "+1", "strand": "+", "nt_start": 1, "nt_end": 444, "protein_length": 147,
                    "has_start_codon": true, "has_stop_codon": true,
                    "protein_sequence": "MVHLTPEEKSAVTALWGKVN...", "protein_md5": "209d686939d0b8d1ea089368150140e2"},
    "alternative_orfs": [{"reading_frame": "-1", "nt_start": 1, "nt_end": 351, "protein_length": 117, "...": "..."}],
    "match_status": "exact_match",
    "total_matches": 3,
    "matches": [
      {"model_id": "AF-P68871-F1", "uniprot_accession": "P68871", "protein_name": "Hemoglobin subunit beta",
       "gene": "HBB", "organism": "Homo sapiens", "taxon_id": 9606, "reviewed": true, "taxon_match": true,
       "oligomeric_state": "MONOMER", "mean_plddt": 97.19,
       "model_url": "https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v6.cif"},
      {"model_id": "AF-P68873-F1", "organism": "Pan troglodytes", "reviewed": true, "taxon_match": false, "mean_plddt": 96.94, "...": "..."},
      {"model_id": "AF-P68872-F1", "organism": "Pan paniscus", "reviewed": true, "taxon_match": false, "mean_plddt": 96.88, "...": "..."}
    ],
    "best_match": {"model_id": "AF-P68871-F1", "...": "..."},
    "next_step": "Call compare_sequence_to_alphafold with uniprot_accession='P68871' to get the pLDDT confidence profile and model files for this protein."
  }
}
```

Human, chimpanzee and bonobo β-globin have identical sequences, so all three share the same MD5 and the human entry is ranked first. With `examples/tp53_R175H_cds.fasta` the result is `"match_status": "no_exact_match"`. That is **not** an error: no AFDB model has exactly that mutant sequence, and `next_step` tells the agent to compare against the wild type instead.

### 5.2 `compare_sequence_to_alphafold`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `uniprot_accession` | string matching the UniProt accession regex | **required** | e.g. `P68871`, `P04637-2`. Case and surrounding whitespace are tolerated. |
| `query_sequence` | string ≤ 100 000 chars or null | null | Protein, or DNA if `query_type="dna"`. Omit it to get only the confidence profile. |
| `query_type` | `"protein"` or `"dna"` | `"protein"` | `dna` translates the longest ORF first. |
| `plddt_threshold` | number 0–100 | 70 | Threshold for "confident segments". |
| `include_per_residue_plddt` | boolean | false | Return the full pLDDT array. |
| `include_alphamissense` | boolean | true | Look up AlphaMissense for each substitution (human canonical proteins only). |

The pLDDT bands follow the AFDB convention: very high ≥ 90, confident 70–90, low 50–70, very low < 50. The fractions this server computes reproduce the `fractionPlddt*` values AFDB publishes (checked for P68871 and P04637). Alignment is semi-global (free end gaps) with BLOSUM62, gap open −10 and extend −0.5, so fragments are handled. `structure_transfer_confidence` is `high` when identity ≥ 95 % and query coverage ≥ 90 %, `medium` when identity ≥ 40 % and coverage ≥ 70 %, and `low` otherwise.

**Example call**: sickle-cell allele vs the wild-type model,
`{"uniprot_accession": "P68871", "query_sequence": "<examples/hbb_sickle_E7V_cds.fasta>", "query_type": "dna"}`.
**Expected output** (abridged; full file: `example_outputs/compare_hbb_sickle_E7V_cds_vs_P68871.json`):

```json
{
  "ok": true,
  "data": {
    "model": {"entry_id": "AF-P68871-F1", "uniprot_id": "HBB_HUMAN", "organism": "Homo sapiens", "model_version": 6,
              "sequence_length": 147, "afdb_reported_mean_plddt": 97.19,
              "urls": {"pdb": "https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v6.pdb", "...": "..."}},
    "confidence": {"mean_plddt": 97.18, "median_plddt": 98.12, "fraction_very_high": 0.973, "fraction_confident": 0.02,
                   "fraction_low": 0.007, "fraction_very_low": 0.0,
                   "confident_segments": [{"start": 2, "end": 147, "length": 146, "mean_plddt": 97.42}],
                   "disordered_segments": [],
                   "interpretation": "Mean pLDDT 97.2 (very high: ...). 99% of residues have pLDDT >= 70; confident segment(s): 2-147."},
    "comparison": {
      "query_source": "dna", "aligned_pairs": 147, "identical": 146, "identity": 0.9932,
      "query_coverage": 1.0, "model_coverage": 1.0, "n_substitutions": 1, "n_insertions": 0, "n_deletions": 0,
      "substitutions": [{"variant": "E7V", "model_position": 7, "reference_residue": "E", "query_residue": "V",
                         "blosum62": -2, "plddt": 96.56, "plddt_band": "very_high",
                         "alphamissense_pathogenicity": 0.2229, "alphamissense_class": "likely_benign"}],
      "structure_transfer_confidence": "high",
      "interpretation": "Query is 99.3% identical to AF-P68871-F1 ... E7V at pLDDT 96.6 (very_high), AlphaMissense 0.223 (likely_benign). ..."
    }
  }
}
```

For TP53 R175H vs `P04637`, the substitution is `R175H` at pLDDT 96.62 in the DNA-binding domain, with AlphaMissense 0.9857 (likely pathogenic). The confidence profile flags residues 1–11, 29–65, 295–307, 309–321 and 358–383 as likely disordered (pLDDT < 50).

**Physics note on the HbS example.** AlphaMissense calls HbS (E7V, written "E6V" in mature-chain numbering) *likely benign*, even though it causes sickle-cell disease. The single-chain fold is unchanged. The disease comes from a new hydrophobic patch that drives intermolecular polymerisation of deoxy-haemoglobin, an aggregation and phase-transition effect that single-chain structure and variant-effect predictors do not model. The tool reports this caveat in its `interpretation` field. It illustrates the limits my project has to keep in mind.

---

## 6. Error handling

Each tool validates its own arguments with the same pydantic models that generate the advertised schema. Wrong types, out-of-range values and missing parameters therefore come back as structured `INVALID_INPUT` results rather than SDK exceptions. No exception escapes a tool call, so the MCP session never crashes.

| `error.code` | Retryable | When | Example message |
|---|---|---|---|
| `INVALID_INPUT` | no | Missing or invalid argument, non-nucleotide characters, < 30 nt, malformed accession, DNA passed as protein, internal stop in a protein | `Missing required argument(s): dna_sequence. ...` |
| `NO_ORF_FOUND` | no | No ORF ≥ `min_protein_length` codons in the searched frames | `No open reading frame of >= 30 codons in frame(s) +1, ... the longest stop-free stretch is 3 codons ...` |
| `NOT_FOUND` | no | AFDB has no model for the accession (HTTP 404) | `AlphaFold DB has no model for UniProt accession O00000 ...` |
| `UPSTREAM_TIMEOUT` | yes | No answer within `AFMCP_HTTP_TIMEOUT` after all retries | `alphafold.ebi.ac.uk did not respond within 20 s (3 attempt(s)) ...` |
| `UPSTREAM_UNAVAILABLE` | yes | Connection failure or HTTP 5xx after retries; also an unrecorded URL in replay mode | `Could not reach rest.uniprot.org: ...` |
| `RATE_LIMITED` | yes | HTTP 429 after retries | `... answered HTTP 429; the service is busy ...` |
| `UPSTREAM_BAD_RESPONSE` | no | Non-JSON or unexpected payload, e.g. after an API change | `... returned HTTP 200 with a body that is not valid JSON.` |
| `INTERNAL_ERROR` | no | Unexpected bug; it is logged with a traceback on stderr | `Unexpected server error (...)` |

Example error result (`example_outputs/error_bad_accession.json`):

```json
{
  "ok": false,
  "tool": "compare_sequence_to_alphafold",
  "data": null,
  "error": {
    "code": "INVALID_INPUT",
    "message": "Invalid arguments - uniprot_accession: must be a UniProtKB accession such as P68871, P69905 or P04637-2 (not a gene or entry name like HBB_HUMAN).",
    "retryable": false,
    "details": {"problems": [{"field": "uniprot_accession", "problem": "must be a UniProtKB accession such as ...", "type": "string_pattern_mismatch"}]}
  },
  "warnings": [],
  "provenance": {"sources": [], "...": "..."}
}
```

Some failures are partial and degrade to a warning instead of an error. If UniProt is down, tool 1 still returns its AFDB matches, without organism data. If AlphaMissense data is missing, for example for non-human proteins, tool 2 still returns the comparison with `null` scores.

---

## 7. Tests

```bash
pytest -v                       # everything (live tests auto-skip if the APIs are unreachable)
pytest -v -m "not live"         # offline only: deterministic, no network needed
pytest -v -m live               # only the real-API tests
```

All calls go through a real MCP client (`mcp.Client`). Every result is validated against the tool's advertised `outputSchema` with `jsonschema` (Draft 2020-12), and the tests check that `isError`, `ok`, `error` and `data` are consistent with each other.

| Test | Tool | Kind |
|---|---|---|
| `test_server_advertises_two_tools_with_schemas` | both | discovery: schemas, annotations |
| `test_find_hbb_cds_returns_human_hemoglobin_first` | 1 | **valid** |
| `test_find_sickle_allele_matches_variant_entry_and_complexes_optional` | 1 | **valid** |
| `test_find_accepts_rna_reverse_strand_lowercase` | 1 | **valid** (frame −1, RNA, FASTA) |
| `test_find_tp53_wild_type_vs_hotspot_mutant` | 1 | **valid** (ranking, no exact match) |
| `test_find_without_identical_model_is_success_with_hint` | 1 | **valid** |
| `test_find_rejects_non_nucleotide_characters` | 1 | **error** `INVALID_INPUT` |
| `test_find_missing_and_out_of_range_arguments` | 1 | **error** `INVALID_INPUT` (missing parameter, enum/range) |
| `test_find_reports_no_orf` | 1 | **error** `NO_ORF_FOUND` |
| `test_compare_confidence_profile_matches_afdb_summary` | 2 | **valid** (bands = AFDB values) |
| `test_compare_sickle_protein_to_wild_type_model` | 2 | **valid** |
| `test_compare_dna_query_tp53_r175h` | 2 | **valid** (DNA query) |
| `test_compare_fragment_with_deletion` | 2 | **valid** (coverage, indel) |
| `test_compare_rejects_malformed_accession_without_network` | 2 | **error** `INVALID_INPUT` |
| `test_compare_unknown_accession_is_not_found` | 2 | **error** `NOT_FOUND` |
| `test_compare_dna_sent_as_protein` | 2 | **error** `INVALID_INPUT` |
| `test_compare_upstream_timeout_is_structured_and_retried` | 2 | **error** `UPSTREAM_TIMEOUT` (simulated API timeout, 3 attempts) |
| `test_find_upstream_server_error_and_bad_json` | 1, 2 | **error** `UPSTREAM_UNAVAILABLE`, `UPSTREAM_BAD_RESPONSE` |
| `test_uniprot_outage_degrades_to_warning` | 1 | partial failure becomes a warning |
| `test_stdio_subprocess_end_to_end` | 2 | real stdio transport (subprocess) |
| `test_live_find_hbb`, `test_live_compare_sickle_dna`, `test_live_errors` | 1, 2 | live AFDB/UniProt (marker `live`) |

The offline tests use the recorded responses in `fixtures/` (see `fixtures/manifest.json`), so they are reproducible and do not depend on the EBI/UniProt servers.

---

## 8. Known failure modes and limitations

1. **Matching is exact, not by homology.** Tool 1 matches the protein's MD5, so a single amino-acid change gives `no_exact_match`. The workaround is to call tool 2 with the wild-type accession. A BLAST or Foldseek homology search is planned for later weeks.
2. **Genomic DNA is not spliced.** Introns, UTR-only or non-coding input either produce `NO_ORF_FOUND` or a short spurious ORF. Pass the CDS or mRNA, and check `alternative_orfs` or force `reading_frame` when in doubt. Only the standard genetic code is supported, so mitochondrial or other alternative codes will mistranslate.
3. **AFDB coverage.** Some proteins have no model: very long proteins, many viral polyproteins, sequences outside UniProt. Isoform models exist only for some accessions. The canonical `AF-<acc>-F1` model is chosen and the others are listed in `other_models_for_accession`.
4. **AlphaMissense** scores exist only for human canonical sequences, and they use UniProt numbering (HbS is E7V here, "E6V" in the clinical literature). Elsewhere the scores are `null` and a warning is added.
5. **What pLDDT measures.** pLDDT is the model's *local* confidence. It does not predict the effect of a mutation on stability, binding or oligomerisation, as the HbS example shows. Domain packing (PAE) is linked (`urls.pae_json`) but not analysed yet.
6. **Upstream availability.** EBI and UniProt outages, rate limits or slow responses give `UPSTREAM_*` / `RATE_LIMITED` errors after 2 retries. Raise `AFMCP_HTTP_TIMEOUT` on slow networks. An API schema change, such as AFDB moving from v6 to v7 file names, would surface as `UPSTREAM_BAD_RESPONSE` or `NOT_FOUND`.
7. **Caching.** Responses are cached in memory for up to 1 h (`AFMCP_CACHE_TTL`), so a newly released AFDB version may take that long to appear. The cache is lost when the server restarts.
8. **Replay mode** only covers the four example sequences and the error cases in `fixtures/`. Every other request fails with `UPSTREAM_UNAVAILABLE` by design, so replay mode never invents data.
9. **Input size.** The limits are 100 000 nt of DNA and 10 000 aa of protein query. Alignment cost scales as length², so very long queries take a few seconds.
10. **SDK version.** `mcp` 1.x will not import this module, because `FastMCP` was renamed in v2. Use the pinned requirements.

---

## 9. Data sources and references

* AlphaFold DB API: <https://alphafold.ebi.ac.uk/api-docs>. Varadi M. *et al.*, "AlphaFold Protein Structure Database in 2024", *Nucleic Acids Res.* 52, D368 (2024).
* AlphaFold: Jumper J. *et al.*, "Highly accurate protein structure prediction with AlphaFold", *Nature* 596, 583 (2021).
* AlphaMissense: Cheng J. *et al.*, "Accurate proteome-wide missense variant effect prediction with AlphaMissense", *Science* 381, eadg7492 (2023).
* UniProt REST API: <https://rest.uniprot.org>. The UniProt Consortium, *Nucleic Acids Res.* (2025).
* Example sequences: NCBI RefSeq NM_000518.5 (HBB) and NM_000546.6 (TP53), coding sequences retrieved via NCBI E-utilities.
* MCP specification: <https://modelcontextprotocol.io/specification>. Python SDK: <https://github.com/modelcontextprotocol/python-sdk>.

AlphaFold DB data are available under CC-BY 4.0, and UniProt data under CC BY 4.0.
