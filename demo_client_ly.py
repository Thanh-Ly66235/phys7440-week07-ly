"""Command-line MCP client for the class demo (PHYS 7440, Assignment 7) - Thanh Ly.

Starts mcp_server_ly.py as a stdio MCP server (exactly how Claude Desktop / Claude Code launch it),
lists its tools, and runs the DNA -> AlphaFold pipeline:

    1. find_alphafold_models_for_dna   (translate the DNA, find identical AlphaFold models)
    2. compare_sequence_to_alphafold   (confidence profile + variant positions vs a reference model)

Examples
    python demo_client_ly.py                                   # sickle-cell HBB vs wild-type P68871 (live APIs)
    python demo_client_ly.py --replay fixtures                 # same, offline from recorded responses
    python demo_client_ly.py --dna examples/tp53_R175H_cds.fasta --reference P04637
    python demo_client_ly.py --errors                          # also show structured error responses
    python demo_client_ly.py --replay fixtures --errors --save example_outputs
"""

import argparse
import json
import os
import sys
from pathlib import Path

import anyio
from mcp import Client, StdioServerParameters

HERE = Path(__file__).resolve().parent


def banner(text: str) -> None:
    print("\n" + "=" * 78 + f"\n{text}\n" + "=" * 78)


async def run(args: argparse.Namespace) -> int:
    server_args = [str(HERE / "mcp_server_ly.py")]
    if args.replay:
        server_args += ["--replay", str(Path(args.replay).resolve())]
    params = StdioServerParameters(
        command=sys.executable, args=server_args, env={**os.environ, "AFMCP_LOG_LEVEL": "WARNING"}, cwd=str(HERE)
    )
    save_dir = Path(args.save) if args.save else None
    if save_dir:
        save_dir.mkdir(parents=True, exist_ok=True)

    async def call(client: Client, name: str, arguments: dict, label: str) -> dict:
        result = await client.call_tool(name, arguments)
        content = result.structured_content or {}
        if save_dir:
            record = {"tool": name, "arguments": arguments, "isError": result.is_error, "structuredContent": content}
            (save_dir / f"{label}.json").write_text(json.dumps(record, indent=2) + "\n")
        if args.json:
            print(json.dumps(content, indent=2))
        return content

    async with Client(params) as client:
        tools = (await client.list_tools()).tools
        banner(f"Connected to MCP server - {len(tools)} tools")
        for t in tools:
            print(f"  * {t.name}: {t.title}")

        dna_text = Path(args.dna).read_text()
        banner(f"1) find_alphafold_models_for_dna  ({Path(args.dna).name})")
        found = await call(
            client, "find_alphafold_models_for_dna", {"dna_sequence": dna_text, "taxon_id": args.taxon},
            f"find_{Path(args.dna).stem}",
        )
        if not found.get("ok"):
            print("ERROR:", found["error"]["code"], "-", found["error"]["message"])
            return 1
        d = found["data"]
        tr = d["translation"]
        print(f"ORF: frame {tr['reading_frame']}, nt {tr['nt_start']}-{tr['nt_end']}, {tr['protein_length']} aa")
        print(f"Protein: {tr['protein_sequence'][:60]}{'...' if tr['protein_length'] > 60 else ''}")
        print(f"AlphaFold DB: {d['match_status']} ({d['total_matches']} identical model(s))")
        for m in d["matches"][:5]:
            print(
                f"  - {m['model_id']:<20} {m['uniprot_accession'] or '':<11} {(m['organism'] or '?'):<24} "
                f"{'Swiss-Prot' if m['reviewed'] else 'TrEMBL':<10} mean pLDDT {m['mean_plddt']}"
            )
        for w in found["warnings"]:
            print("  warning:", w)

        reference = args.reference or (d["best_match"] or {}).get("uniprot_accession")
        if not reference:
            print("\nNo identical model; pass --reference <UniProt accession> to compare with a known protein.")
            return 0

        banner(f"2) compare_sequence_to_alphafold  (query vs {reference})")
        compared = await call(
            client, "compare_sequence_to_alphafold",
            {"uniprot_accession": reference, "query_sequence": dna_text, "query_type": "dna"},
            f"compare_{Path(args.dna).stem}_vs_{reference}",
        )
        if not compared.get("ok"):
            print("ERROR:", compared["error"]["code"], "-", compared["error"]["message"])
            return 1
        c = compared["data"]
        print(f"Model {c['model']['entry_id']} (v{c['model']['model_version']}, {c['model']['organism']}): {c['model']['urls']['pdb']}")
        print(c["confidence"]["interpretation"])
        print(c["comparison"]["interpretation"])
        for s in c["comparison"]["substitutions"][:10]:
            am = (
                f"AlphaMissense {s['alphamissense_pathogenicity']} ({s['alphamissense_class']})"
                if s["alphamissense_pathogenicity"] is not None
                else "AlphaMissense n/a"
            )
            print(f"  {s['variant']:<8} pLDDT {s['plddt']} ({s['plddt_band']}), BLOSUM62 {s['blosum62']}, {am}")

        if args.errors:
            banner("3) Error handling - every failure is a structured result, never a crash")
            cases = [
                ("error_invalid_dna", "find_alphafold_models_for_dna", {"dna_sequence": "THIS IS NOT A DNA SEQUENCE AT ALL!!"}),
                ("error_missing_parameter", "find_alphafold_models_for_dna", {"reading_frame": "+1"}),
                ("error_no_orf", "find_alphafold_models_for_dna", {"dna_sequence": "TTAG" * 25}),
                ("error_bad_accession", "compare_sequence_to_alphafold", {"uniprot_accession": "HBB_HUMAN!"}),
                ("error_not_found", "compare_sequence_to_alphafold", {"uniprot_accession": "O00000"}),
            ]
            for label, name, arguments in cases:
                out = await call(client, name, arguments, label)
                err = out.get("error") or {}
                print(f"  {label:<24} -> {err.get('code')}: {err.get('message', '')[:110]}")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Demo MCP client for mcp_server_ly.py")
    parser.add_argument("--dna", default=str(HERE / "examples" / "hbb_sickle_E7V_cds.fasta"), help="FASTA/plain DNA file")
    parser.add_argument("--reference", help="UniProt accession to compare against (default: best exact match)")
    parser.add_argument("--taxon", type=int, default=9606, help="preferred NCBI taxon id (default 9606, human)")
    parser.add_argument("--replay", metavar="DIR", help="run the server offline from recorded responses (e.g. fixtures)")
    parser.add_argument("--errors", action="store_true", help="also demonstrate structured error responses")
    parser.add_argument("--json", action="store_true", help="print the raw structured JSON of every call")
    parser.add_argument("--save", metavar="DIR", help="write every tool call + result to DIR/*.json")
    args = parser.parse_args()
    if not args.reference and Path(args.dna).name == "hbb_sickle_E7V_cds.fasta":
        args.reference = "P68871"  # compare the sickle allele with the wild-type model by default
    sys.exit(anyio.run(run, args))


if __name__ == "__main__":
    main()
