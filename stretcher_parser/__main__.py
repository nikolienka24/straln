import sys
import os
import argparse

import pandas as pd

from . import parser, reformat
from .utils import save_statistics, get_offsets_from_file
from analysis import find_alternative_mutations


def handle_parse(args):
    """Logic for the 'parse' command"""
    output_folder = args.output_folder or "./stralln_parse_results"
    os.makedirs(output_folder, exist_ok=True)

    # Handle offsets
    if args.offset1 is not None and args.offset2 is not None:
        offset1, offset2 = args.offset1, args.offset2
    else:
        offset1, offset2 = get_offsets_from_file(args.aln_input_file)
        # Manual override if only one was provided
        offset1 = args.offset1 if args.offset1 is not None else offset1
        offset2 = args.offset2 if args.offset2 is not None else offset2
        print(f"[*] Using offsets: {offset1}, {offset2}")

    # 1. Run Parser
    parsed_file = os.path.join(output_folder, "parsed.bedpe")
    length, mismatches, gaps = parser.run(
        args.aln_input_file, parsed_file, args.seq_name1, args.seq_name2, offset1, offset2
    )

    # 2. Save Stats
    output_stats = os.path.join(output_folder, "stats.txt")
    save_statistics(args.seq_name1, args.seq_name2, length, mismatches, gaps, output_stats)

    # 3. Join rows
    parsed_joined_file = os.path.join(output_folder, "parsed.joined.bedpe")
    reformat.join_consecutive_rows(parsed_file, parsed_joined_file)

    # 4. Convert to BED
    bed1 = os.path.join(output_folder, f"{args.seq_name1}.bed")
    bed2 = os.path.join(output_folder, f"{args.seq_name2}.bed")
    reformat.bedpe_to_bed(parsed_joined_file, bed1, bed2)

    print(f"[✔] Parsing complete. Files saved to: {output_folder}")


def handle_overlap(args):
    """Logic for the 'overlap' command"""
    output_folder = args.output_folder or "./stralln_alternative_mutations"
    os.makedirs(output_folder, exist_ok=True)
    
    find_alternative_mutations.find(
        args.vcf, args.bedpe, output_folder, args.chrom, args.distance
    )
    print(f"[✔] Overlap analysis complete. Results in: {output_folder}")


def handle_snap(args):
    """Logic for the 'snap' command (Placeholder for your IGV logic)"""
    print(f"[*] Snapping screenshots using {args.config}...")
    # Your IGV snapping logic here
    pass


def main():
    parser_main = argparse.ArgumentParser(
        description="stralln: Stretcher Alignment Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser_main.add_subparsers(dest="command", help="Available commands")

    # ----- COMMAND: parse -----
    p_parse = subparsers.add_parser("parse", help="Convert .aln to BEDPE/BED and merge consecutive rows.")
    p_parse.add_argument("aln_input_file", help="Path to EMBOSS .aln file")
    p_parse.add_argument("-o", "--output_folder", default="./stralln_results", help="Output directory")
    p_parse.add_argument("-s1", "--seq_name1", default="seq1", help="Label for seq1")
    p_parse.add_argument("-s2", "--seq_name2", default="seq2", help="Label for seq2")
    p_parse.add_argument("-off1", "--offset1", type=int, help="Manual offset for seq1")
    p_parse.add_argument("-off2", "--offset2", type=int, help="Manual offset for seq2")
    p_parse.set_defaults(func=handle_parse)

    # ----- COMMAND: overlap -----
    p_overlap = subparsers.add_parser("overlap", help="Find alternative mutations from VCF file.")
    p_overlap.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    p_overlap.add_argument("-b", "--bedpe", required=True, help="Input parsed BEDPE file")
    p_overlap.add_argument("-c", "--chrom", required=True, help="Chromosome")
    p_overlap.add_argument("-d", "--distance", type=int, default=100, help="Window size (bp)")
    p_overlap.add_argument("-o", "--output_folder", default=".", help="Output directory")
    p_overlap.set_defaults(func=handle_overlap)

    # ----- COMMAND: snap -----
    p_snap = subparsers.add_parser("snap", help="Generate IGV snapshots")
    p_snap.add_argument("config", help="Path to setup.json")
    p_snap.add_argument("-o", "--outdir", default="./snapshots", help="Output directory")
    p_snap.set_defaults(func=handle_snap)

    args = parser_main.parse_args()

    if not args.command:
        parser_main.print_help()
        sys.exit(1)

    # Route to the appropriate function
    args.func(args)


if __name__ == "__main__":
    main()