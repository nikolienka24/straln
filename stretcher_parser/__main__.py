import sys
import os
import argparse
import json

import pandas as pd

from . import parser, reformat
from .utils import save_statistics, get_offsets_from_file
from analysis import find_alternative_mutations, create_igv_batch, swap_bedpe_columns


def handle_parse(args):
    """Logic for the 'parse' command"""
    output_folder = args.output_folder or "./straln_parse_results"
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
    output_folder = args.output_folder
    os.makedirs(output_folder, exist_ok=True)

    bedpe = None
    if args.aln:
        print(f"[*] Argument -a (--aln) detected. Running parser first...")

        # temporary bedpe files
        parsed_file = os.path.join(output_folder, "tmp_parsed.bedpe")
        parsed_joined_file = os.path.join(output_folder, "tmp_parsed.joined.bedpe")
        
        offset1, offset2 = get_offsets_from_file(args.aln)
        parser.run(args.aln, parsed_file, "seq1", "seq2", offset1, offset2)
        reformat.join_consecutive_rows(parsed_file, parsed_joined_file)
        
        bedpe = parsed_joined_file

    elif args.bedpe:
        bedpe = args.bedpe

    print(f"[*] Running overlap analysis with: {bedpe}")
    find_alternative_mutations.find(
        args.vcf, bedpe, output_folder, args.chrom, args.distance, args.percentual_identity
    )

    if args.aln:
        os.remove(parsed_file)
    
    print(f"[✔] Overlap analysis complete. Results in: {output_folder}")


def handle_snap(args):
    """Logic for the 'snap' command: Generates IGV batch script from JSON"""
    print(f"[*] Reading configuration from {args.config}...")
    
        # Use the first argument as the JSON config path
    config_file = args.config

    # Load configuration from JSON
    with open(config_file, 'r') as f:
        config = json.load(f)

    # Extracting data from the new JSON structure
    # We use .get() to provide safe defaults
    paths = config.get('paths', {})
    settings = config.get('igv_settings', {})

    # Generate batch script
    batch_script = create_igv_batch.create_igv_batch_script(
        regions_file=paths.get('regions_file'),
        output_dir=paths.get('output_dir'),
        genome_fasta=paths.get('genome_fasta'),
        tracks=paths.get('tracks', []), # This can now be a list of multiple files
        window_padding=settings.get('window_padding', 50),
        max_panel_height=settings.get('max_panel_height', 2000)
    )

    # Write batch script
    batch_file = paths.get('script_file', 'igv_batch_script.bat')
    with open(batch_file, 'w') as f:
        f.write(batch_script)

    print(f"[✔] IGV batch script created: {batch_file}")
    print(f"[*] How to run: check the documentation")


def handle_swap(args):
    """Logic for the 'swap' command"""
    output_path = args.output_file
    output_dir = os.path.dirname(output_path)
    
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print(f"[*] Swapping columns in: {args.input_file}")
    
    swap_bedpe_columns.swap(args.input_file, output_path)
    
    print(f"[✔] Swap complete. Result saved to: {output_path}")



def main():
    parser_main = argparse.ArgumentParser(
        description="straln: Stretcher Alignment Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser_main.add_subparsers(dest="command", help="Available commands")

    # ----- COMMAND: parse -----
    p_parse = subparsers.add_parser("parse", help="Convert .aln to BEDPE/BED and merge consecutive rows.")
    p_parse.add_argument("aln_input_file", help="Path to EMBOSS .aln file")
    p_parse.add_argument("-o", "--output_folder", default="./straln_results", help="Output directory")
    p_parse.add_argument("-s1", "--seq_name1", default="seq1", help="Label for seq1")
    p_parse.add_argument("-s2", "--seq_name2", default="seq2", help="Label for seq2")
    p_parse.add_argument("-off1", "--offset1", type=int, help="Manual offset for seq1")
    p_parse.add_argument("-off2", "--offset2", type=int, help="Manual offset for seq2")
    p_parse.set_defaults(func=handle_parse)

    # ----- COMMAND: overlap -----
    p_overlap = subparsers.add_parser("overlap", help="Find alternative mutations from VCF file.")
    
    input_group = p_overlap.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-b", "--bedpe", help="Input parsed BEDPE file")
    input_group.add_argument("-a", "--aln", help="Input alignment file")

    p_overlap.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    p_overlap.add_argument("-c", "--chrom", required=True, help="Chromosome")
    p_overlap.add_argument("-d", "--distance", type=int, default=100, help="Window size (bp)")
    p_overlap.add_argument("-o", "--output_folder", default="straln_overlap_result", help="Path to output directory")
    p_overlap.add_argument("-p", "--percentual_identity", type=float, default=99, help="Percentual identity of alternative mutations to the mutation in vcf (float)")
    p_overlap.set_defaults(func=handle_overlap)

    # ----- COMMAND: snap -----
    p_snap = subparsers.add_parser("snap", help="Generate IGV snapshots")
    p_snap.add_argument("config", help="Path to setup.json")
    p_snap.set_defaults(func=handle_snap)

    # ----- COMMAND: swap -----
    p_swap = subparsers.add_parser("swap", help="Swap alignment files.")
    p_swap.add_argument("input_file", help="Path to input BEDPE")
    p_swap.add_argument("-o", "--output_file", default="./swapped.bedpe", help="Path to output file")
    p_swap.set_defaults(func=handle_swap)

    args = parser_main.parse_args()

    if not args.command:
        parser_main.print_help()
        sys.exit(1)

    # Route to the appropriate function
    args.func(args)


if __name__ == "__main__":
    main()