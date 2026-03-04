import sys
import os
import argparse

import vcf
import pandas as pd

from . import parser, reformat
from .utils import save_statistics, get_offsets_from_file
from analysis import find_alternative_mutations


def parse_arguments() -> argparse.Namespace:
    parser_arg = argparse.ArgumentParser(
        description="""
    straln: Stretcher Alignment & Alternative Mutation Finder
    ---------------------------------------------------------
    Parses .aln files from EMBOSS stretcher into BEDPE/BED formats and optionally 
    compares them with VCF files to identify alternative genomic variations.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # ----- POSITIONAL ARGUMENTS -----
    parser_arg.add_argument("aln_input_file", help="path to the .aln alignment file from EMBOSS stretcher (.aln format)")

    # ----- OPTIONAL ARGUMENTS -----
    parser_arg.add_argument("-o", "--output_folder", type=str, default=None, help="path to output folder (default: current folder)")

    # ----- MUTATION ANALYSIS GROUP (OPTIONAL) -----
    mutation = parser_arg.add_argument_group('Alternative Mutation Analysis (Optional)')
    mutation.add_argument(
        "-vcf", "--vcf_input_file",
        help="input VCF file, providing this triggers the search for alternative mutations"
    )
    parser_arg.add_argument(
        "-c", "--chromosome", 
        help="target chromosome (1-22, X, Y)")
    mutation.add_argument(
        "-d", "--distance", type=int, default=100,
        help="window size (bp) for finding nearby alternative mutations (default: 100)"
    )

    # ----- ALIGNMENT METADATA GROUP -----
    metadata = parser_arg.add_argument_group('Alignment Metadata & Offsets (Optional)')
    metadata.add_argument("-s1", "--seq_name1", default="Seq1", help="label for the first sequence")
    metadata.add_argument("-s2", "--seq_name2", default="Seq2", help="label for the second sequence")
    metadata.add_argument("-off1", "--offset1", type=int, help="starting position for seq1 (if not provided, the tool automatically extracts it from the .aln header)")
    metadata.add_argument("-off2", "--offset2", type=int, help="starting position for seq2 (if not provided, the tool automatically extracts it from the .aln header)")

    return parser_arg.parse_args()


def main() -> None:
    args = parse_arguments()

    # 0. parse input files =============================================================
    aln_input_file = args.aln_input_file

    vcf_input_file = None
    chromosome = None
    if args.vcf_input_file:
        vcf_input_file = args.vcf_input_file
        chromosome = args.chromosome

    output_folder = "."
    if args.output_folder:
        output_folder = args.output_folder
        os.makedirs(output_folder, exist_ok=True)

    seq_name1, seq_name2 = "seq1", "seq2"
    if args.seq_name1:
        seq_name1 = args.seq_name1
    if args.seq_name2:
        seq_name2 = args.seq_name2

    distance_threshold = None
    if args.distance:
        distance_threshold = args.distance

    if args.offset1 is not None and args.offset2 is not None:
        offset1 = args.offset1
        offset2 = args.offset2
    else:
        # If one or both are missing, we fetch them from the alignment file
        offset1, offset2 = get_offsets_from_file(aln_input_file)

        # If the user provided one but not the other, the manual one takes precedence
        if args.offset1 is not None:
            offset1 = args.offset1
        if args.offset2 is not None:
            offset2 = args.offset2

        print(f"Using offsets extracted from alignment header: {offset1}, {offset2}")

    # 1. stretcher parser =================================================================
    parsed_file = output_folder + "/parsed.bedpe"
    length, mismatches, gaps = parser.run(aln_input_file, parsed_file, seq_name1, seq_name2, offset1, offset2)

    output_stats = output_folder + "/stats.txt"
    save_statistics(seq_name1, seq_name2, length, mismatches, gaps, output_stats)

    # 2. join consecutive rows ============================================================
    parsed_joined_file = output_folder + "/parsed.joined.bedpe"
    reformat.join_consecutive_rows(parsed_file, parsed_joined_file)

    # 3. convert bedpe to bed two separate bed files =======================================================
    bed1, bed2 = output_folder + "/seq1.bed", output_folder + "/seq2.bed"                                       
    reformat.bedpe_to_bed(parsed_joined_file, bed1, bed2)

    # 4. find alternative mutations =======================================================
    if vcf_input_file:
        find_alternative_mutations.find(vcf_input_file, parsed_joined_file, output_folder, chromosome, distance_threshold)

    # 5. print final information =======================================================
    print("Output files saved to " + output_folder)


if __name__ == "__main__":
    main()