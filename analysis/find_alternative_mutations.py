from typing import List, Optional, Union
import vcf, re, os, tempfile, math
import pandas as pd
from difflib import SequenceMatcher
from Levenshtein import ratio


def filter_vcf_by_chromosome(vcf_path: str, chromosome: str) -> List[str]:
    """
    Filters VCF for lines that contain mutations for a specific chromosome using regex
    """
    filtered_records = []
    pattern = re.compile(rf"(?:^|[^0-9]){chromosome}(?![0-9])", re.IGNORECASE)

    with open(vcf_path, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            # Extract the CHROM column (the first one)
            columns = line.split('\t')
            if not columns:
                continue

            # Check if our pattern exists in that specific column
            chrom_column = columns[0]
            if pattern.search(chrom_column):
                filtered_records.append(line)

    return filtered_records


def is_99_percent_match(vcf_ref: str, vcf_alt: str, aln_seq1: str, aln_seq2: str) -> bool:
    """
    Checks if the VCF mutation matches the alignment sequences using Levenshtein ratio.
    This handles shifts and small length differences better than SequenceMatcher.
    """
    # Ensure inputs are strings and uppercase
    vcf_ref, vcf_alt = str(vcf_ref).upper(), str(vcf_alt).upper()
    aln_seq1, aln_seq2 = str(aln_seq1).upper(), str(aln_seq2).upper()

    # Compare Ref to Sequence1 and Alt to Sequence2
    # Levenshtein ratio is: (sum of lengths - edit_distance) / sum of lengths
    ref_sim = ratio(vcf_ref, aln_seq1)
    alt_sim = ratio(vcf_alt, aln_seq2)

    # Both alleles must meet the 99% threshold
    return ref_sim >= 0.99 and alt_sim >= 0.99


def find(vcf_input: str, parsed_bedpe: str, output_folder: str,
        chromosome: str, threshold: Optional[int] = None) -> None:
    """
        Uses filter_vcf_by_chromosome to get relevant records, then maps them
        to the alignment to find nearby mutations, saving each to its own file.
    """

    # 1. Set up the specific subfolder for individual mutation files
    alt_mut_dir = os.path.join(output_folder, "alternative_mutations")
    os.makedirs(alt_mut_dir, exist_ok=True)

    # 2. Use your existing filtering function to get the relevant VCF lines
    vcf_lines = filter_vcf_by_chromosome(vcf_input, chromosome)
    if not vcf_lines:
        print(f"No VCF records found for chromosome {chromosome}.")
        return

    # 3. Load the alignment differences (BEDPE)
    try:
        df_diffs = pd.read_csv(parsed_bedpe, sep='\t')
    except Exception as e:
        print(f"Error reading BEDPE file: {e}")
        return

    if df_diffs.empty:
        print("The alignment file is empty. No mutations to compare.")
        return

    found_count = 0

    # 4. Process each filtered VCF record
    for line in vcf_lines:
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue

        vcf_start = int(parts[1]) - 1
        vcf_ref = parts[3]
        vcf_alt = parts[4]

        # 5. Calculate distance from this VCF site to every difference in the alignment
        # We compare VCF position (Seq1) to alignment position (start1)
        distances = (df_diffs['start1'] - vcf_start).abs()
        if threshold:
            mask_dist = (distances > 0) & (distances <= threshold)
        else:
            mask_dist = distances > 0

        if mask_dist.any():
            # Get candidates within distance
            candidates = df_diffs[mask_dist].copy()
            candidates['dist_to_vcf'] = distances[mask_dist]

            # 6. Apply the 99% Similarity Filter
            # We filter the 'nearby' mutations to only those that match the VCF alleles
            match_mask = candidates.apply(
                lambda row: is_99_percent_match(
                    vcf_ref, vcf_alt,
                    row['seq1'], row['seq2']
                ), axis=1
            )

            nearby_muts = candidates[match_mask]

            if not nearby_muts.empty:
                # Save the individual file
                file_name = f"mut_pos_{vcf_start}.tsv"
                file_path = os.path.join(alt_mut_dir, file_name)

                # Add metadata for context
                output_df = nearby_muts.copy()
                output_df['vcf_ref'] = vcf_ref
                output_df['vcf_alt'] = vcf_alt
                output_df['vcf_start'] = vcf_start

                output_df.to_csv(file_path, sep='\t', index=False)
                found_count += 1

    print(f"Found {len(vcf_lines)} records on chromosome {chromosome}.")
    print(f"Created {found_count} individual mutation files in: {alt_mut_dir}")
