from typing import Tuple, Optional


def save_statistics(seq_name1: str, seq_name2: str,
                    length: int, mismatches: int, gaps: int,
                    output_file: str) -> None:
    """
    Creates a summary report. It calculates the error percentages (mismatch and gap rates)
    and saves all the final numbers into a text file.
    """

    with open(output_file, "w") as f:
        f.write(f"Alignment statistics between '{seq_name1}' and '{seq_name2}':\n")
        f.write(f"  Total aligned positions (excluding ignored ends): {length}\n")
        f.write(f"  Base mismatches (including gaps): {mismatches}\n")
        f.write(f"  Gaps detected: {gaps}\n")
        if length > 0:
            mismatch_rate = (mismatches + gaps) / length * 100
            gap_rate = gaps / length * 100
            f.write(f"  Mismatch rate (including gaps): {mismatch_rate:.8f}%\n")
            f.write(f"  Gap rate: {gap_rate:.8f}%\n")


def get_offsets_from_file(aln_file_path: str) -> Tuple[int, int]:
    """
    Reads the header of the Stretcher file to find the starting genomic positions.
    It automatically figures out where your sequences begin on the chromosome.
    """

    offset_seq1 = None
    offset_seq2 = None

    with open(aln_file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("# 1:"):
                # Format: "# 1: 87413227-97757117"
                parts = line.split(":")[1].strip().split("-")
                offset_seq1 = int(parts[0])
            elif line.startswith("# 2:"):
                parts = line.split(":")[1].strip().split("-")
                offset_seq2 = int(parts[0])

            if offset_seq1 is not None and offset_seq2 is not None:
                break

    if offset_seq1 is None or offset_seq2 is None:
        raise ValueError(f"Could not parse offsets from Stretcher file: {aln_file_path}")

    return offset_seq1, offset_seq2

