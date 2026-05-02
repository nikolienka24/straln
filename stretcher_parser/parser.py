from typing import Tuple, List, Optional, Union
from stretcher_parser.utils import get_offsets_from_file


def _check_sequences(sequence_name1, sequence_name2,
                     seq1_data, seq2_data,
                     found_start, buffer,
                     prefix_a, prefix_b):
    output = ""
    p1, p2 = seq1_data[1], seq2_data[1]
    seq1_chunk, seq2_chunk = seq1_data[0], seq2_data[0]

    line_len, line_len_curr = 0, 0
    line_mismatches, line_gaps = 0, 0

    for a, b in zip(seq1_chunk, seq2_chunk):
        # 1. Pozície PRED spracovaním písmena
        current_p1_start = p1
        current_p2_start = p2

        # 2. Posun počítadla (ak nie je gap)
        if a != "-": p1 += 1
        if b != "-": p2 += 1

        # 3. Hľadanie začiatku reálneho alignmentu
        if not found_start:
            if a != "-" and b != "-":
                found_start = True
            else:
                continue

        line_len_curr += 1

        # 4. Ak nájdeme nesúlad (SNP alebo Indel)
        if a != b:
            if a == "-" or b == "-":
                # INDEL
                s1, e1 = current_p1_start - 1, p1
                s2, e2 = current_p2_start - 1, p2
                out1 = (prefix_a if prefix_a else "") + (a if a != "-" else "")
                out2 = (prefix_b if prefix_b else "") + (b if b != "-" else "")
            else:
                # SNP
                s1, e1 = current_p1_start, p1
                s2, e2 = current_p2_start, p2
                out1, out2 = a, b

            buffer.append(f"{sequence_name1}\t{s1}\t{e1}\t{sequence_name2}\t{s2}\t{e2}\t{out1}\t{out2}\n")
            line_mismatches += 1
            if a == "-" or b == "-": line_gaps += 1

        # 5. DÔLEŽITÉ: Pri zhode (A==A) vyprázdnime buffer do outputu
        if a == b and a != "-":
            if buffer:
                output += "".join(buffer)
                buffer.clear()

            line_len += line_len_curr
            line_len_curr = 0

        # Aktualizácia kotvy pre Indely
        if a != "-": prefix_a = a
        if b != "-": prefix_b = b

    return output, line_len, line_mismatches, line_gaps, found_start, buffer, p1, p2, prefix_a, prefix_b

def _parse(input_file: str, output_file: str,
           seq_name1: str, seq_name2:str,
           offset_seq1: int = 0, offset_seq2: int = 0) -> Tuple[int, int, int]:
    """
    The main engine that opens the alignment file, skips the technical headers,
    and reads the sequences line-by-line to find differences.
    """

    with open(input_file) as infile, open(output_file, 'w') as outfile:
        # BEDPE header with 2 additional columns
        outfile.write("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tnucleotide1\tnucleotide2\n")

        found_start = False
        buffer = []
        seq_len, mismatches, gaps = 0, 0, 0

        # read header
        line = infile.readline()
        while line.startswith("#"):
            line = infile.readline()

        line = infile.readline()
        while line.startswith("#"):
            line = infile.readline()

        infile.readline()  # skip blank line at the end of header
        # end read header

        curr1, curr2 = offset_seq1, offset_seq2
        prefix_a, prefix_b = None, None
        while True:
            if line.startswith("#") or not line:
                break

            seq1_line = infile.readline().strip().split()
            infile.readline()  # match line
            seq2_line = infile.readline().strip().split()
            infile.readline()  # lower pos
            infile.readline()  # blank line

            if seq1_line[0].startswith("#") or seq2_line[0].startswith("#"):
                break

            seq1_seq = seq1_line[1]
            seq2_seq = seq2_line[1]

            if not seq2_seq:
                break

            out, line_length, mm_add, gaps_add, found_start, buffer, curr1, curr2, prefix_a, prefix_b = _check_sequences(
                seq_name1, seq_name2,
                [seq1_seq, curr1],
                [seq2_seq, curr2],
                found_start, buffer,
                prefix_a, prefix_b
            )
            outfile.write(out)

            seq_len += line_length
            mismatches += mm_add
            gaps += gaps_add

            line = infile.readline()
            if not line:
                break

        return seq_len, mismatches, gaps


def run(in_file: str, out_file: str,
        seq_name1: str, seq_name2: str,
        offset1: int, offset2: int) -> Tuple[int, int, int]:
    """
    Starts the parsing process and returns a summary of the results,
    including how long the sequences are and how many errors were found.
    """
    length, mismatches, gaps = _parse(in_file, out_file, seq_name1, seq_name2, offset1, offset2)
    return length, mismatches, gaps
