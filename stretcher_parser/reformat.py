from typing import List, TextIO

### FORMAT CONVERSIONS ===========================================================
def bedpe_to_bed(infile: str, out1: str, out2: str) -> None:
    """
    Splits the combined alignment file into two separate BED files,
    one for each sequence, so they can be viewed individually.
    """

    with open(infile) as f, open(out1, "w") as bed1, open(out2, "w") as bed2:
        header = next(f)
        bed1.write("chrom\tstart\tend\tseq\n")
        bed2.write("chrom\tstart\tend\tseq\n")

        for line in f:
            chrom1, start1, end1, chrom2, start2, end2, seq1, seq2 = line.rstrip().split("\t")
            bed1.write(f"{chrom1}\t{start1}\t{end1}\t{seq1}\n")
            bed2.write(f"{chrom2}\t{start2}\t{end2}\t{seq2}\n")


### REFORMATTING =================================================================
def _flush_buffer(buf: List[List[str]], out_fh: TextIO) -> None:
    """
    Merge consecutive positions in buffer and write a single BEDPE row.
    - Coordinates: start = first position, end = last position
    - Nucleotides: concatenate all nucleotides in the run, ignoring gaps
    """
    if not buf:
        return

    # assume that chrom1 and chrom2 are in the whole file same
    chrom1 = buf[0][0]
    chrom2 = buf[0][3]

    # Get start/end coordinates for ref and alt
    seq1, seq2 = "", ""
    start1, end1 = None, None
    start2, end2 = None, None
    for idx, line in enumerate(buf):
        if idx == 0:
            start1 = line[1]
            start2 = line[4]
            seq1, seq2 = line[6], line[7]
            continue

        if len(line[6]) == 1 and len(line[7]) == 1:
             seq1 += line[6]
             seq2 += line[7]
        elif len(line[6]) == 2 and len(line[7]) == 1:
            seq1 += line[6][1]
        elif len(line[6]) == 1 and len(line[7]) == 2:
            seq2 += line[7][1]

    end1, end2 = buf[-1][2] + 1, buf[-1][5] + 1

    out_fh.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{seq1}\t{seq2}\n")


def _is_consecutive(prev: List, curr: List) -> bool:
    """
    Riadky nadväzujú, ak sa ich súradnice prekrývajú o 1 (anchor)
    alebo na seba plynule nadväzujú (0-based).
    """
    # Rozdiel medzi štartom aktuálneho a koncom predchádzajúceho
    # Pri indeli s kotvou bude tento rozdiel -1 (prekryv)
    # Pri SNP bez kotvy bude tento rozdiel 0
    diff1 = curr[1] - prev[2]
    diff2 = curr[4] - prev[5]

    # Povolený rozsah je -1 (prekryv o 1 bázu) až 0 (plynulé nadviazanie)
    # Zároveň kontrolujeme, či aspoň jedna sekvencia ostáva na mieste (Indel)
    # alebo obe pokračujú (SNP)
    conn1 = -1 <= diff1 <= 0
    conn2 = -1 <= diff2 <= 0

    return conn1 and conn2


def join_consecutive_rows(input_file: str, output_file: str) -> None:
    """
    Reads through the alignment differences and merges nearby changes
    into single rows. This makes the data much cleaner and easier to read.
    """

    with open(input_file, 'r') as in_fh, open(output_file, 'w') as out_fh:
        out_fh.write("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsequence1\tsequence2\n")

        buffer = []
        found_header = False
        for line in in_fh:
            if not found_header or not line.strip():
                found_header = True
                continue

            row = line.strip().split("\t")
            if len(row) < 8:
                continue

            # BEDPE row: chrom1, start1, end1, chrom2, start2, end2, nt1, nt2
            chrom1, start1, end1, chrom2, start2, end2, nt1, nt2 = row
            row_arr = [chrom1, int(start1), int(end1) - 1, chrom2, int(start2), int(end2) - 1, nt1, nt2]

            if not buffer:
                buffer.append(row_arr)
                continue

            if _is_consecutive(buffer[-1], row_arr):
                buffer.append(row_arr)
            else:
                _flush_buffer(buffer, out_fh)
                buffer = [row_arr]

        # Flush remaining buffer
        _flush_buffer(buffer, out_fh)