import pandas as pd

def swap(input_file, output_file):
    cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'seq1', 'seq2']
    df = pd.read_csv(input_file, sep='\t', names=cols, header=0)

    new_order = ['chrom2', 'start2', 'end2', 'chrom1', 'start1', 'end1', 'seq2', 'seq1']
    df_swapped = df[new_order]
    
    df_swapped.columns = cols
    
    df_swapped.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"Hotovo! Dáta boli prehodené a uložené do {output_file}")