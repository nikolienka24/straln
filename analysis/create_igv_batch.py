import sys
import os
import json

def create_igv_batch_script(regions_file, output_dir, genome_fasta, tracks, window_padding=50, max_panel_height=2000):
    """Generate IGV batch script from BED/BEDPE file"""
    batch_commands = []

    # Initial setup
    batch_commands.append("new")
    batch_commands.append(f"genome {genome_fasta}")
    
    # Load all tracks (BAMs, VCFs, or BEDs) defined in the JSON
    for track in tracks:
        batch_commands.append(f"load {track}")
        
    batch_commands.append(f"snapshotDirectory {os.path.abspath(output_dir)}")
    batch_commands.append(f"maxPanelHeight {max_panel_height}")

    # Read positions from file
    with open(regions_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            # Handle both BED (3+ cols) and BEDPE (6+ cols)
            if len(parts) < 3:
                continue

            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            region_id = parts[3] if len(parts) > 3 else f"region_{start}"
            
            # IGV uses 1-based coordinates for 'goto', but BED is 0-based
            # Adding padding for visual context
            window_start = max(1, start + 1 - window_padding)
            window_end = end + window_padding

            region = f"{chrom}:{window_start}-{window_end}"
            snapshot_name = f"{region_id}_{chrom}_{start}.png"

            batch_commands.append(f"goto {region}")
            batch_commands.append("sort base")
            batch_commands.append(f"snapshot {snapshot_name}")

    batch_commands.append("exit")
    return '\n'.join(batch_commands)