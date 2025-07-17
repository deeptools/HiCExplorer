import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two loop files (BEDGRAPH) and their associated matrices."
    )
    parser.add_argument(
        '--target-loop-file', required=True, help='Target loop file in BEDGRAPH format'
    )
    parser.add_argument(
        '--target-matrix', required=True, help='Matrix associated with target loop file'
    )
    parser.add_argument(
        '--control-loop-file', required=True, help='Control loop file in BEDGRAPH format'
    )
    parser.add_argument(
        '--control-matrix', required=True, help='Matrix associated with control loop file'
    )
    parser.add_argument(
        '--threads', '-t', type=int, default=1, help='Number of threads to use (default: 1)'
    )
    parser.add_argument(
        '--p-value', type=float, default=0.05, help='P-value threshold for significance (default: 0.05)'
    )
    parser.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}',
        help='Show program\'s version number and exit'
    )
    parser.add_argument(
        '--window-size', type=int, default=5,
        help='Window size (in bins) around the loop center to extract (default: 10)'
    )
    parser.add_argument(
        '--peak-size', type=int, default=2,
        help='Peak size (in bins) for the central region (default: 2)'
    )
    parser.add_argument(
        '--output-prefix', required=True, help='Prefix for output files (e.g., results will be saved as <prefix>_target.tsv, etc.)'
    )
    # parser.add_argument(
    #     '-h', '--help', action='help', default=argparse.SUPPRESS,
    #     help='Show this help message and exit'
    # )
    return parser.parse_args()

# Normalize all regions
def normalize(arr):
    arr = np.array(arr)
    return (arr - arr.min()) / (arr.max() - arr.min()) if arr.max() > arr.min() else arr
    
def extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size, peak_size):
    indices_x = hic_matrix_target.getRegionBinRange(chrom1, start1, end1)[0]
    indices_y = hic_matrix_target.getRegionBinRange(chrom2, start2, end2)[0]
    
    start_idx_x = max(indices_x - window_size, 0)
    end_idx_x = indices_x + window_size + 1
    # For control matrix
    start_idx_y = max(indices_y - window_size, 0)
    end_idx_y = indices_y + window_size + 1
    pixels_target = hic_matrix_target.matrix[start_idx_x:end_idx_x, start_idx_y:end_idx_y].toarray()

    pixels_control = hic_matrix_control.matrix[start_idx_x:end_idx_x, start_idx_y:end_idx_y].toarray()
    
    # print(f"Target pixels shape: {pixels_target.shape}")
    # print(f"Control pixels shape: {pixels_control.shape}")
    center = window_size


    peak_start = center - peak_size
    peak_end = center + peak_size + 1
    peak_region_target = pixels_target[peak_start:peak_end, peak_start:peak_end].flatten()

    peak_start = center - peak_size
    peak_end = center + peak_size + 1
    peak_region_control = pixels_control[peak_start:peak_end, peak_start:peak_end].flatten()

    
    # # Donut test data extraction for target
    # horizontal_target = []
    # # top middle

    # # print("Pixels target shape:", pixels_target.shape)
    # # print("Pixels control shape:", pixels_control.shape)
    # # print("Center index:", center)
    # # print("Window size:", window_size)
    # # print("Peak size:", peak_size)

    # # print("Horizontal extraction (target) I : {}".format(pixels_target[:center - peak_size, center - peak_size:center + peak_size + 1]))
    # # print("Horizontal extraction (target) II : {}".format(pixels_target[center + peak_size + 1:, center - peak_size:center + peak_size + 1]))
    # horizontal_target.extend(pixels_target[:center - peak_size, center - peak_size:center + peak_size + 1].flatten())
    # # bottom middle
    # horizontal_target.extend(pixels_target[center + peak_size + 1:, center - peak_size:center + peak_size + 1].flatten())
    # horizontal_target = np.array(horizontal_target).flatten()


    # # print("Vertical extraction (target) I : {}".format(pixels_target[center - peak_size:center + peak_size + 1, :center - peak_size]))
    # # print("Vertical extraction (target) II : {}".format(pixels_target[center - peak_size:center + peak_size + 1, center + peak_size + 1:]))
    # vertical_target = []
    # # left
    # vertical_target.extend(pixels_target[center - peak_size:center + peak_size + 1, :center - peak_size].flatten())
    # # right
    # vertical_target.extend(pixels_target[center - peak_size:center + peak_size + 1, center + peak_size + 1:].flatten())
    # vertical_target = np.array(vertical_target).flatten()

    # # print("Vertical extraction (target) shape:", vertical_target.shape)
    # # print("Horizontal extraction (target) shape:", horizontal_target.shape)

    # # print("Peak region shape (target):", peak_region_target.shape)

    # # Donut test data extraction for control
    # horizontal_control = []
    # horizontal_control.extend(pixels_control[:center - peak_size, center - peak_size:center + peak_size + 1].flatten())
    # horizontal_control.extend(pixels_control[center + peak_size + 1:, center - peak_size:center + peak_size + 1].flatten())
    # horizontal_control = np.array(horizontal_control).flatten()

    # vertical_control = []
    # vertical_control.extend(pixels_control[center - peak_size:center + peak_size + 1, :center - peak_size].flatten())
    # vertical_control.extend(pixels_control[center - peak_size:center + peak_size + 1, center + peak_size + 1:].flatten())
    # vertical_control = np.array(vertical_control).flatten()


    
    # # peak_region_target = normalize(peak_region_target)
    # # peak_region_control = normalize(peak_region_control)
    # # horizontal_target = normalize(horizontal_target)
    # # horizontal_control = normalize(horizontal_control)
    # vertical_target = normalize(vertical_target)
    # vertical_control = normalize(vertical_control)

    # Compute ranksums for all regions
    stat_peak, pval_peak = ranksums(peak_region_target, peak_region_control)
    stat_region, pval_region = ranksums(peak_region_target, peak_region_control)
    # stat_horizontal, pval_horizontal = ranksums(horizontal_target, horizontal_control)
    # stat_vertical, pval_vertical = ranksums(vertical_target, vertical_control)

    # Combine results (example: return as tuple)
    stat = (stat_peak, stat_region)#, stat_horizontal, stat_vertical)
    pval = (pval_peak, pval_peak)#, pval_horizontal, pval_vertical)

    return stat, pval

def detect_differential_loops(target_loop_file, target_matrix, control_loop_file, control_matrix, threads, p_value, window_size, peak_size, output_prefix):
    # Placeholder for the actual implementation
    print(f"Detecting differential loops between {target_loop_file} and {control_loop_file} using matrices {target_matrix} and {control_matrix}.")
    print(f"Using {threads} cores with a p-value threshold of {p_value}.")

    loops_target = pd.read_csv(target_loop_file, sep='\t', header=None)
    loops_control = pd.read_csv(control_loop_file, sep='\t', header=None)

    # Read matrices
    hic_matrix_target = hm.hiCMatrix(target_matrix)
    hic_matrix_control = hm.hiCMatrix(control_matrix)

    print(f"Loaded {len(loops_target)} loops from {target_loop_file}")
    print(f"Loaded {len(loops_control)} loops from {control_loop_file}")
    print(f"Target matrix shape: {hic_matrix_target.matrix.shape}")
    print(f"Control matrix shape: {hic_matrix_control.matrix.shape}")

    # Find exact matching loops between the two files
    matching_loops = pd.merge(
        loops_target, loops_control,
        on=[0, 1, 2, 3, 4, 5],  # Assuming BEDGRAPH columns: chrom1, start1, chrom2, start2
        how='inner'
    )
    # print(f"Found {len(matching_loops)} exact matching loops between the two files.")

    # Get non-matching loops from each file
    non_matching_loops_target = pd.merge(
        loops_target, loops_control,
        on=[0, 1, 2, 3, 4, 5],
        how='left', indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])

    non_matching_loops_control = pd.merge(
        loops_control, loops_target,
        on=[0, 1, 2, 3, 4, 5],
        how='left', indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])

    print(f"Found {len(non_matching_loops_target)} loops unique to {target_loop_file}")
    print(f"Found {len(non_matching_loops_control)} loops unique to {control_loop_file}")

    results_target = []

    # Iterate through non-matching loops in target file
    for idx, loop in non_matching_loops_target.iterrows():
        chrom1, start1, end1, chrom2, start2, end2 = loop[0], loop[1], loop[2], loop[3], loop[4], loop[5]
        # print(f"Unique loop in target file: {chrom1}:{start1}-{chrom2}:{start2}")

        stat, pval = extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size, peak_size)
        results = []
        results_target.append({
            'chrom1': chrom1,
            'start1': start1,
            'end1': end1,
            'chrom2': chrom2,
            'start2': start2,
            'end2': end2,
            'statistic_peak': stat[0],
            'statistic_region': stat[1],
            'pvalue_peak': pval[0],
            'pvalue_region': pval[1]
        })

    # After the loop, create a DataFrame
    results_target_df = pd.DataFrame(results_target)
    # print(results_target_df)

    results_control = []

    # Iterate through non-matching loops in control file
    for idx, loop in non_matching_loops_control.iterrows():
        chrom1, start1, end1, chrom2, start2, end2 = loop[0], loop[1], loop[2], loop[3], loop[4], loop[5]
        # print(f"Unique loop in target file: {chrom1}:{start1}-{chrom2}:{start2}")

        stat, pval = extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size, peak_size)
        results_control.append({
            'chrom1': chrom1,
            'start1': start1,
            'end1': end1,
            'chrom2': chrom2,
            'start2': start2,
            'end2': end2,
            'statistic_peak': stat[0],
            'statistic_region': stat[1],
            'pvalue_peak': pval[0],
            'pvalue_region': pval[1]
            # 'pvalue_vertical': pval[2]
        })

    results_control_df = pd.DataFrame(results_control)

    results_matching = []

    # Iterate through non-matching loops in control file
    for idx, loop in matching_loops.iterrows():
        chrom1, start1, end1, chrom2, start2, end2 = loop[0], loop[1], loop[2], loop[3], loop[4], loop[5]
        # print(f"Unique loop in target file: {chrom1}:{start1}-{chrom2}:{start2}")

        stat, pval = extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size, peak_size)
        results_matching.append({
            'chrom1': chrom1,
            'start1': start1,
            'end1': end1,
            'chrom2': chrom2,
            'start2': start2,
            'end2': end2,
            'statistic_peak': stat[0],
            'statistic_region': stat[1],
            'pvalue_peak': pval[0],
            'pvalue_region': pval[1]
        })

    results_matching_df = pd.DataFrame(results_matching)

    # print(results_control_df)

    # Count significant loops for each p-value type
    # Count significant loops for each p-value type
    num_significant_target_peak = (results_target_df['pvalue_peak'] < p_value).sum()
    num_significant_target_region = (results_target_df['pvalue_region'] < p_value).sum()

    num_significant_control_peak = (results_control_df['pvalue_peak'] < p_value).sum()
    num_significant_control_region = (results_control_df['pvalue_region'] < p_value).sum()

    num_significant_matching_peak = (results_matching_df['pvalue_peak'] < p_value).sum()
    num_significant_matching_region = (results_matching_df['pvalue_region'] < p_value).sum()



    # Apply FDR correction to each p-value type
    results_target_df['fdr_peak'] = multipletests(results_target_df['pvalue_peak'], method='fdr_bh', alpha=p_value)[1]
    results_target_df['fdr_region'] = multipletests(results_target_df['pvalue_region'], method='fdr_bh', alpha=p_value)[1]

    results_control_df['fdr_peak'] = multipletests(results_control_df['pvalue_peak'], method='fdr_bh', alpha=p_value)[1]
    results_control_df['fdr_region'] = multipletests(results_control_df['pvalue_region'], method='fdr_bh', alpha=p_value)[1]

    results_matching_df['fdr_peak'] = multipletests(results_matching_df['pvalue_peak'], method='fdr_bh', alpha=p_value)[1]
    results_matching_df['fdr_region'] = multipletests(results_matching_df['pvalue_region'], method='fdr_bh', alpha=p_value)[1]

    num_significant_target_fdr_peak = (results_target_df['fdr_peak'] < p_value).sum()
    num_significant_target_fdr_region = (results_target_df['fdr_region'] < p_value).sum()

    num_significant_control_fdr_peak = (results_control_df['fdr_peak'] < p_value).sum()
    num_significant_control_fdr_region = (results_control_df['fdr_region'] < p_value).sum()

    num_significant_matching_fdr_peak = (results_matching_df['fdr_peak'] < p_value).sum()
    num_significant_matching_fdr_region = (results_matching_df['fdr_region'] < p_value).sum()

    print(f"Number of unique target loops with FDR peak < {p_value}: {num_significant_target_fdr_peak}")
    print(f"Number of unique target loops with FDR region < {p_value}: {num_significant_target_fdr_region}")

    print(f"Number of unique control loops with FDR peak < {p_value}: {num_significant_control_fdr_peak}")
    print(f"Number of unique control loops with FDR region < {p_value}: {num_significant_control_fdr_region}")

    print(f"Number of matching loops with FDR peak < {p_value}: {num_significant_matching_fdr_peak}")
    print(f"Number of matching loops with FDR region < {p_value}: {num_significant_matching_fdr_region}")

    print("\n\n")

    num_significant_target_fdr_all = ((results_target_df['fdr_peak'] < p_value) &
                                      (results_target_df['fdr_region'] < p_value)).sum()

    num_significant_control_fdr_all = ((results_control_df['fdr_peak'] < p_value) &
                                       (results_control_df['fdr_region'] < p_value)).sum()

    num_significant_matching_fdr_all = ((results_matching_df['fdr_peak'] < p_value) &
                                        (results_matching_df['fdr_region'] < p_value)).sum()

    print(f"Number of unique target loops with both FDRs < {p_value}: {num_significant_target_fdr_all}")
    print(f"Number of unique control loops with both FDRs < {p_value}: {num_significant_control_fdr_all}")
    print(f"Number of matching loops with both FDRs < {p_value}: {num_significant_matching_fdr_all}")

    print("\n\n")
    print(f"Total number of differential loops FDR: {num_significant_target_fdr_all + num_significant_control_fdr_all + num_significant_matching_fdr_all}")

    def chrom_sort_key(row):
        def parse_chrom(chrom):
            chrom = str(chrom)
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            # Try to convert to int, fallback to string for X/Y/M
            try:
                return (0, int(chrom))
            except ValueError:
                # X, Y, M, etc. Sort X=23, Y=24, M=25, else string
                special = {'X': 23, 'Y': 24, 'M': 25}
                return (1, special.get(chrom.upper(), chrom))
        return (
            parse_chrom(row['chrom1']),
            int(row['start1']),
            int(row['end1']),
            parse_chrom(row['chrom2']),
            int(row['start2']),
            int(row['end2'])
        )

    results_target_df = results_target_df.sort_values(
        by=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'],
        key=lambda cols: [chrom_sort_key(row) for _, row in results_target_df.iterrows()]
    ).reset_index(drop=True)

    results_control_df = results_control_df.sort_values(
        by=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'],
        key=lambda cols: [chrom_sort_key(row) for _, row in results_control_df.iterrows()]
    ).reset_index(drop=True)

    results_matching_df = results_matching_df.sort_values(
        by=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'],
        key=lambda cols: [chrom_sort_key(row) for _, row in results_matching_df.iterrows()]
    ).reset_index(drop=True)
    # Write results to files using the output prefix
    results_target_df.to_csv(f"{output_prefix}_target.tsv", sep='\t', index=False)
    results_control_df.to_csv(f"{output_prefix}_control.tsv", sep='\t', index=False)
    results_matching_df.to_csv(f"{output_prefix}_matching.tsv", sep='\t', index=False)

    # Write filtered results to three separate files
    results_target_df_filtered = results_target_df[
        (results_target_df['fdr_peak'] < p_value) &
        (results_target_df['fdr_region'] < p_value)
    ]
    results_target_df_filtered.to_csv(f"{output_prefix}_target_fdr.tsv", sep='\t', index=False, header=True)

    results_control_df_filtered = results_control_df[
        (results_control_df['fdr_peak'] < p_value) &
        (results_control_df['fdr_region'] < p_value)
    ]
    results_control_df_filtered.to_csv(f"{output_prefix}_control_fdr.tsv", sep='\t', index=False, header=True)

    results_matching_df_filtered = results_matching_df[
        (results_matching_df['fdr_peak'] < p_value) &
        (results_matching_df['fdr_region'] < p_value)
    ]
    results_matching_df_filtered.to_csv(f"{output_prefix}_matching_fdr.tsv", sep='\t', index=False, header=True)

    # Write out three files with only chrom1, start1, end1, chrom2, start2, end2 and no header
    results_target_df_filtered[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].to_csv(
        f"{output_prefix}_target_fdr_loops.bed", sep='\t', index=False, header=False
    )
    results_control_df_filtered[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].to_csv(
        f"{output_prefix}_control_fdr_loops.bed", sep='\t', index=False, header=False
    )
    results_matching_df_filtered[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].to_csv(
        f"{output_prefix}_matching_fdr_loops.bed", sep='\t', index=False, header=False
    )

def main(args=None):
    args = parse_args()
    detect_differential_loops(
        args.target_loop_file,
        args.target_matrix,
        args.control_loop_file,
        args.control_matrix,
        args.threads,
        args.p_value,
        args.window_size,
        args.peak_size,
        args.output_prefix)
