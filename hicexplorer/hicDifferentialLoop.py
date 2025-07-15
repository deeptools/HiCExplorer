import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

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
    # parser.add_argument(
    #     '-h', '--help', action='help', default=argparse.SUPPRESS,
    #     help='Show this help message and exit'
    # )
    return parser.parse_args()


def extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size):
    indices_target = hic_matrix_target.getRegionBinRange(chrom1, start1, end1)[0]
    indices_control = hic_matrix_control.getRegionBinRange(chrom2, start2, end2)[0]

    # Extract a 21x21 window (+/-10 pixels) around the central bin in both matrices
    # window_size = 10
    # For target matrix
    start_idx_target = max(indices_target - window_size, 0)
    end_idx_target = indices_target + window_size + 1
    pixels_target = hic_matrix_target.matrix[start_idx_target:end_idx_target, start_idx_target:end_idx_target].toarray()
    # For control matrix
    start_idx_control = max(indices_control - window_size, 0)
    end_idx_control = indices_control + window_size + 1
    pixels_control = hic_matrix_control.matrix[start_idx_control:end_idx_control, start_idx_control:end_idx_control].toarray()
    
    flat_target = pixels_target.flatten()
    flat_control = pixels_control.flatten()

    flat_target = (flat_target - flat_target.min()) / (flat_target.max() - flat_target.min()) if flat_target.max() > flat_target.min() else flat_target
    flat_control = (flat_control - flat_control.min()) / (flat_control.max() - flat_control.min()) if flat_control.max() > flat_control.min() else flat_control
    stat, pval = ranksums(flat_target, flat_control)


    # print(f"Rank-sum test statistic: {stat}, p-value: {pval}")
    return stat, pval

def detect_differential_loops(target_loop_file, target_matrix, control_loop_file, control_matrix, threads, p_value, window_size):
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

        stat, pval = extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size)
        results = []
        results_target.append({
            'chrom1': chrom1,
            'start1': start1,
            'end1': end1,
            'chrom2': chrom2,
            'start2': start2,
            'end2': end2,
            'statistic': stat,
            'pvalue': pval
        })

    # After the loop, create a DataFrame
    results_target_df = pd.DataFrame(results_target)
    # print(results_target_df)

    results_control = []

    # Iterate through non-matching loops in control file
    for idx, loop in non_matching_loops_control.iterrows():
        chrom1, start1, end1, chrom2, start2, end2 = loop[0], loop[1], loop[2], loop[3], loop[4], loop[5]
        # print(f"Unique loop in target file: {chrom1}:{start1}-{chrom2}:{start2}")

        stat, pval = extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size)
        results_control.append({
            'chrom1': chrom1,
            'start1': start1,
            'end1': end1,
            'chrom2': chrom2,
            'start2': start2,
            'end2': end2,
            'statistic': stat,
            'pvalue': pval
        })

    results_control_df = pd.DataFrame(results_control)

    results_matching = []

    # Iterate through non-matching loops in control file
    for idx, loop in matching_loops.iterrows():
        chrom1, start1, end1, chrom2, start2, end2 = loop[0], loop[1], loop[2], loop[3], loop[4], loop[5]
        # print(f"Unique loop in target file: {chrom1}:{start1}-{chrom2}:{start2}")

        stat, pval = extract_and_apply_test(hic_matrix_target, hic_matrix_control, chrom1, start1, end1, chrom2, start2, end2, window_size)
        results_matching.append({
            'chrom1': chrom1,
            'start1': start1,
            'end1': end1,
            'chrom2': chrom2,
            'start2': start2,
            'end2': end2,
            'statistic': stat,
            'pvalue': pval
        })

    results_matching_df = pd.DataFrame(results_matching)

    # print(results_control_df)

    num_significant_target = (results_target_df['pvalue'] < p_value).sum()
    num_significant_control = (results_control_df['pvalue'] < p_value).sum()
    num_significant_matching = (results_matching_df['pvalue'] < p_value).sum()
    print(f"Number of unique target loops with p-value < {p_value}: {num_significant_target}")
    print(f"Number of unique control loops with p-value < {p_value}: {num_significant_control}")
    print(f"Number of matching loops with p-value < {p_value}: {num_significant_matching}")


    # Apply FDR correction to p-values
    results_target_df['fdr'] = multipletests(results_target_df['pvalue'], method='fdr_bh', alpha=p_value)[1]
    results_control_df['fdr'] = multipletests(results_control_df['pvalue'], method='fdr_bh', alpha=p_value)[1]
    results_matching_df['fdr'] = multipletests(results_matching_df['pvalue'], method='fdr_bh', alpha=p_value)[1]

    num_significant_target_fdr = (results_target_df['fdr'] < p_value).sum()
    num_significant_control_fdr = (results_control_df['fdr'] < p_value).sum()
    num_significant_matching_fdr = (results_matching_df['fdr'] < p_value).sum()
    print(f"Number of unique target loops with FDR < {p_value}: {num_significant_target_fdr}")
    print(f"Number of unique control loops with FDR < {p_value}: {num_significant_control_fdr}")
    print(f"Number of matching loops with FDR < {p_value}: {num_significant_matching_fdr}")
def main(args=None):
    args = parse_args()
    detect_differential_loops(
        args.target_loop_file,
        args.target_matrix,
        args.control_loop_file,
        args.control_matrix,
        args.threads,
        args.p_value,
        args.window_size)
