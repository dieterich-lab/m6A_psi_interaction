import argparse
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def configure_matplotlib(dpi=1200, font_size=8):
    mpl.rcParams['figure.dpi'] = dpi
    mpl.rcParams['savefig.dpi'] = dpi
    mpl.rcParams['font.size'] = font_size
    mpl.rcParams['legend.fontsize'] = font_size - 2
    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size
    mpl.rcParams['xtick.major.size'] = 4
    mpl.rcParams['ytick.major.size'] = 4
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['font.family'] = 'Arial'

def get_longest_isoform_optimized(df):
    """
    Filters DataFrame to keep only the longest isoform for each unique genomic site (chr, coord).
    """
    # Calculate transcript length
    df['tx_len'] = df['utr5_size'] + df['cds_size'] + df['utr3_size']
    
    # Sort by tx_len descending
    # If multiple transcripts have same max length, determination is arbitrary (first one)
    df_sorted = df.sort_values('tx_len', ascending=False)
    
    # Drop duplicates by chr/coord, keeping first (longest)
    return df_sorted.drop_duplicates(subset=['chr', 'coord'])

def main():
    parser = argparse.ArgumentParser(description="Plot combined metagene density for m6A and pseudoU")
    parser.add_argument("--m6a_input", required=True, help="Input TSV for m6A")
    parser.add_argument("--psi_input", required=True, help="Input TSV for PseudoU")
    parser.add_argument("--output", required=True, help="Output plot file (e.g. .png)")
    parser.add_argument("--use_longest_isoform", action="store_true", help="Filter to longest isoform")
    
    args = parser.parse_args()
    
    configure_matplotlib()
    
    print(f"Reading m6A data: {args.m6a_input}")
    df_m6a = pd.read_csv(args.m6a_input, sep='\t')
    df_m6a['type'] = 'm6A'
    
    print(f"Reading PseudoU data: {args.psi_input}")
    df_psi = pd.read_csv(args.psi_input, sep='\t')
    df_psi['type'] = 'PseudoU'
    
    if args.use_longest_isoform:
        print("Filtering m6A to longest isoforms...")
        df_m6a = get_longest_isoform_optimized(df_m6a)
        print("Filtering PseudoU to longest isoforms...")
        df_psi = get_longest_isoform_optimized(df_psi)
    
    # Combine
    df_combined = pd.concat([df_m6a, df_psi])
    
    print("Plotting...")
    plt.figure(figsize=(6, 4))
    
    # Use distinct colors. m6A = Red, PseudoU = Blue
    # The user specifies "m6a project uses red and blue".
    # We assign Red to the primary modification (m6A) and Blue to the secondary (PseudoU)
    # mirroring the common usage in this project's volcano plots (pos=Red, neg=Blue).
    palette = {'m6A': 'tab:red', 'PseudoU': 'tab:blue'}
    
    # KDE Plot
    sns.kdeplot(data=df_combined, x='rel_location', hue='type', 
                common_norm=False, fill=False, linewidth=2,
                palette=palette)
    
    # Annotate regions
    plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=2, color='black', linestyle='--', alpha=0.5)
    
    # Add region labels (Reverting to the informative layout)
    # We need to get current ylim to position text
    current_ylim = plt.ylim()
    
    plt.text(0.5, current_ylim[1]*0.95, "5' UTR", ha='center')
    plt.text(1.5, current_ylim[1]*0.95, "CDS", ha='center')
    plt.text(2.5, current_ylim[1]*0.95, "3' UTR", ha='center')
    
    plt.xlim(0, 3)
    plt.xlabel('Metagene Coordinate')
    plt.ylabel('Density')
    plt.title('Metagene Distribution of Modifications')
    
    plt.tight_layout()
    plt.savefig(args.output)
    print(f"Saved plot to {args.output}")

if __name__ == "__main__":
    main()
