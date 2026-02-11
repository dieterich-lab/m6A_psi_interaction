import gzip
import argparse
import sys
import os

def process_bedmethyl(input_file, mod_code, min_cov, min_pct, output_file):
    print(f"Processing {input_file} extracting mod '{mod_code}' (cov>={min_cov}, pct>={min_pct})...")
    
    count_pass = 0
    count_total = 0
    
    with gzip.open(input_file, 'rt') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#') or line.startswith('track'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
                
            count_total += 1
            
            # Columns (0-based):
            # 0: chrom
            # 1: start
            # 2: end
            # 3: mod
            # 4: score
            # 5: strand
            # ...
            # 9: coverage
            # 10: percent_modified
            
            try:
                f_mod = parts[3]
                f_cov = int(parts[9])
                f_pct = float(parts[10]) # percent often 0-100 or 0.0-100.0
                
                if f_mod == mod_code:
                    if f_cov >= min_cov and f_pct >= min_pct:
                        # Output BED6
                        # Use percent as score? Or original score?
                        # Original score (col 4 in file) is effectively integer.
                        # We keep original columns 0-5.
                        
                        fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[10]}\t{parts[5]}\n")
                        count_pass += 1
                        
            except ValueError:
                continue
                
    print(f"Finished. Total lines: {count_total}. Kept: {count_pass}. Output: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract sites from bedmethyl file")
    parser.add_argument("-i", "--input", required=True, help="Input bedmethyl.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument("-m", "--mod", required=True, help="Modification code (e.g. 'a' or '17802')")
    parser.add_argument("--min_cov", type=int, default=10, help="Minimum coverage")
    parser.add_argument("--min_pct", type=float, default=5.0, help="Minimum modification percentage")
    
    args = parser.parse_args()
    
    process_bedmethyl(args.input, args.mod, args.min_cov, args.min_pct, args.output)

if __name__ == "__main__":
    main()
