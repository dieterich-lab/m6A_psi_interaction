import argparse
import pandas as pd
import numpy as np
from pybedtools import BedTool
import os
import sys

def detect_bed_chrom_format(bed_file):
    """
    Detects if the BED file uses 'chr' prefix.
    Returns True if 'chr' is present, False otherwise.
    """
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            parts = line.split('\t')
            if len(parts) > 0:
                return parts[0].startswith('chr')
    return False

def parse_genepred(genepred_file, use_chr_prefix=True):
    """
    Parses a GenePred file and returns a dictionary of transcripts.
    use_chr_prefix: ensure chrom names match this preference.
    """
    transcripts = {}
    print(f"Reading GenePred file: {genepred_file} (Target chrom format: {'chr' if use_chr_prefix else 'no-chr'})")
    
    # GenePred often has a header line or bin column. We'll try to detect it.
    # UCSC genePred usually: bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds ...
    
    with open(genepred_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            
            # Flexible parsing: assume standard genePred if 10+ columns
            # We need: name(1), chrom(2), strand(3), txStart(4), txEnd(5), cdsStart(6), cdsEnd(7), exonStarts(9), exonEnds(10)
            # If bin column is present (often first column numeric), indices are shifted by 1.
            # Let's check if the first column is likely 'bin' (integer).
             
            has_bin = False
            try:
                int(parts[0])
                # If name is in col 1, it's string.
                # If bin is col 0, it's int.
                # However, name can be number-like.
                # Standard check: if columns are enough.
                if len(parts) >= 10:
                    # Check if col 2 is a chrom (chr...)
                    if parts[2].startswith('chr'):
                         has_bin = True
                    elif parts[1].startswith('chr'):
                         has_bin = False
            except ValueError:
                has_bin = False

            offset = 1 if has_bin else 0
            
            try:
                name = parts[0 + offset]
                chrom = parts[1 + offset]
                strand = parts[2 + offset]
                txStart = int(parts[3 + offset])
                txEnd = int(parts[4 + offset])
                cdsStart = int(parts[5 + offset])
                cdsEnd = int(parts[6 + offset])
                # exonCount = int(parts[7 + offset])
                exonStarts = [int(x) for x in parts[8 + offset].split(',') if x]
                exonEnds = [int(x) for x in parts[9 + offset].split(',') if x]
                name2 = parts[11 + offset] if len(parts) > 11 + offset else name 
                
                # Normalize chromosome
                if use_chr_prefix and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                elif not use_chr_prefix and chrom.startswith('chr'):
                    chrom = chrom[3:]
                    
            except (ValueError, IndexError):
                continue

            transcripts[name] = {
                'name': name,
                'name2': name2,
                'chrom': chrom,
                'strand': strand,
                'txStart': txStart,
                'txEnd': txEnd,
                'cdsStart': cdsStart,
                'cdsEnd': cdsEnd,
                'exonStarts': exonStarts,
                'exonEnds': exonEnds
            }
    
    print(f"Loaded {len(transcripts)} transcripts.")
    return transcripts

def get_transcript_features(tx, site_pos):
    """
    Calculates features (region lengths, relative position) for a site in a transcript.
    site_pos is 0-based genomic coordinate.
    """
    
    # 1. Calculate transcriptomic coordinates of exons
    # Flatten exons into a list of genomic positions included in transcript
    # and map genomic to transcriptomic.
    
    # However, we only need to map specific points: site_pos, cdsStart, cdsEnd.
    # And we need to know the total length of regions.
    
    chrom = tx['chrom']
    strand = tx['strand']
    exonStarts = tx['exonStarts']
    exonEnds = tx['exonEnds']
    cdsStart = tx['cdsStart']
    cdsEnd = tx['cdsEnd']
    
    # Check if site is in any exon
    site_in_exon = False
    for s, e in zip(exonStarts, exonEnds):
        if s <= site_pos < e:
            site_in_exon = True
            break
    if not site_in_exon:
        return None

    # Calculate lengths of boolean regions (UTR5, CDS, UTR3)
    # We iterate exons and sum up lengths of parts that fall into each region.
    
    len_5utr = 0
    len_cds = 0
    len_3utr = 0
    
    # Variables to track site position in transcriptomic terms
    dist_from_tx_start = 0 
    
    # For strand processing:
    # If '+', 5' is txStart (lowest coord). UTR5 is < cdsStart.
    # If '-', 5' is txEnd (highest coord). UTR5 is > cdsEnd.
    
    # Strategy: Calculate length of each region. Calculate distance of site from start of each region.
    
    # Genomic lengths in transcript order
    sorted_exons = list(zip(exonStarts, exonEnds)) # Always low to high
    
    # Calculate transcriptomic vectors for 5utr, cds, 3utr
    # We will accumulate lengths.
    
    # Iterate exons low to high
    
    # Define genomic regions for CDS
    # CDS is from cdsStart to cdsEnd.
    
    # Calculate per-exon contribution
    
    total_len_low_to_high = 0
    site_len_low_to_high = 0
    
    for estart, eend in sorted_exons:
        # Exon length
        elen = eend - estart
        
        # Site pos check
        if estart <= site_pos < eend:
            site_offset_in_exon = site_pos - estart
            site_len_low_to_high = total_len_low_to_high + site_offset_in_exon
        
        total_len_low_to_high += elen
        
        # Region contributions
        # Overlap with UTR/CDS regions based on genomic coordinates
        
        # 3 Segments genomic: (-inf, cdsStart), [cdsStart, cdsEnd), [cdsEnd, inf)
        
        # Segment 1: < cdsStart
        seg1_s = estart
        seg1_e = min(eend, cdsStart)
        if seg1_e > seg1_s:
            if strand == '+': len_5utr += (seg1_e - seg1_s)
            else: len_3utr += (seg1_e - seg1_s)
            
        # Segment 2: CDS
        seg2_s = max(estart, cdsStart)
        seg2_e = min(eend, cdsEnd)
        if seg2_e > seg2_s:
            len_cds += (seg2_e - seg2_s)
            
        # Segment 3: > cdsEnd
        seg3_s = max(estart, cdsEnd)
        seg3_e = eend
        if seg3_e > seg3_s:
            if strand == '+': len_3utr += (seg3_e - seg3_s)
            else: len_5utr += (seg3_e - seg3_s)

    if len_cds == 0:
        return None # Non-coding transcript or site in non-coding? Filter it.

    # Calculate site position relative to regions
    if strand == '+':
        site_tx_pos = site_len_low_to_high
        
        # Calculate transcriptomic start of CDS
        # Iterate again? Or simpler: 
        # CDS start is determined by sum of lengths < cdsStart
        cds_start_tx = 0
        for estart, eend in sorted_exons:
            seg_e = min(eend, cdsStart)
            if seg_e > estart:
                cds_start_tx += (seg_e - estart)
        
        cds_end_tx = cds_start_tx + len_cds
        
        dist_from_cds_start = site_tx_pos - cds_start_tx
        dist_from_cds_end = site_tx_pos - cds_end_tx
        
    else: # strand == '-'
        # Transcript starts at highest coordinate.
        # site_len_low_to_high is distance from genomic start (lowest coord).
        # genomic_len = total_len_low_to_high (after loop)
        
        total_transcript_len = total_len_low_to_high
        # site_tx_pos (distance from 5' end ie high genomic)
        # site (genomic) is estart + offset.
        # dist from high end = (total_len) - (dist_from_low_end + 1)? 
        # site_pos is 0-based.
        # let's map site_pos to 1-based index from low end?
        # site_len_low_to_high is 0-based index from low end.
        
        site_tx_pos = (total_transcript_len - 1) - site_len_low_to_high
        
        # CDS start is at cdsEnd (genomic).
        # Calculate length of regions > cdsEnd (genomic) to get 5'UTR length (transcriptomic)
        cds_start_tx = 0
        for estart, eend in sorted_exons:
            seg_s = max(estart, cdsEnd)
            if eend > seg_s:
                cds_start_tx += (eend - seg_s)
        
        cds_end_tx = cds_start_tx + len_cds
        
        dist_from_cds_start = site_tx_pos - cds_start_tx
        dist_from_cds_end = site_tx_pos - cds_end_tx

    # Determine relative location [0-3]
    rel_location = None
    
    if site_tx_pos < cds_start_tx: # 5' UTR
        if len_5utr > 0:
            rel_location = site_tx_pos / len_5utr
        else:
            rel_location = 0 # Should not happen if site is there
            
    elif site_tx_pos < cds_end_tx: # CDS
        if len_cds > 0:
            rel_location = 1.0 + (dist_from_cds_start / len_cds)
            
    else: # 3' UTR
        if len_3utr > 0:
            rel_location = 2.0 + (dist_from_cds_end / len_3utr)
        else:
            rel_location = 3.0 # At the very end?

    return {
        'rel_location': rel_location,
        'dist_from_cds_start': dist_from_cds_start,
        'dist_from_cds_end': dist_from_cds_end,
        'utr5_size': len_5utr,
        'cds_size': len_cds,
        'utr3_size': len_3utr
    }

def main():
    parser = argparse.ArgumentParser(description="Generate metagene coordinates from BED file and GenePred")
    parser.add_argument("-b", "--bed", required=True, help="Input BED file with sites")
    parser.add_argument("-g", "--genepred", required=True, help="Input GenePred annotation file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    
    args = parser.parse_args()
    
    # Detect BED format
    bed_has_chr = detect_bed_chrom_format(args.bed)
    print(f"Detected BED chromosome format: {'chr' if bed_has_chr else 'no-chr'}")

    # 1. Parse GenePred
    transcripts = parse_genepred(args.genepred, use_chr_prefix=bed_has_chr)
    
    # 2. Intersect BED with GenePred using pybedtools
    # We first make a temporary BED6 file from GenePred to allow intersection
    # Or simplified: BED12. 
    # But pybedtools from_dataframe is easier if we construct a dataframe
    
    print("Preparing annotation for intersection...")
    # Create bed12 string or file
    bed_lines = []
    for tx in transcripts.values():
        # bed12: chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
        # thickStart = cdsStart, thickEnd = cdsEnd
        
        blockSizes = []
        blockStarts = []
        
        # Verify order
        # exonStarts/Ends are sorted.
        e_starts = tx['exonStarts']
        e_ends = tx['exonEnds']
        
        for s, e in zip(e_starts, e_ends):
            blockSizes.append(str(e - s))
            blockStarts.append(str(s - tx['txStart']))
            
        bed_line = [
            tx['chrom'],
            str(tx['txStart']),
            str(tx['txEnd']),
            tx['name'],
            "0",
            tx['strand'],
            str(tx['cdsStart']),
            str(tx['cdsEnd']),
            "0",
            str(len(e_starts)),
            ",".join(blockSizes),
            ",".join(blockStarts)
        ]
        bed_lines.append("\t".join(bed_line))
    
    annot_bed = BedTool("\n".join(bed_lines), from_string=True)
    sites_bed = BedTool(args.bed)
    
    print("Intersecting sites with annotation...")
    # -wo: Write the original A and B entries plus the number of base pairs of overlap.
    # We want site info and transcript info.
    # sites_bed could be any format. genePred mapped to BED12 is named by transcript ID.
    # intersection: site mapping to full transcript span?
    # Yes, -wa -wb.
    
    intersection = sites_bed.intersect(annot_bed, wa=True, wb=True)
    
    results = []
    
    print("Processing intersections...")
    # Intersection fields depend on input BED cols.
    # If input is BED6, first 6 cols are site.
    # Then BED12 cols (12 cols).
    # We need to know where the site is.
    
    # We can iterate the intersection.
    # 'interval' object of pybedtools
    # interval.fields gives all fields
    
    for interval in intersection:
        # fields is a list of strings
        fields = interval.fields
        
        # Assume input BED has at least 3 cols: chr, start, end.
        site_chrom = fields[0]
        site_start = int(fields[1])
        site_end = int(fields[2]) # Half-open
        
        # Site position: usually the start for single nucleotide (or mid?)
        # For m6A sites, usually single base.
        site_pos = site_start
        
        # Annotation starts after site fields.
        # How many cols in input? pybedtools handles this.
        # interval[0]...interval[N-1] are site fields?
        # The fields from the 'B' file (annotation) start after 'A' fields.
        # But we don't know N easily unless we check.
        # The 'name' of the transcript is in column 4 of the BED12 (0-indexed).
        # So in the combined fields, look for the transcript name.
        
        # Actually, let's just match using the transcript name which we put in the bed12.
        # The overlap gives us the transcript name.
        
        # Finding the transcript part:
        # Since we generated the partial BED, we know the format.
        # BED12 has 12 columns.
        # The name is the 4th column of the BED12 part.
        # It's better to use `wb=True` and know that the last 12 columns are the annotation.
        
        annot_fields = fields[-12:]
        tx_name = annot_fields[3]
        
        if tx_name not in transcripts:
            continue
            
        tx = transcripts[tx_name]
        
        # Calculate features
        feats = get_transcript_features(tx, site_pos)
        
        if feats and feats['rel_location'] is not None:
            results.append({
                'chr': site_chrom,
                'coord': site_pos,
                'rel_location': feats['rel_location'],
                'dist_from_cds_start': feats['dist_from_cds_start'],
                'dist_from_cds_end': feats['dist_from_cds_end'],
                'utr5_size': feats['utr5_size'],
                'cds_size': feats['cds_size'],
                'utr3_size': feats['utr3_size'],
                'tx_id': tx['name'],
                'gene_name': tx.get('name2', tx['name'])
            })
            
    if not results:
        print("No valid sites mapped to transcripts found.")
        sys.exit(0)
        
    df = pd.DataFrame(results)
    
    print(f"Writing {len(df)} mapped sites to {args.output}")
    # Write TSV
    df.to_csv(args.output, sep='\t', index=False)
    
if __name__ == "__main__":
    main()
