# Adapted from ISMapper
# Copyright (c) 2014, Jane Hawkey, Kathryn Holt

import os
import csv
import logging
from Bio import SeqIO

import bisect

def parse_genbank(gbk_file, qualifier='product'):
    """
    Parse a genbank file to extract all CDS/tRNA/rRNA features 
    to annotate adjacent genes.
    """
    features = []
    for record in SeqIO.parse(gbk_file, 'genbank'):
        for feature in record.features:
            if feature.type in ['CDS', 'tRNA', 'rRNA']:
                locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
                description = feature.qualifiers.get(qualifier, ['unknown'])[0]
                features.append({
                    'start': int(feature.location.start),
                    'end': int(feature.location.end),
                    'strand': feature.location.strand,
                    'locus_tag': locus_tag,
                    'description': description
                })
    # Sort features by start position for fast bisect search
    features.sort(key=lambda x: x['start'])
    return features

def find_closest_gene(pos, features):
    """Find the closest gene feature to a given position."""
    if pos == -1 or not features:
        return 'unknown', 'unknown'
        
    starts = [f['start'] for f in features]
    idx = bisect.bisect_left(starts, pos)
    
    closest_dist = float('inf')
    closest_feature = None
    
    # Check a narrow window around the insertion to handle overlapping or long preceding genes
    for i in range(max(0, idx - 10), min(len(features), idx + 10)):
        f = features[i]
        # Distance to start or end
        dist = min(abs(pos - f['start']), abs(pos - f['end']))
        if f['start'] <= pos <= f['end']:
            dist = 0 # Interrupting the gene
            
        if dist < closest_dist:
            closest_dist = dist
            closest_feature = f
            
    if closest_feature:
        return closest_feature['locus_tag'], closest_feature['description']
    return 'unknown', 'unknown'

def parse_coverage(cov_file):
    """Parse a bedgraph coverage file into a memory dictionary."""
    cov_data = {}
    if not os.path.exists(cov_file):
        return cov_data
        
    with open(cov_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4: continue
            
            chrom = parts[0]
            c_start = int(parts[1])
            c_end = int(parts[2])
            c_depth = int(parts[3])
            
            if chrom not in cov_data:
                cov_data[chrom] = []
            cov_data[chrom].append((c_start, c_end, c_depth))
    return cov_data

def calculate_depth_stats(cov_data, chrom, start, end):
    """Calculate mean depth, median depth, and IQR string from pre-parsed coverage data."""
    depths = []
    if chrom not in cov_data:
        return 0, 0, "0-0"
        
    for c_start, c_end, c_depth in cov_data[chrom]:
        # Optimization: cov_data is sorted by start, so we can break early
        if c_start >= end:
            break
        if c_end <= start:
            continue
            
        # Find overlap
        overlap_start = max(start, c_start)
        overlap_end = min(end, c_end)
        if overlap_start < overlap_end:
            # Add depth for each base in overlap
            depths.extend([c_depth] * (overlap_end - overlap_start))
                 
    if not depths:
        return 0, 0, "0-0"
        
    mean_val = sum(depths) / len(depths)
    
    # Calculate median and IQR (linear interpolation)
    sorted_depths = sorted(depths)
    n = len(sorted_depths)
    
    def get_val(p):
        idx = (n - 1) * p
        idx_f = int(idx)
        idx_c = idx_f + 1 if idx_f + 1 < n else idx_f
        weight = idx - idx_f
        return sorted_depths[idx_f] * (1.0 - weight) + sorted_depths[idx_c] * weight

    median_val = get_val(0.5)
    q1_val = get_val(0.25)
    q3_val = get_val(0.75)
    
    iqr_str = f"{round(q1_val, 2)}-{round(q3_val, 2)}"
    
    return mean_val, median_val, iqr_str

def get_peaks(merged_bed, cov_file):
    """Load peaks from BED and calculate depths."""
    peaks = []
    if not os.path.exists(merged_bed): return peaks
    
    # Parse coverage into memory ONCE
    cov_data = parse_coverage(cov_file)
    
    with open(merged_bed, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0]
                start, end = int(parts[1]), int(parts[2])
                mean_cov, median_cov, iqr_str = calculate_depth_stats(cov_data, chrom, start, end)
                peaks.append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'mean': mean_cov,
                    'median': median_cov,
                    'iqr': iqr_str,
                    'chimera': False,
                    'paired': False
                })
    return peaks

from .mapping import run_command

def check_overlap(p1, p2):
    start = max(p1['start'], p2['start'])
    end = min(p1['end'], p2['end'])
    return max(0, end - start)

def is_partial_overlap(p1, p2, max_overlap=100):
    ov = check_overlap(p1, p2)
    if ov == 0:
        return False
    len1 = p1['end'] - p1['start']
    len2 = p2['end'] - p2['start']
    # Not a full containment
    if ov >= len1 - 5 or ov >= len2 - 5:
        return False
    return ov < max_overlap

def is_full_overlap(p1, p2):
    ov = check_overlap(p1, p2)
    if ov == 0:
        return False
    len1 = p1['end'] - p1['start']
    len2 = p2['end'] - p2['start']
    return ov >= min(len1, len2) * 0.9

def get_flanking_features(x, y, features):
    """
    Find flanking features for insertion at [x, y].
    If x is -1 or y is -1, handle appropriately.
    """
    left_feat = None
    right_feat = None
    
    if x == -1 and y == -1:
        return None, None
        
    if not features:
        return None, None
        
    starts = [f['start'] for f in features]
    
    # 1. Resolve left feature (starts <= x)
    if x != -1:
        idx_l = bisect.bisect_left(starts, x)
        for i in range(idx_l, -1, -1):
            if i < len(features):
                f = features[i]
                if f['start'] <= x:
                    left_feat = f
                    break
        if not left_feat:
            left_feat = features[0]
            
    # 2. Resolve right feature (starts >= y)
    if y != -1:
        idx_r = bisect.bisect_left(starts, y)
        if idx_r - 1 >= 0 and features[idx_r - 1]['start'] <= y <= features[idx_r - 1]['end']:
            right_feat = features[idx_r - 1]
        elif idx_r < len(features):
            right_feat = features[idx_r]
        else:
            right_feat = features[-1]
            
    # Fallback to other feature if one is not resolved
    if not left_feat and right_feat:
        left_feat = right_feat
    if not right_feat and left_feat:
        right_feat = left_feat
        
    return left_feat, right_feat

def scan_known_is_positions(targets_fasta, ref_fasta, threads, tmp_dir):
    """
    Scan the reference genome for the full sequence of the IS elements.
    Returns a list of dictionaries: [{'chr': chrom, 'start': start, 'end': end, 'name': name, 'strand': strand}, ...]
    """
    known_is = []
    if not os.path.exists(targets_fasta):
        logging.warning(f"Targets file {targets_fasta} not found. Skipping reference scan.")
        return known_is

    full_seqs = {}
    for record in SeqIO.parse(targets_fasta, "fasta"):
        # Header format: >ST16:IS1R_IS1:FULL or similar
        if ":FULL" in record.id:
            full_seqs[record.id] = record

    if not full_seqs:
        logging.warning("No :FULL target sequences found in targets.fasta. Skipping reference scan.")
        return known_is

    os.makedirs(tmp_dir, exist_ok=True)
    tmp_targets_fa = os.path.join(tmp_dir, "tmp_targets_full.fasta")
    with open(tmp_targets_fa, "w") as f:
        for record in full_seqs.values():
            SeqIO.write(record, f, "fasta")

    # Resolve ref_fasta to FASTA format if it is a GenBank file
    mapping_ref = ref_fasta
    if ref_fasta.endswith('.gb') or ref_fasta.endswith('.gbk'):
        ref_base_name = os.path.basename(ref_fasta).rsplit('.', 1)[0]
        converted_fa = os.path.join(tmp_dir, ref_base_name + '.fasta')
        if os.path.exists(converted_fa):
            mapping_ref = converted_fa

    from .mapping import bwa_index
    bwa_cmd = bwa_index(mapping_ref)
    
    sam_file = os.path.join(tmp_dir, "targets_mapped.sam")
    # bwa mem with -a option to output all alignments (repeats copies)
    cmd = f"{bwa_cmd} mem -a -t {threads} {mapping_ref} {tmp_targets_fa} > {sam_file}"
    logging.info(f"Mapping target sequences to reference genome: {cmd}")
    run_command(cmd, shell=True)

    import re
    def get_align_len(cigar):
        return sum(int(val) for val, op in re.findall(r'(\d+)([MIDNX=H])', cigar) if op in 'MDN=X')

    if os.path.exists(sam_file):
        with open(sam_file, "r") as f:
            for line in f:
                if line.startswith("@"): continue
                parts = line.strip().split("\t")
                if len(parts) < 6: continue
                flag = int(parts[1])
                if flag & 4: continue # Unmapped
                
                name = parts[0]
                chrom = parts[2]
                pos = int(parts[3])
                cigar = parts[5]
                
                align_len = get_align_len(cigar)
                start = pos
                end = pos + align_len
                strand = '-' if (flag & 16) else '+'
                
                known_is.append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'strand': strand
                })
        logging.info(f"Identified {len(known_is)} known IS element copies on the reference genome.")
    return known_is

def parse_bed_hits(left_merged, right_merged, left_cov, right_cov, features, is_length, flank_len, targets_fasta, ref_fasta, threads):
    """
    Parse the bedtools output to find paired hits, chimeras, and single flanks.
    Uses reference target guide scan and coordinate overlap matching rules.
    """
    left_peaks = get_peaks(left_merged, left_cov)
    right_peaks = get_peaks(right_merged, right_cov)
    
    for p in left_peaks:
        p['paired'] = False
    for p in right_peaks:
        p['paired'] = False
        
    hits = []
    
    # 1. Full-Overlap Artifact Filtering (Chimeras)
    for l_peak in left_peaks:
        for r_peak in right_peaks:
            if l_peak['chr'] != r_peak['chr']: continue
            
            # Check overlap
            if is_full_overlap(l_peak, r_peak):
                ratio = max(l_peak['mean'], r_peak['mean']) / max(0.001, min(l_peak['mean'], r_peak['mean']))
                if ratio > 5.0:
                    # Chimera
                    if l_peak['mean'] < r_peak['mean']:
                        l_peak['chimera'] = True
                    else:
                        r_peak['chimera'] = True
                            
    # Filter chimeras out of the active pools
    left_peaks = [p for p in left_peaks if not p['chimera']]
    right_peaks = [p for p in right_peaks if not p['chimera']]

    # Scan the reference genome for known IS positions
    tmp_dir = os.path.dirname(left_cov)
    known_is = scan_known_is_positions(targets_fasta, ref_fasta, threads, tmp_dir)

    total_known = len(known_is)
    detected_known = 0

    # ==========================================
    # STAGE 1: Resolve Known IS Elements
    # ==========================================
    for is_el in known_is:
        chrom = is_el['chr']
        start = is_el['start']
        end = is_el['end']
        strand = is_el['strand']
        
        l_peak = None
        r_peak = None
        
        if strand == '+':
            # HEAD is at start, TAIL is at end
            for lp in left_peaks:
                if lp['paired']: continue
                if lp['chr'] != chrom: continue
                if not (lp['end'] < start - flank_len or lp['start'] > start + 100):
                    l_peak = lp
                    break
            for rp in right_peaks:
                if rp['paired']: continue
                if rp['chr'] != chrom: continue
                if not (rp['end'] < end - 100 or rp['start'] > end + flank_len):
                    r_peak = rp
                    break
        else: # strand == '-'
            # TAIL is at start, HEAD is at end
            for rp in right_peaks:
                if rp['paired']: continue
                if rp['chr'] != chrom: continue
                if not (rp['end'] < start - flank_len or rp['start'] > start + 100):
                    l_peak = rp
                    break
            for lp in left_peaks:
                if lp['paired']: continue
                if lp['chr'] != chrom: continue
                if not (lp['end'] < end - 100 or lp['start'] > end + flank_len):
                    r_peak = lp
                    break
                    
        if l_peak and r_peak:
            l_peak['paired'] = True
            r_peak['paired'] = True
            detected_known += 1
            lp_hit = l_peak if strand == '+' else r_peak
            rp_hit = r_peak if strand == '+' else l_peak
            
            hits.append({
                'type': 'Known Pair',
                'l_peak': lp_hit,
                'r_peak': rp_hit,
                'orientation': strand
            })

    logging.info(f"Number of Known IS on reference sequences: {total_known}")
    logging.info(f"Number of detected (paired) Known IS: {detected_known}")

    # ==========================================
    # STAGE 2: Resolve Novel Tandems and Single Insertions
    # ==========================================
    left_clean = [p for p in left_peaks if not p['paired']]
    right_clean = [p for p in right_peaks if not p['paired']]

    # 2.1: Novel Same-Direction Tandem (++ or --)
    # central peaks: r1 (TAIL) and l2 (HEAD) fully overlap.
    # outer peaks: l1 (HEAD) and r2 (TAIL) partially overlap with this central pair.
    for rp1 in list(right_clean):
        if rp1['paired']: continue
        for lp2 in list(left_clean):
            if lp2['paired']: continue
            if rp1['chr'] != lp2['chr']: continue
            
            if is_full_overlap(rp1, lp2):
                # Search for l1 (HEAD) upstream
                lp1 = None
                for p in left_clean:
                    if p['paired']: continue
                    if p['chr'] != rp1['chr']: continue
                    if is_partial_overlap(p, rp1) or is_partial_overlap(p, lp2):
                        if p['start'] < rp1['start']:
                            lp1 = p
                            break
                            
                # Search for r2 (TAIL) downstream
                rp2 = None
                for p in right_clean:
                    if p['paired']: continue
                    if p['chr'] != rp1['chr']: continue
                    if is_partial_overlap(p, rp1) or is_partial_overlap(p, lp2):
                        if p['start'] > rp1['start']:
                            rp2 = p
                            break
                            
                if lp1 and rp2:
                    lp1['paired'] = True
                    rp1['paired'] = True
                    lp2['paired'] = True
                    rp2['paired'] = True
                    
                    hits.append({'type': 'Novel Pair', 'l_peak': lp1, 'r_peak': rp1})
                    hits.append({'type': 'Novel Pair', 'l_peak': lp2, 'r_peak': rp2})
                    hits.append({'type': 'Tandem Gap (same direction)', 'l_peak': lp2, 'r_peak': rp1})
                    break # Break inner loop for rp1

    # Update active clean pools
    left_clean = [p for p in left_clean if not p['paired']]
    right_clean = [p for p in right_clean if not p['paired']]

    # 2.2: Novel Opposite-Direction Tandem (+-)
    # One central TAIL peak (rpN) partially overlapping with two distinct HEAD peaks (lp1, lp2)
    for rpN in list(right_clean):
        if rpN['paired']: continue
        overlapping_heads = []
        for lp in left_clean:
            if lp['paired']: continue
            if lp['chr'] == rpN['chr'] and is_partial_overlap(lp, rpN):
                overlapping_heads.append(lp)
                
        if len(overlapping_heads) == 2:
            lp1, lp2 = overlapping_heads
            lp1['paired'] = True
            lp2['paired'] = True
            rpN['paired'] = True
            
            hits.append({'type': 'Tandem Pair (+-)', 'l_peak': lp1, 'r_peak': rpN})
            hits.append({'type': 'Tandem Gap (opposite - HEAD)', 'l_peak': lp2, 'r_peak': None})

    # Update active clean pools
    left_clean = [p for p in left_clean if not p['paired']]
    right_clean = [p for p in right_clean if not p['paired']]

    # 2.3: Novel Opposite-Direction Tandem (-+)
    # One central HEAD peak (lpN) partially overlapping with two distinct TAIL peaks (rp1, rp2)
    for lpN in list(left_clean):
        if lpN['paired']: continue
        overlapping_tails = []
        for rp in right_clean:
            if rp['paired']: continue
            if rp['chr'] == lpN['chr'] and is_partial_overlap(rp, lpN):
                overlapping_tails.append(rp)
                
        if len(overlapping_tails) == 2:
            rp1, rp2 = overlapping_tails
            rp1['paired'] = True
            rp2['paired'] = True
            lpN['paired'] = True
            
            hits.append({'type': 'Tandem Pair (-+)', 'l_peak': lpN, 'r_peak': rp1})
            hits.append({'type': 'Tandem Gap (opposite - TAIL)', 'l_peak': None, 'r_peak': rp2})

    # Update active clean pools
    left_clean = [p for p in left_clean if not p['paired']]
    right_clean = [p for p in right_clean if not p['paired']]

    # 2.4: Novel Single Insertion
    # Remaining l_peak (HEAD) and r_peak (TAIL) on same contig that partially overlap (overlap < 100 bp)
    # or are separated by <= 100 bp.
    for lp in list(left_clean):
        if lp['paired']: continue
        closest_dist = float('inf')
        closest_r = None
        
        for rp in right_clean:
            if rp['paired']: continue
            if lp['chr'] != rp['chr']: continue
            
            ov = check_overlap(lp, rp)
            if ov > 0:
                dist = 0 if ov < 100 else float('inf')
            else:
                if lp['end'] < rp['start']:
                    dist = rp['start'] - lp['end']
                elif rp['end'] < lp['start']:
                    dist = lp['start'] - rp['end']
                else:
                    dist = 0
            
            if dist > 100:
                if lp['end'] < rp['start'] and dist > 100:
                    break
                continue
                
            if dist < closest_dist:
                closest_dist = dist
                closest_r = rp
                
        if closest_r:
            lp['paired'] = True
            closest_r['paired'] = True
            if is_full_overlap(lp, closest_r):
                hits.append({'type': 'Ambiguous Full Overlap', 'l_peak': lp, 'r_peak': closest_r})
            else:
                ov = check_overlap(lp, closest_r)
                if ov > 0:
                    hits.append({'type': 'Novel Pair (TSD)', 'l_peak': lp, 'r_peak': closest_r})
                else:
                    hits.append({'type': 'Novel Pair', 'l_peak': lp, 'r_peak': closest_r})

    # Update active clean pools for singletons
    left_clean = [p for p in left_peaks if not p['paired']]
    right_clean = [p for p in right_peaks if not p['paired']]

    # ==========================================
    # STAGE 3: Resolve Remaining Single Flanks
    # ==========================================
    for lp in left_clean:
        hits.append({'type': 'HEAD-only', 'l_peak': lp, 'r_peak': None})
    for rp in right_clean:
        hits.append({'type': 'TAIL-only', 'l_peak': None, 'r_peak': rp})

    # 4. Resolve Orientation and Formatting
    final_hits = []
    for h in hits:
        if h['type'] == 'Ambiguous Full Overlap':
            h['type'] = 'Off-Target Amplicon (Noise)'
        orientation = '?'
        left_pos_val = 'N/A'
        right_pos_val = 'N/A'
        left_depth_median = 'N/A'
        left_depth_iqr = 'N/A'
        right_depth_median = 'N/A'
        right_depth_iqr = 'N/A'
        
        lp = h['l_peak']
        rp = h['r_peak']
        
        if 'orientation' in h:
            orientation = h['orientation']
        elif lp and rp:
            # Check overlap or standard orientation to resolve orientation
            ov = check_overlap(lp, rp)
            if ov > 0:
                orientation = '+' if lp['start'] < rp['start'] else '-'
            else:
                orientation = '+' if lp['end'] < rp['start'] else '-'
                
        # Resolve Coordinates and Depths
        if lp and rp:
            ov = check_overlap(lp, rp)
            if ov > 0:
                if lp['start'] < rp['start']:
                    left_pos_val = f"{lp['start']}-{lp['end']}"
                    left_depth_median = round(lp['median'], 2)
                    left_depth_iqr = lp['iqr']
                    right_pos_val = f"{rp['start']}-{rp['end']}"
                    right_depth_median = round(rp['median'], 2)
                    right_depth_iqr = rp['iqr']
                else:
                    left_pos_val = f"{rp['start']}-{rp['end']}"
                    left_depth_median = round(rp['median'], 2)
                    left_depth_iqr = rp['iqr']
                    right_pos_val = f"{lp['start']}-{lp['end']}"
                    right_depth_median = round(lp['median'], 2)
                    right_depth_iqr = lp['iqr']
            else:
                if lp['end'] < rp['start']: # lp is left
                    left_pos_val = f"{lp['start']}-{lp['end']}"
                    left_depth_median = round(lp['median'], 2)
                    left_depth_iqr = lp['iqr']
                    right_pos_val = f"{rp['start']}-{rp['end']}"
                    right_depth_median = round(rp['median'], 2)
                    right_depth_iqr = rp['iqr']
                else: # rp is left
                    left_pos_val = f"{rp['start']}-{rp['end']}"
                    left_depth_median = round(rp['median'], 2)
                    left_depth_iqr = rp['iqr']
                    right_pos_val = f"{lp['start']}-{lp['end']}"
                    right_depth_median = round(lp['median'], 2)
                    right_depth_iqr = lp['iqr']
                    
        elif lp: # HEAD-only
            orientation = '?'
            left_pos_val = f"{lp['start']}-{lp['end']}"
            left_depth_median = round(lp['median'], 2)
            left_depth_iqr = lp['iqr']
            
        elif rp: # TAIL-only
            orientation = '?'
            right_pos_val = f"{rp['start']}-{rp['end']}"
            right_depth_median = round(rp['median'], 2)
            right_depth_iqr = rp['iqr']
            
        # Get coordinates and gap
        x_val = 'N/A'
        y_val = 'N/A'
        gap_val = 'N/A'
        
        if left_pos_val != 'N/A' and right_pos_val != 'N/A':
            l_start, l_end = map(int, left_pos_val.split('-'))
            r_start, r_end = map(int, right_pos_val.split('-'))
            if l_end < r_start:
                x_val = l_end
                y_val = r_start
                gap_val = y_val - x_val
            elif r_end < l_start:
                x_val = r_end
                y_val = l_start
                gap_val = y_val - x_val
            else: # overlapping
                overlap_start = max(l_start, r_start)
                overlap_end = min(l_end, r_end)
                x_val = overlap_start
                y_val = overlap_end
                gap_val = -(y_val - x_val)
        elif left_pos_val != 'N/A':
            l_start, l_end = map(int, left_pos_val.split('-'))
            x_val = l_start
            y_val = l_end
        elif right_pos_val != 'N/A':
            r_start, r_end = map(int, right_pos_val.split('-'))
            x_val = r_start
            y_val = r_end
            
        # Get Flanking Genes Details
        l_coord = int(left_pos_val.split('-')[1]) if left_pos_val != 'N/A' else -1
        r_coord = int(right_pos_val.split('-')[0]) if right_pos_val != 'N/A' else -1
        
        left_feat, right_feat = get_flanking_features(l_coord, r_coord, features)
        
        left_gene = 'N/A'
        left_description = 'N/A'
        left_strand = 'N/A'
        left_distance = 0
        left_interrupted = 'False'
        
        if left_feat:
            left_gene = left_feat['locus_tag']
            left_description = left_feat['description']
            left_strand = str(left_feat['strand'])
            if l_coord != -1:
                if left_feat['start'] <= l_coord <= left_feat['end']:
                    left_distance = 0
                    left_interrupted = 'True'
                else:
                    left_distance = abs(left_feat['end'] - l_coord)
                    
        right_gene = 'N/A'
        right_description = 'N/A'
        right_strand = 'N/A'
        right_distance = 0
        right_interrupted = 'False'
        
        if right_feat:
            right_gene = right_feat['locus_tag']
            right_description = right_feat['description']
            right_strand = str(right_feat['strand'])
            if r_coord != -1:
                if right_feat['start'] <= r_coord <= right_feat['end']:
                    right_distance = 0
                    right_interrupted = 'True'
                else:
                    right_distance = abs(right_feat['start'] - r_coord)
                    
        gene_interruption = 'True' if (left_interrupted == 'True' or right_interrupted == 'True') else 'False'
        
        ismap_orientation = '?'
        if orientation == '+':
            ismap_orientation = 'F'
        elif orientation == '-':
            ismap_orientation = 'R'
            
        call_val = h['type']
        if call_val == 'Known Pair':
            call_val = 'known'
        elif call_val == 'Novel Pair':
            call_val = 'novel'
        elif call_val == 'Novel Pair (TSD)':
            call_val = 'novel (TSD)'
            
        contig = lp['chr'] if lp else rp['chr']
        
        final_hits.append({
            'type': h['type'],
            'contig': contig,
            'orientation': ismap_orientation,
            'x': x_val,
            'y': y_val,
            'gap': gap_val,
            'call': call_val,
            'left_pos': left_pos_val,
            'right_pos': right_pos_val,
            'left_depth_median': left_depth_median,
            'left_depth_iqr': left_depth_iqr,
            'right_depth_median': right_depth_median,
            'right_depth_iqr': right_depth_iqr,
            'left_gene': left_gene,
            'left_description': left_description,
            'left_strand': left_strand,
            'left_distance': left_distance,
            'right_gene': right_gene,
            'right_description': right_description,
            'right_strand': right_strand,
            'right_distance': right_distance,
            'gene_interruption': gene_interruption
        })
        
    # Sort hits by contig and position
    def get_sort_key(hit):
        pos = 0
        if hit['x'] != 'N/A':
            pos = int(hit['x'])
        elif hit['y'] != 'N/A':
            pos = int(hit['y'])
        return (hit['contig'], pos)
        
    final_hits.sort(key=get_sort_key)
    return final_hits

def create_typing_output(left_merged, right_merged, left_cov, right_cov, ref_fasta, out_file, is_length=4000, flank_len=300, targets_fasta="config/targets.fasta", threads=1):
    """
    Generate final summary table.
    """
    logging.info(f"Generating mapping report: {out_file}")
    
    # Check if reference is genbank
    features = []
    if ref_fasta.endswith('.gbk') or ref_fasta.endswith('.gb'):
        features = parse_genbank(ref_fasta)
    else:
        logging.warning("Reference is not GenBank. Gene annotation will be skipped.")
        
    hits = parse_bed_hits(left_merged, right_merged, left_cov, right_cov, features, is_length, flank_len, targets_fasta, ref_fasta, threads)
    
    # Separate into main table and unpaired singletons
    main_hits = []
    unpaired_hits = []
    
    main_types = ['Known Pair', 'Novel Pair', 'Novel Pair (TSD)', 'Tandem Pair (+-)', 'Tandem Pair (-+)', 'Tandem Gap (same direction)', 'Tandem Gap (opposite - HEAD)', 'Tandem Gap (opposite - TAIL)']
    
    main_idx = 1
    unpaired_idx = 1
    for h in hits:
        if h['type'] in main_types:
            h['region'] = f"region_{main_idx}"
            main_idx += 1
            main_hits.append(h)
        else:
            h['region'] = f"unpaired_{unpaired_idx}"
            unpaired_idx += 1
            unpaired_hits.append(h)
            
    fieldnames = [
        'region', 'contig', 'orientation', 'x', 'y', 'gap', 'call',
        'left_pos', 'right_pos',
        'left_depth_median', 'left_depth_iqr', 
        'right_depth_median', 'right_depth_iqr', 
        'left_gene', 'left_description', 'left_strand', 'left_distance',
        'right_gene', 'right_description', 'right_strand', 'right_distance',
        'gene_interruption'
    ]
                  
    # Write main hits to table.tsv
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for hit in main_hits:
            row = {k: hit[k] for k in fieldnames if k in hit}
            writer.writerow(row)
            
    logging.info(f"Report successfully written to {out_file}.")
    
    # Write unpaired/noise hits to unpaired.tsv (omit left_gene & right_gene information)
    unpaired_fieldnames = [
        'region', 'contig', 'orientation', 'x', 'y', 'gap', 'call',
        'left_pos', 'right_pos',
        'left_depth_median', 'left_depth_iqr', 
        'right_depth_median', 'right_depth_iqr'
    ]
    out_unpaired = out_file.replace('_table.tsv', '_unpaired.tsv')
    with open(out_unpaired, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=unpaired_fieldnames, delimiter='\t')
        writer.writeheader()
        for hit in unpaired_hits:
            row = {k: hit[k] for k in unpaired_fieldnames if k in hit}
            writer.writerow(row)
            
    logging.info(f"Unpaired singletons successfully written to {out_unpaired}.")
