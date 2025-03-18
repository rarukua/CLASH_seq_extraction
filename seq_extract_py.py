#!/usr/bin/env python
# coding: utf-8

import json
import pandas as pd
import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract piRNA/mRNA matching information from SAM and dictionary files."
    )
    # Accept one or more SAM file paths.
    parser.add_argument(
        "--sam_file",
        type=str,
        nargs='+',
        required=True,
        help="Path(s) to the SAM CSV file(s) (e.g., CLASH_sam1_cleaned.csv CLASH_sam2_cleaned.csv, etc.)"
    )
    # Accept one or more JSON dictionary files.
    parser.add_argument(
        "--dic_file",
        type=str,
        nargs='+',
        required=True,
        help="Path(s) to the JSON file(s) with transcript-to-sequence mapping (e.g., dic_enst2seq_UTR3.json dic_enst2seq_CDS.json)"
    )
    # The piRNA file (default value can be overridden).
    parser.add_argument(
        "--piRNA_file",
        type=str,
        default="piRNA_gold.csv",
        help="Path to the piRNA CSV file (default: piRNA_gold.csv)"
    )
    return parser.parse_args()

def find_substring_with_one_mismatch(substring, full_string, max_mismatch=1):
    """
    Return True if 'substring' occurs in 'full_string' allowing up to 'max_mismatch' mismatches.
    Otherwise, return False.
    """
    len_sub = len(substring)
    len_full = len(full_string)
    
    # If substring is longer than full_string, no match is possible
    if len_sub > len_full:
        return False
    
    # Check every possible starting position in full_string
    for start_idx in range(len_full - len_sub + 1):
        mismatch_count = 0
        # Compare characters
        for i in range(len_sub):
            if full_string[start_idx + i] != substring[i]:
                mismatch_count += 1
                if mismatch_count > max_mismatch:
                    # Too many mismatches; break and check next start position
                    break
        
        # If we never exceeded max_mismatch, we have a valid match
        if mismatch_count <= max_mismatch:
            return True
    
    # If no valid start index was found, return False
    return False


def find_approx_matches(short_seq, long_seq,max_mismatch=2):
    """
    Slide 'short_seq' over 'long_seq' and return all windows 
    where the mismatch count is exactly 0 or exactly 2.
    
    Returns a list of dicts like:
      [
        {
          "start_idx": <start position in long_seq>,
          "window_seq": <substring of long_seq>,
          "mismatch_count": 0 or 2
        },
        ...
      ]
    """
    results = []
    len_short = len(short_seq)
    len_long = len(long_seq)
    
    # If short_seq is longer than long_seq, no valid alignment
    if len_short > len_long:
        return results
    half_length = len_long // 2
    
    # Slide 'short_seq' across 'long_seq'
    for start_idx in range(len_long - len_short + 1):
        window_seq = long_seq[start_idx : start_idx + len_short]
        mismatch_count = count_mismatches(short_seq, window_seq)
        
        # Keep only windows with exactly 0 or 2 mismatches
        if mismatch_count <= max_mismatch and start_idx <= half_length:
            results.append({
                "start_idx": start_idx,
                "window_seq": window_seq,
                "mismatch_count": mismatch_count
            })
    
    return results

def find_approx_matches_tail_to_head(short_seq, long_seq, max_mismatch=2):
    """
    Slide 'short_seq' along 'long_seq' from tail to head,
    returning all windows that have <= max_mismatch mismatches.
    
    Unlike a typical left-to-right search, we start from the end of long_seq
    and move toward index 0. We do NOT reverse short_seq in any way.
    
    Args:
        short_seq (str): The unmatched piRNA portion (tail).
        long_seq (str): The piRNA or mRNA we’re searching in.
        max_mismatch (int): The maximum number of mismatches allowed (defaults to 2).
    
    Returns:
        A list of dicts, each with:
            {
              "start_idx": <start position in long_seq>,
              "window_seq": <substring of long_seq>,
              "mismatch_count": <number of mismatches>,
            }
    """
    results = []
    len_short = len(short_seq)
    len_long = len(long_seq)
    
    # If the piRNA tail (short_seq) is longer than long_seq, no valid alignment
    if len_short > len_long:
        return results
    half_length = len_long // 2

    # We'll start from the last possible alignment index
    # and move toward 0, effectively searching tail→head.
    for start_idx in range(len_long - len_short, -1, -1):
        window_seq = long_seq[start_idx : start_idx + len_short]
        
        # Count mismatches
        mismatch_count = sum(a != b for a, b in zip(short_seq, window_seq))
        
        if mismatch_count <= max_mismatch and start_idx >= half_length:
            results.append({
                "start_idx": start_idx,
                "window_seq": window_seq,
                "mismatch_count": mismatch_count
            })
    
    return results


def count_mismatches(seq1, seq2):
    """Count the number of positions at which seq1 and seq2 differ."""
    return sum(a != b for a, b in zip(seq1, seq2))


def process_sam_file(sam_file, dic_enst2seq, piRNA_seq):
    """
    Process a single SAM file and return two DataFrames:
    one for piRNA matches and one for mRNA matches.
    """
    # Read the SAM file.
    sam = pd.read_csv(sam_file)
    # (Optional) Extract columns if needed (assumed: V3, V10, V4).
    # Select rows with long sequences (length > 50)
    unique_long_strings = set(sam['V10'][sam['V10'].str.len() > 50])
    long_string_rows = sam[sam['V10'].isin(unique_long_strings)]
    
    # ----- First Matching: piRNA_mRNA_matches from long_string_rows -----
    piRNA_mRNA_matches = {}
    grouped = long_string_rows.groupby('V3')
    for target_string_name, group in grouped:
        if target_string_name in dic_enst2seq:
            # Unpack the tuple from the dictionary.
            long_string, binding_region, binding_start, binding_end = dic_enst2seq[target_string_name]
            for _, row in group.iterrows():
                target_string = row['V10']
                target_start = int(row['V4'])  # Get the start position hint from column V4
                if not (binding_start <= target_start <= binding_end):
                    continue
                for i in range(len(target_string)):
                    substring = target_string[i:]
                    if find_substring_with_one_mismatch(substring, long_string, max_mismatch=1):
                        matched = substring
                        unmatched_part = target_string[:i]
                        match_with_unmatched = unmatched_part + substring
                        binding = ""
                        if binding_start <= target_start <= binding_end:
                            binding = binding_region  # e.g., "UTR5:"
                        if target_string_name not in piRNA_mRNA_matches:
                            piRNA_mRNA_matches[target_string_name] = []
                        piRNA_mRNA_matches[target_string_name].append({
                            "target_string": target_string,
                            "matched": substring,
                            "unmatched": unmatched_part,
                            "combined": match_with_unmatched,
                            "binding_region": binding
                        })
                        break

    # Filter matches: keep only entries with matched part ≥ 11 and unmatched part ≥ 6.
    piRNA_mRNA_matches_filtered = [
        (gene, [entry for entry in entries 
                if len(entry['matched']) >= 30 and len(entry['unmatched']) >= 10])
        for gene, entries in piRNA_mRNA_matches.items()
        if any(len(entry['matched']) >= 30 and len(entry['unmatched']) >= 10 for entry in entries)
    ]
    
    # Also, filter matches for second matching: entries with matched part < 30.
    filtered_matches = [
        (gene, [entry for entry in entries if len(entry['matched']) < len(entry['unmatched'])])
        for gene, entries in piRNA_mRNA_matches.items()
        if any(len(entry['matched']) < len(entry['unmatched']) for entry in entries)
    ]
    mRNA_piRNA_transcripts = [gene for gene, _ in filtered_matches]
    filtered_sam = long_string_rows[long_string_rows['V3'].isin(mRNA_piRNA_transcripts)]
    
    # ----- Second Matching: mRNA_piRNA_matches from filtered_sam -----
    mRNA_piRNA_matches = {}
    grouped = filtered_sam.groupby('V3')
    for target_string_name, group in grouped:
        if target_string_name in dic_enst2seq:
            long_string, binding_region, binding_start, binding_end = dic_enst2seq[target_string_name]
            for _, row in group.iterrows():
                target_string = row['V10']
                target_start = int(row['V4'])
                if not (binding_start <= target_start <= binding_end):
                    continue
                for i in range(len(target_string)):
                    substring = target_string[:-i-1]
                    if find_substring_with_one_mismatch(substring, long_string, max_mismatch=1):
                        matched = substring
                        unmatched_part = target_string[-i-1:]
                        match_with_unmatched = unmatched_part + matched
                        binding = ""
                        if binding_start <= target_start <= binding_end:
                            binding = binding_region
                        if target_string_name not in mRNA_piRNA_matches:
                            mRNA_piRNA_matches[target_string_name] = []
                        mRNA_piRNA_matches[target_string_name].append({
                            "target_string": target_string,
                            "matched": substring,
                            "unmatched": unmatched_part,
                            "combined": match_with_unmatched,
                            "binding_region": binding
                        })
                        break
    mRNA_piRNA_matches_filtered = [
        (gene, [entry for entry in entries 
                if len(entry['matched']) >= 30 and len(entry['unmatched']) >= 10])
        for gene, entries in mRNA_piRNA_matches.items()
        if any(len(entry['matched']) >= 30 and len(entry['unmatched']) >= 10 for entry in entries)
    ]
    
    # ----- Process piRNA Matches -----
    results_piRNA = []
    for transcript_id, match_list in piRNA_mRNA_matches_filtered:
        for match_info in match_list:
            matched_seq = match_info["matched"]
            unmatched_seq = match_info["unmatched"]
            binding_region=match_info['binding_region']
            for idx, piRNA_row in piRNA_seq.iterrows():
                piRNA_name = piRNA_row["piRNA"]
                piRNA_full_seq = piRNA_row["seq"]
                approx_matches = find_approx_matches_tail_to_head(unmatched_seq, piRNA_full_seq)
                if approx_matches:
                    for m in approx_matches:
                        results_piRNA.append({
                            "transcript_id": transcript_id,
                            "target_string": match_info["target_string"],
                            "mRNA_part": matched_seq,
                            "unmatched_seq": unmatched_seq,
                            "piRNA_name": piRNA_name,
                            "piRNA_sequence": piRNA_full_seq,
                            "matched_part": m["window_seq"],
                            "start_index_in_piRNA": m["start_idx"],
                            "mismatch_count": m["mismatch_count"],
                            "binding_region": binding_region
                        })
    df_piRNA = pd.DataFrame(results_piRNA)
    df_piRNA['piRNA_mRNA_order'] = 1

    # ----- Process mRNA Matches -----
    results_mRNA = []
    for transcript_id, match_list in mRNA_piRNA_matches_filtered:
        for match_info in match_list:
            matched_seq = match_info["matched"]
            unmatched_seq = match_info["unmatched"]
            binding_region = match_info["binding_region"]
            for idx, piRNA_row in piRNA_seq.iterrows():
                piRNA_name = piRNA_row["piRNA"]
                piRNA_full_seq = piRNA_row["seq"]
                approx_matches = find_approx_matches(unmatched_seq, piRNA_full_seq)
                if approx_matches:
                    for m in approx_matches:
                        results_mRNA.append({
                            "transcript_id": transcript_id,
                            "target_string": match_info["target_string"],
                            "mRNA_part": matched_seq,
                            "unmatched_seq": unmatched_seq,
                            "piRNA_name": piRNA_name,
                            "piRNA_sequence": piRNA_full_seq,
                            "matched_part": m["window_seq"],
                            "start_index_in_piRNA": m["start_idx"],
                            "mismatch_count": m["mismatch_count"],
                            "binding_region": binding_region
                        })
    df_mRNA = pd.DataFrame(results_mRNA)
    df_mRNA['piRNA_mRNA_order'] = 0

    return df_piRNA, df_mRNA

def main():
    args = parse_arguments()
    
    # Merge multiple JSON files into a single dictionary.
    dic_enst2seq = {}
    for dic_file in args.dic_file:
        with open(dic_file, 'r') as f:
            temp_dic = json.load(f)
        dic_enst2seq.update(temp_dic)
    
    # Load the piRNA file once.
    piRNA_seq = pd.read_csv(args.piRNA_file)
    dic_basename = "_".join([os.path.splitext(os.path.basename(f))[0] for f in args.dic_file])
    
    # Process each SAM file individually.
    for sam_file in args.sam_file:
        sam_basename = os.path.splitext(os.path.basename(sam_file))[0]
        print(f"Processing {sam_file} ...")
        df_piRNA, df_mRNA = process_sam_file(sam_file, dic_enst2seq, piRNA_seq)
        
        # Create dynamic output filenames using the input file's base name.
        piRNA_output_filename = f"./second_batch/result/piRNA_mRNA_pairing_results_{sam_basename}_{dic_basename}.csv"
        mRNA_output_filename = f"./second_batch/result/mRNA_piRNA_pairing_results_{sam_basename}_{dic_basename}.csv"
        
        df_piRNA.to_csv(piRNA_output_filename, index=False)
        df_mRNA.to_csv(mRNA_output_filename, index=False)
        print(f"Saved output to:\n  {piRNA_output_filename}\n  {mRNA_output_filename}")

if __name__ == "__main__":
    main()
