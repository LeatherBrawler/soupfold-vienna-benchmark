import subprocess
import re
import time
import pandas as pd
import traceback
import numpy as np

def run_multifold(sequences):
    start_time = time.time()
    multifold_proc = subprocess.run(
        ["RNAmultifold", "-d0"], 
        input=sequences, 
        capture_output=True, 
        text=True
    )
    end_time = time.time()
    runtime = end_time - start_time
    
    match = re.search(r"^([\.()&]+)\s+\(\s*([-+]?\d*\.\d+|\d+)\s*\)", multifold_proc.stdout, re.MULTILINE)
    if match:
        vienna_structure = match.group(1)
        vienna_mfe = float(match.group(2))
    else:
        vienna_structure = None
        vienna_mfe = None

    return vienna_structure, vienna_mfe, runtime


def run_soupfold(distinct_seqs, m_total):    
    cmd = ["./strand_soup.exe", "--hetero", str(m_total)] + distinct_seqs
    start_time = time.time()
    
    try:
        soup_proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        error_msg = f"CRASH: {e}\nSTDERR:\n{e.stderr}\nSTDOUT:\n{e.stdout}"
        raise RuntimeError(error_msg)

    end_time = time.time()
    runtime = end_time - start_time
    output = soup_proc.stdout
    # print(output)
    possibilities = []
    blocks = output.split("Starting point:")[1:]
    
    for block in blocks:
        # print("block:\n", block)
        mfe_match = re.search(r"E_total=([-+]?\d*\.?\d+)\s*kcal/mol", block)
        # print(mfe_match)
        soup_mfe = float(mfe_match.group(1)) if mfe_match else None
        
        order_match = re.search(r"Ordered list of sequences:\s*([\d\s]+)", block)
        # print(order_match)
        ordered_list = []
        if order_match:
            ordered_list = [int(idx) for idx in order_match.group(1).strip().split()]

        seq_match = re.search(r"^([A-Z&]+)$", block, re.MULTILINE)
        sequences = seq_match.group(1) if seq_match else None
            
        struct_match = re.search(r"^([\.()&]+)$", block, re.MULTILINE)
        soup_structure = struct_match.group(1) if struct_match else None
        
        possibilities.append({
            "soup_mfe": soup_mfe,
            "ordered_list": ordered_list,
            "soup_structure": soup_structure,
            "sequences": sequences
        })
    
    return possibilities, runtime, output

def create_test_cases(k_range = 1):
    # create different test cases, 
    # first randomize number of triplets, draw a graph of differences in mfe with increasing number of triplets, keeping m and k constant
    # then randomize the triplet choices, each with random k repeats, with higher ranges, eg k in (1, 10), (6, 15), (10, 20) 
    # then test on increasing m number of strands in the structure, 
    m_range = np.arange(2, 21, 3)
    k_ranges = [np.arange(1, 10), np.arange(6, 15), np.arange(10, 20)] # np.arange(20, 50)]
    triplets = ["CUG", "CAG", "GAU", "CGG", "UAG", "CCG"]
    number_of_triplets = np.arange(2, 12, 3)
    test_cases = []
    for m in m_range:
        for k_range in k_ranges:
            for num_triplets in number_of_triplets:
                chosen_triplets = np.random.choice(triplets, size=num_triplets, replace=True)
                distinct_seqs = []
                for triplet in chosen_triplets:
                    k = np.random.choice(k_range)
                    distinct_seqs.append(str(triplet) * int(k))
                test_cases.append((distinct_seqs, m))
    return test_cases

def reverse_sequence(seq):
    seq = seq[::-1]
    rev = ""
    for i in range(len(seq)):
        c = seq[i]
        if c == ')':
            rev += '('
        elif c == '(':
            rev += ')'
        else:
            rev += c
    return rev 

def find_first_pairing(seq):
    for i in range(len(seq)):
        if seq[i] == '(':
            return i

def shift_and_compare(seq1, seq2):
    n1 = find_first_pairing(seq1)
    n2 = find_first_pairing(seq2)
    n_diff = n1 - n2
    if n_diff == 0: # never reach this case if wrapped in compare_structure below
        return seq1 == seq2
    if n_diff > 0: # n1 > n2, add n_diff dots to n2 at start and remove n_diff dots at the end (just do it)
        new_seq2 = "." * n_diff + seq2[:-n_diff]
        # print(f"Shifted seq2 (Vienna): {new_seq2}")
        return seq1 == new_seq2
    else:
        new_seq1 = "." * (-n_diff) + seq1[:n_diff]
        # print(f"Shifted seq1(Soupfold): {new_seq1}")
        return seq2 == new_seq1

def fast_count_base_pairs(seq, seq1):
    """seq is nucleotide sequence, seq1 is dot bracket notation"""
    stack = []
    dict1 = {"GC": 0, "AU": 0, "GU": 0, "CG": 0, "UA": 0, "UG": 0, "Other": 0}
    for i in range(len(seq1)):
        if seq1[i] == '(':
            stack.append(i)
        elif seq1[i] == ')':
            if stack:
                j = stack.pop()
                bp = seq[j] + seq[i]
                if bp in dict1:
                    dict1[bp] += 1
                else:                    
                    dict1["Other"] += 1
    final_counts = {
        "GC": dict1["GC"] + dict1["CG"],
        "AU": dict1["AU"] + dict1["UA"],
        "GU": dict1["GU"] + dict1["UG"],
        "Other": dict1["Other"]
    }
    return final_counts

def compare_structure(seq1, seq2):
    # 0 is different, 1 is same structure, 2 is same with reverse, 3 is same with shifting, 4 is same with shifting and reverse
    if seq1 == seq2:
        return 1
    rev = reverse_sequence(seq2)
    if seq1 == rev:
        return 2
    if shift_and_compare(seq1, seq2):
        return 3
    # print("Going in shift for rev")
    if shift_and_compare(seq1, rev):
        return 4
    return 0

def main():
    # triplets = ["CUG", "CAG", "GAU", "CGG", "UAG", "CCG"]
    test_cases_eg = [
        # (distinct_sequences_list, m_total)
        (["CAGCAGCAG", "CGGCGG", "CAGCAGCAGCAG", "CUGCUGCUGCUGCUGCUG", "GAUGAUGAUGAUGAU"], 11),
        (["CAGCAGCAG", "CGGCGG", "CAGCAGCAGCAG"], 11),
        (["CAGCAGCAG", "CGGCGG", "CAGCAGCAGCAG", "CUGCUGCUGCUGCUGCUG", "GAUGAUGAUGAUGAU"], 5),
        (["CAGCAGCAG", "CGGCGG", "CAGCAGCAGCAG"], 5),
    ]
    
    test_cases = create_test_cases()

    results = []
    with open("Results/Output/error.txt", "w") as f:
        f.write(f"Starting new test..............\n")
        f.close()


    test_cases = test_cases_eg if not test_cases else test_cases
    num_big_difference = 0
    for idx, (distinct_seqs, m_total) in enumerate(test_cases):
        with open("Results/Output/error.txt", "a") as f:
            f.write(f"Parameters: m={m_total}, seqs={distinct_seqs}\n")
            f.close()

        print(f"\n--- Running Test {idx+1}: m={m_total}, seqs={distinct_seqs} ---")
        
        try:
            # 1. Run SoupFold
            possibilities, soup_time, raw_output = run_soupfold(distinct_seqs, m_total)
            
            if not possibilities:
                raise ValueError("SoupFold did not return any structures.")
            
            print(f"Found {len(possibilities)} possibilities. Running RNAmultifold for each...")
            
            # 2. Loop through every possibility found by SoupFold
            for p_idx, poss in enumerate(possibilities):
                soup_mfe = poss["soup_mfe"]
                ordered_indices = poss["ordered_list"]
                soup_struct = poss["soup_structure"]
                sequences = poss["sequences"]
                if not ordered_indices:
                    raise ValueError(f"Possibility {p_idx+1} missing ordered list.")
                
                vienna_struct, vienna_mfe, vienna_time = run_multifold(sequences)
                
                mfe_diff = soup_mfe - vienna_mfe if (soup_mfe is not None and vienna_mfe is not None) else None
                match2 = False
                # print(type(vienna_struct))
                match = compare_structure(soup_struct, vienna_struct)
                is_struct_match = (not (not match)) # !!match
                match2 = is_struct_match

                match_type = [False, "Same structure", "Same with inverse", "Same with shift", "Same with inverse shift"]
                print("==========================================================================================")
                print(f"  [{p_idx+1}/{len(possibilities)}] Soup: {soup_mfe} | Vienna: {vienna_mfe} | Struct Match? {match_type[match]} | mfe_diff: {mfe_diff:.2f}")
                if not is_struct_match:
                    with open("Results/Output/error.txt", "a") as f:
                        f.write(f"========== Difference in structure ==========\n")
                        f.write(f" Soup struct, Vienna struct, with sequence:\n")
                        f.write(f"{soup_struct}\n")
                        f.write(f"{vienna_struct}\n")
                        f.write(f"{sequences}\n")
                        dict1 = fast_count_base_pairs(sequences, soup_struct)
                        dict2 = fast_count_base_pairs(sequences, vienna_struct)
                        f.write(f"Base pair counts in Soup structure: {dict1}\n")
                        f.write(f"Base pair counts in Vienna structure: {dict2}\n")

                        if dict1 == dict2:
                            print("Same number of base pairs of each type")
                            match2 = True
                            f.write("Same number of base pairs of each type\n")
                num_big_difference += 1 if not match2 else 0
                results.append({
                    "Test_ID": idx + 1,
                    "Possibility_ID": p_idx + 1,           # Tracks which permutation this is
                    "Total_Possibilities": len(possibilities),
                    "m": m_total,
                    "Distinct_Seqs": ",".join(distinct_seqs),
                    "Used_Sequence_Order": " ".join(map(str, ordered_indices)),
                    "Full_Constructed_Seq": sequences,
                    "Soup_MFE": soup_mfe,
                    "Vienna_MFE": vienna_mfe,
                    "MFE_Diff": mfe_diff,
                    "Struct_Match": is_struct_match,
                    "Soup_Time(s)": soup_time,             # Total time taken by SoupFold for all possibilities 
                    "Vienna_Time(s)": vienna_time,          # Time taken for this specific RNAmultifold run
                    "Structure_Soup": soup_struct,
                    "Structure_Vienna": vienna_struct if match != 1 else 0, # meaning actually matching structure
                    "Base_Pair_Counts_Soup": dict1 if not is_struct_match else 0,
                    "Base_Pair_Counts_Vienna": dict2 if not is_struct_match else 0,
                    "Same Base_Pair_Counts": dict1 == dict2 if not is_struct_match else 0,
                    "Struct_Match_Type": match_type[match]
                })
            
        except Exception as e:
            print(f"FAILED. Writing to error.txt...")
            with open("Results/Output/error.txt", "a") as f:
                f.write(f"========== ERROR IN TEST {idx+1} ==========\n")
                f.write(f"Parameters: m={m_total}, seqs={distinct_seqs}\n")
                f.write(f"Exception Type: {type(e).__name__}\n")
                f.write(f"Details:\n{str(e)}\n")
                f.write(traceback.format_exc())
                f.write("===========================================\n\n")
            
            results.append({
                "Test_ID": idx + 1,
                "m": m_total,
                "Distinct_Seqs": ",".join(distinct_seqs),
                "Error": "Crashed - See error.txt"
            })

    df = pd.DataFrame(results)
    csv_filename = "Results/heterogeneous_benchmark_results.csv"
    df.to_csv(csv_filename, index=False)
    # print(f"\nBenchmarking complete. Data saved to {csv_filename}")
    print(f"Number of big differences: {num_big_difference} out of {len(results)} possibilities.")


if __name__ == "__main__":
    main()
