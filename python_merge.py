import os
import re

import pandas as pd
def merge(file_path_a,file_path_b,results_path,ARRANGE):
    # Load data
    source_path_a=file_path_a
    source_path_b=file_path_b
    filename_a=os.path.basename(source_path_a)
    filename_b=os.path.basename(source_path_b)
    df_a = pd.read_csv(source_path_a)
    df_b = pd.read_csv(source_path_b)

    # Define dictionary
    pc_group_dict = {
        r"301\.[^,]+?184": "EPA-PC",
        r"\[327\.[^,]+?184": "DHA-PC",
        r"301\..+?327.+?184": "EPA-PC;DHA-PC"

    }

    # Determine table type
    if "-PC" in filename_b:
        active_pattern_group = pc_group_dict
        value_threshold = 43.9909
        match_mode = "pc"

    elif "-PE" in filename_b:

        value_threshold = 2.0146
        match_mode = "pe"

    results = []

    # Core logic
    for i in range(len(df_a)):
        g_left = df_a.iloc[i]['precursor']
        e_left = df_a.iloc[i]['RT']

        for j in range(len(df_b)):
            g_right = df_b.iloc[j]['precursor']
            e_right = df_b.iloc[j]['RT']

            # Step 1: G-value difference test
            g_diff_abs = abs(g_left - g_right)
            #2.0146
            #43.9909
            if not (value_threshold-0.005 <= g_diff_abs <= value_threshold+0.005):
                continue


            # Step 2: E-value difference test
            e_diff_abs = abs(e_left - e_right)
            if 0 <= e_diff_abs <= 0.2:
                results.append({
                    'a_index': i,
                    'b_index': j,
                    'e_diff': e_diff_abs ,
                    'a_row': i + 2,
                'b_row': j + 2
                })

    results.sort(key=lambda x:x["e_diff"])
    final_results=[]
    used_a=set()
    used_b=set()
    for res in results:
        a_idx=res["a_index"]
        b_idx=res["b_index"]
        if a_idx not in used_a and b_idx not in used_b:
            final_results.append((res["a_row"],res["b_row"]))
            used_a.add(a_idx)
            used_b.add(b_idx)


    results=final_results
    # Output results

    print(f"\nTotal matches found: {len(results)} groups")


    if results:
        result_rows = []
        for idx, (left_row, right_row) in enumerate(results):
            row_a = df_a.iloc[left_row - 2].copy()
            row_b = df_b.iloc[right_row - 2].copy()

            result_rows.append(row_a)
            result_rows.append(row_b)

        # Create results DataFrame
        result_df = pd.DataFrame(result_rows)

        # Enter Phase 2
        result_df.index = range(2, 2 + len(result_df))

        result_df['sort_key'] = result_df[ARRANGE].copy()
        result_df['Resource'] = ''
        result_df['Assignment'] = ''


        index_list = result_df.index.tolist()
        for i in range(0, len(index_list), 2):
            if i + 1 < len(index_list):
                even_idx = index_list[i]
                odd_idx = index_list[i + 1]

                result_df.loc[odd_idx, 'sort_key'] = result_df.loc[even_idx, 'sort_key']
                result_df.loc[even_idx, 'Resource'] = filename_a
                result_df.loc[odd_idx, 'Resource'] = filename_b

                key_concat = str(result_df.loc[even_idx, 'keyion_observe']) + str(result_df.loc[odd_idx, 'keyion_observe'])

                belong_result = None
                if match_mode == "pc":
                    for pattern, belong in active_pattern_group.items():
                        if re.search(pattern, key_concat):
                            belong_result = belong
                            break
                else :
                    complex_pe_scenarios = re.findall(r"(\d{3})\.", key_concat)
                    groups=complex_pe_scenarios

                    # Convert last two capture groups to int and validate absolute difference is 141
                    last1 = int(groups[-1])
                    last2 = int(groups[-2])
                    if abs(last1 - last2) != 141:
                        exit(1)

                    # Extract remaining capture groups (excluding last two) into a set
                    other_groups = set(groups[:-2])

                    # Validate set elements
                    DHA = "327" in other_groups
                    EPA = "301" in other_groups
                    has_both = DHA and EPA
                    if DHA and not has_both:
                        belong_result="DHA-PE"
                    elif EPA and not has_both:
                        belong_result="EPA-PE"
                    elif has_both:
                        belong_result="EPA-PE;DHA-PE"
                result_df.loc[even_idx, 'Assignment'] = belong_result or 'unclassified'
                result_df.loc[odd_idx, 'Assignment'] = belong_result or 'unclassified'

        # Sort DataFrame and drop specified columns
        result_df['sort_key'] = (result_df['sort_key'] + result_df.groupby('sort_key').cumcount() * 0.000001).round(6)

        result_df = result_df.sort_values(by=['sort_key', 'Resource'], ascending=[True, True])
        result_df = result_df.drop(columns=['sort_key'])

        os.makedirs(os.path.dirname(results_path), exist_ok=True)
        result_df.to_excel(results_path, index=False)
        print("OK")

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    paths = {}
    with open('python_merge.conf', 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                paths[key] = value.strip()

    results_path = os.path.join(base_dir, "merge_results", f"{paths.get('results_path')}.xlsx")
    file_path_a = os.path.join(base_dir, f"{paths.get('file_path_a')}.csv")
    file_path_b = os.path.join(base_dir, f"{paths.get('file_path_b')}.csv")
    ARRANGE = paths.get('ARRANGE')
    merge(file_path_a, file_path_b, results_path,ARRANGE)


