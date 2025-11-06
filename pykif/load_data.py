import os
import pandas as pd
from MSFileReader import MSFileReader
import numpy as np
import pickle

def load_corydalis_decumbens():
    # Define the base path
    Path = r'C:\work\msclassifier\data'

    # Define the Excel file path
    excel_file = os.path.join(Path, '夏天无151个成分表-朱金凤-20250711.xlsx')

    # Read the Excel file
    try:
        df = pd.read_excel(excel_file)
        df = df.dropna(subset=['数据'])
    except FileNotFoundError:
        print(f"Error: Excel file not found at {excel_file}")
        exit()

    # Initialize a list to hold all the data to be pickled
    all_spectra_data = []

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        # Extract information from the Excel row
        scanNum = int(row['二级图编号']) # Convert scanNum to an integer
        fileName_excel = row['数据'] + '.raw' # Assuming '数据' column contains the raw file name without extension
        classL = row['大类']
        subclassL = row['小类']
        filePath = os.path.join(Path, fileName_excel)

        try:
            raw_file = MSFileReader(filePath)
        except Exception as e:
            print(f"Error opening raw file {fileName_excel}: {e}")
            continue

        # Get label data (mass and intensity lists)
        try:
            label_data = raw_file.GetLabelData(scanNum)
            if label_data and len(label_data) > 0:
                a = label_data[0]
                mz_values = a[0]
                intensity_values = a[1]
            else:
                print(f"No label data found for scan number {scanNum} in {fileName_excel}")
                continue
        except Exception as e:
            print(f"Error getting label data for scan number {scanNum} in {fileName_excel}: {e}")
            continue

        # Prepare the data dictionary for the current spectrum
        current_spectrum_data = {
            'peaks': {
                'mz': np.array(mz_values, dtype=np.float32),
                'intensities': np.array(intensity_values, dtype=np.float32)
            },
            'metaData': [
                {'name': 'scanNum', 'value': scanNum},
                {'name': 'raw filename', 'value': fileName_excel}
            ],
            'classification_L:':[
                {'name': 'class', 'value': classL},
                {'name': 'sub_class', 'value': subclassL}
            ]
        }
        all_spectra_data.append(current_spectrum_data)
        print(index)
        raw_file.Close()
    # Define the output pickle file name for the consolidated data
    output_pickle_file = os.path.join(Path, "xtw.pickle")

    # Save all the collected data to a single pickle file
    try:
        with open(output_pickle_file, 'wb') as f:
            pickle.dump(all_spectra_data, f)
        print(f"\nSuccessfully saved all spectra data to {output_pickle_file}")
    except Exception as e:
        print(f"Error saving consolidated pickle file {output_pickle_file}: {e}")

# --- 新增函数: process_and_save_spectra_from_excel ---
def process_iris_domestica_1_data():
    """
    Reads spectrum scan information from 'Table S1-YWY-20220419.xlsx' sheet '终版'.
    It groups scans by their corresponding .raw files to minimize file I/O,
    reads spectrum data using MSFileReader, assigns fixed classification labels
    (flavonoid, isoflavone as per the provided snippet's behavior),
    and saves positive and negative ion mode data to separate pickle files:
    'Iris_domestica_1_P.pickle' and 'Iris_domestica_1_N.pickle'.
    """
    excel_path = r'C:\work\msclassifier\data\Table S1\Table S1-YWY-20220419.xlsx'
    raw_base_path = r'C:\work\msclassifier\data\Table S1'
    
    # New output pickle file names for positive and negative ions
    output_pickle_file_p = r'C:\work\msclassifier\data\Iris_domestica_1_positive.pickle'
    output_pickle_file_n = r'C:\work\msclassifier\data\Iris_domestica_1_negative.pickle'
    
    sheet_name = '终版'

    fixed_class = 'flavonoid'
    fixed_sub_class = 'isoflavone'

    print(f"--- Starting processing from Excel: {excel_path}, Sheet: {sheet_name} ---")

    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except FileNotFoundError:
        print(f"Error: Excel file not found at {excel_path}")
        return
    except KeyError:
        print(f"Error: Sheet '{sheet_name}' not found in {excel_path}")
        return
    except Exception as e:
        print(f"Error reading Excel file {excel_path}: {e}")
        return

    # Initialize separate lists for positive and negative spectra
    all_spectra_data_p = []
    all_spectra_data_n = []

    # Helper function to extract scan number from string like "1正07-11316"
    def extract_scan_number(scan_str):
        try:
            return int(str(scan_str).split('-')[-1])
        except (ValueError, IndexError):
            print(f"Warning: Could not extract scan number from '{scan_str}'. Skipping.")
            return None

    # Define mappings for prefixes to raw file paths and assigned raw filenames
    prefix_map = {
        '1正': {'raw_file_sub_path': os.path.join('P', '苷元-P.raw'), 'assigned_filename': 'P_ID_01.raw', 'ion_mode_type': 'positive'},
        '2正': {'raw_file_sub_path': os.path.join('P', 'ID_02.raw'), 'assigned_filename': 'P_ID_02.raw', 'ion_mode_type': 'positive'},
        '1负': {'raw_file_sub_path': os.path.join('N', 'ID_01.raw'), 'assigned_filename': 'N_ID_01.raw', 'ion_mode_type': 'negative'},
        '2负': {'raw_file_sub_path': os.path.join('N', 'ID_02.raw'), 'assigned_filename': 'N_ID_02.raw', 'ion_mode_type': 'negative'},
    }

    # Group scans by the raw file they belong to
    # {full_raw_file_path: [(scan_num, assigned_filename, ion_mode_type)]}
    grouped_scans_by_raw_file = {} 

    for index, row in df.iterrows():
        scan_string = row.get('scan') 
        if not pd.isna(scan_string) and isinstance(scan_string, str):
            scan_prefix = scan_string[:2]
            scan_num = extract_scan_number(scan_string)

            if scan_num is None:
                continue

            if scan_prefix in prefix_map:
                mapping = prefix_map[scan_prefix]
                raw_file_full_path = os.path.join(raw_base_path, mapping['raw_file_sub_path'])
                assigned_raw_filename = mapping['assigned_filename']
                ion_mode_type = mapping['ion_mode_type']

                if raw_file_full_path not in grouped_scans_by_raw_file:
                    grouped_scans_by_raw_file[raw_file_full_path] = []
                
                grouped_scans_by_raw_file[raw_file_full_path].append(
                    (scan_num, assigned_raw_filename, ion_mode_type)
                )
            # else: # No need to print warning for unrecognized prefixes unless debugging
            #     pass

    processed_count_p = 0
    processed_count_n = 0

    # Process each raw file
    for raw_file_full_path, scans_to_process in grouped_scans_by_raw_file.items():
        if not os.path.exists(raw_file_full_path):
            print(f"Warning: Raw file not found at {raw_file_full_path}. Skipping all scans for this file.")
            continue
        
        raw_file_reader = None
        try:
            raw_file_reader = MSFileReader(raw_file_full_path)
            print(f"Processing scans from {raw_file_full_path}...")

            for scan_num, assigned_raw_filename, ion_mode_type in scans_to_process:
                try:
                    label_data = raw_file_reader.GetLabelData(scan_num)
                    collision_energy = raw_file_reader.GetCollisionEnergyForScanNum(scan_num, 2)
                    scan_start_time = raw_file_reader.GetScanHeaderInfoForScanNum(scan_num)['StartTime']
                    tic = raw_file_reader.GetScanHeaderInfoForScanNum(scan_num)['TIC']
                    print(tic)
                    if label_data and len(label_data) > 0:
                        mz_values = label_data[0][0]
                        intensity_values = label_data[0][1]
                        
                        current_spectrum_data = {
                            'peaks': {
                                'mz': np.array(mz_values, dtype=np.float32),
                                'intensities': np.array(intensity_values, dtype=np.float32)
                            },
                            'metaData': [
                                {'name': 'scanNum', 'value': scan_num},
                                {'name': 'raw filename', 'value': assigned_raw_filename},
                                {'name': 'tic', 'value': tic},
                                {'name': 'collision_energy', 'value': collision_energy},
                                {'name': 'scan_start_time', 'value': scan_start_time}
                            ],
                            'classification_L:': [ # Hardcoded as per the provided snippet's behavior
                                {'name': 'class', 'value': fixed_class},
                                {'name': 'sub_class', 'value': fixed_sub_class},
                            ]
                        }
                        
                        if ion_mode_type == 'positive':
                            all_spectra_data_p.append(current_spectrum_data)
                            processed_count_p += 1
                        elif ion_mode_type == 'negative':
                            all_spectra_data_n.append(current_spectrum_data)
                            processed_count_n += 1
                    else:
                        print(f"No label data found for scan number {scan_num} in {raw_file_full_path}. Skipping.")
                except Exception as e:
                    print(f"Error processing scan {scan_num} from {raw_file_full_path}: {e}. Skipping.")
        except Exception as e:
            print(f"Error opening raw file {raw_file_full_path}: {e}")
        finally:
            if raw_file_reader:
                raw_file_reader.Close()

    # Save positive ion data to its own pickle file
    try:
        output_dir_p = os.path.dirname(output_pickle_file_p)
        if output_dir_p and not os.path.exists(output_dir_p):
            os.makedirs(output_dir_p, exist_ok=True)
        with open(output_pickle_file_p, 'wb') as f:
            pickle.dump(all_spectra_data_p, f)
        print(f"\n--- Successfully saved {processed_count_p} positive spectra data to {output_pickle_file_p} ---")
    except Exception as e:
        print(f"Error saving positive pickle file {output_pickle_file_p}: {e}")

    # Save negative ion data to its own pickle file
    try:
        output_dir_n = os.path.dirname(output_pickle_file_n)
        if output_dir_n and not os.path.exists(output_dir_n):
            os.makedirs(output_dir_n, exist_ok=True)
        with open(output_pickle_file_n, 'wb') as f:
            pickle.dump(all_spectra_data_n, f)
        print(f"--- Successfully saved {processed_count_n} negative spectra data to {output_pickle_file_n} ---")
    except Exception as e:
        print(f"Error saving negative pickle file {output_pickle_file_n}: {e}")

def process_iris_domestica_2_data():
    """
    Reads scan information from 'Table S2-YWY-20220419.xlsx' sheet '终版',
    filters rows based on 'Ion mode' (+ or -), reads corresponding spectrum data
    from specific .raw files, assigns fixed classification labels,
    and saves positive and negative ion mode data to separate pickle files:
    'Iris_domestica_2_P.pickle' and 'Iris_domestica_2_N.pickle'.
    """
    excel_path = r'C:\work\msclassifier\data\Table S2\Table S2-YWY-20220419.xlsx'
    raw_base_path = r'C:\work\msclassifier\data\Table S2'
    
    # New output pickle file names
    output_pickle_file_p = r'C:\work\msclassifier\data\Iris_domestica_2_positive.pickle'
    output_pickle_file_n = r'C:\work\msclassifier\data\Iris_domestica_2_negative.pickle'
    
    sheet_name = '终版'

    fixed_class = 'flavonoid'
    fixed_sub_class = 'isoflavone'

    print(f"--- Starting processing from Excel: {excel_path}, Sheet: {sheet_name} ---")

    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except FileNotFoundError:
        print(f"Error: Excel file not found at {excel_path}")
        return
    except KeyError:
        print(f"Error: Sheet '{sheet_name}' not found in {excel_path}")
        return
    except Exception as e:
        print(f"Error reading Excel file {excel_path}: {e}")
        return

    # Initialize separate lists for positive and negative spectra
    all_spectra_data_p = []
    all_spectra_data_n = []

    # Helper function to extract scan number from string like "#5087"
    def extract_scan_number(scan_str_with_hash):
        try:
            # Remove '#' and convert to integer
            return int(str(scan_str_with_hash).replace('#', ''))
        except (ValueError, AttributeError):
            print(f"Warning: Could not extract scan number from '{scan_str_with_hash}'. Skipping.")
            return None

    processed_count_p = 0
    processed_count_n = 0

    # Process positive ion mode rows
    positive_rows = df[df['Ion mode'].astype(str).str.endswith('+')].copy() # Ensure 'Ion mode' is string
    raw_file_p = os.path.join(raw_base_path, '20201204_TG_TG_5_P.raw')
    assigned_filename_p = '20201204_TG_TG_5_P.raw'

    if not os.path.exists(raw_file_p):
        print(f"Warning: Positive raw file not found at {raw_file_p}. Skipping positive ion mode processing.")
    else:
        raw_file_reader_p = None
        try:
            raw_file_reader_p = MSFileReader(raw_file_p)
            print(f"Processing positive ion mode from {raw_file_p}...")
            for index, row in positive_rows.iterrows():
                scan_num_str = row.get('二级') # Assuming '二级' is the column with scan numbers
                scan_num = extract_scan_number(scan_num_str)

                if scan_num is None:
                    continue
                try:
                    label_data = raw_file_reader_p.GetLabelData(scan_num)
                    collision_energy = raw_file_reader_p.GetCollisionEnergyForScanNum(scan_num, 2)
                    scan_start_time = raw_file_reader_p.GetScanHeaderInfoForScanNum(scan_num)['StartTime']
                    tic = raw_file_reader_p.GetScanHeaderInfoForScanNum(scan_num)['TIC']
                    if label_data and len(label_data) > 0:
                        mz_values = label_data[0][0]
                        intensity_values = label_data[0][1]
                        
                        current_spectrum_data = {
                            'peaks': {
                                'mz': np.array(mz_values, dtype=np.float32),
                                'intensities': np.array(intensity_values, dtype=np.float32)
                            },
                            'metaData': [
                                {'name': 'scanNum', 'value': scan_num},
                                {'name': 'raw filename', 'value': assigned_filename_p},
                                {'name': 'tic', 'value': tic},
                                {'name': 'collision_energy', 'value': collision_energy},
                                {'name': 'scan_start_time', 'value': scan_start_time}
                            ],
                            'classification_L:': [
                                {'name': 'class', 'value': fixed_class},
                                {'name': 'sub_class', 'value': fixed_sub_class}
                            ]
                        }
                        all_spectra_data_p.append(current_spectrum_data)
                        processed_count_p += 1
                        # print('P:', processed_count_p) # Keep for debugging if needed, remove for final
                    else:
                        print(f"No label data found for scan number {scan_num} in {raw_file_p}. Skipping.")
                except Exception as e:
                    print(f"Error processing positive scan {scan_num} from {raw_file_p}: {e}. Skipping.")
        except Exception as e:
            print(f"Error opening positive raw file {raw_file_p}: {e}")
        finally:
            if raw_file_reader_p:
                raw_file_reader_p.Close()


    # Process negative ion mode rows
    negative_rows = df[df['Ion mode'].astype(str).str.endswith('-')].copy() # Ensure 'Ion mode' is string
    raw_file_n = os.path.join(raw_base_path, '20201204_TG_TG_5_N.raw') # Note: user specified P folder for N raw file
    assigned_filename_n = '20201204_TG_TG_5_N.raw'

    if not os.path.exists(raw_file_n):
        print(f"Warning: Negative raw file not found at {raw_file_n}. Skipping negative ion mode processing.")
    else:
        raw_file_reader_n = None
        try:
            raw_file_reader_n = MSFileReader(raw_file_n)
            print(f"Processing negative ion mode from {raw_file_n}...")
            for index, row in negative_rows.iterrows():
                scan_num_str = row.get('二级')
                scan_num = extract_scan_number(scan_num_str)

                if scan_num is None:
                    continue

                try:
                    label_data = raw_file_reader_n.GetLabelData(scan_num)
                    collision_energy = raw_file_reader_n.GetCollisionEnergyForScanNum(scan_num, 2)
                    scan_start_time = raw_file_reader_n.GetScanHeaderInfoForScanNum(scan_num)['StartTime']
                    tic = raw_file_reader_n.GetScanHeaderInfoForScanNum(scan_num)['TIC']
                    if label_data and len(label_data) > 0:
                        mz_values = label_data[0][0]
                        intensity_values = label_data[0][1]
                        
                        current_spectrum_data = {
                            'peaks': {
                                'mz': np.array(mz_values, dtype=np.float32),
                                'intensities': np.array(intensity_values, dtype=np.float32)
                            },
                            'metaData': [
                                {'name': 'scanNum', 'value': scan_num},
                                {'name': 'raw filename', 'value': assigned_filename_n},
                                {'name': 'tic', 'value': tic},
                                {'name': 'collision_energy', 'value': collision_energy},
                                {'name': 'scan_start_time', 'value': scan_start_time}
                            ],
                            'classification_L:': [
                                {'name': 'class', 'value': fixed_class},
                                {'name': 'sub_class', 'value': fixed_sub_class}
                            ]
                        }
                        all_spectra_data_n.append(current_spectrum_data)
                        processed_count_n += 1
                        # print('N:', processed_count_n) # Keep for debugging if needed, remove for final
                    else:
                        print(f"No label data found for scan number {scan_num} in {raw_file_n}. Skipping.")
                except Exception as e:
                    print(f"Error processing negative scan {scan_num} from {raw_file_n}: {e}. Skipping.")
        except Exception as e:
            print(f"Error opening negative raw file {raw_file_n}: {e}")
        finally:
            if raw_file_reader_n:
                raw_file_reader_n.Close()

    # Save positive ion data to its own pickle file
    try:
        output_dir_p = os.path.dirname(output_pickle_file_p)
        if output_dir_p and not os.path.exists(output_dir_p):
            os.makedirs(output_dir_p, exist_ok=True)
        with open(output_pickle_file_p, 'wb') as f:
            pickle.dump(all_spectra_data_p, f)
        print(f"\n--- Successfully saved {processed_count_p} positive spectra data to {output_pickle_file_p} ---")
    except Exception as e:
        print(f"Error saving positive pickle file {output_pickle_file_p}: {e}")

    # Save negative ion data to its own pickle file
    try:
        output_dir_n = os.path.dirname(output_pickle_file_n)
        if output_dir_n and not os.path.exists(output_dir_n):
            os.makedirs(output_dir_n, exist_ok=True)
        with open(output_pickle_file_n, 'wb') as f:
            pickle.dump(all_spectra_data_n, f)
        print(f"--- Successfully saved {processed_count_n} negative spectra data to {output_pickle_file_n} ---")
    except Exception as e:
        print(f"Error saving negative pickle file {output_pickle_file_n}: {e}")

def process_glycyrrhiza_uralensis_data():
    """
    Reads scan information from '甘草-35个.xlsx' sheet 'Sheet1',
    reads corresponding spectrum data from specified .raw files in the same directory,
    assigns classification labels based on '大类' and '小类' columns,
    and saves positive and negative ion mode data to separate pickle files:
    'Glycyrrhiza_uralensis_P.pickle' and 'Glycyrrhiza_uralensis_N.pickle'.
    """
    excel_path = r'C:\work\msclassifier\data\白芷\20251105-白芷化合物分类表.xlsx'
    # The .raw files are in the same directory as the Excel file
    raw_base_path = os.path.dirname(excel_path)

    # New output pickle file names for positive and negative data
    output_pickle_file_p = os.path.join(raw_base_path, 'bz_positive.pickle')
    output_pickle_file_n = os.path.join(raw_base_path, 'bz_negative.pickle')

    sheet_name = 'Sheet1'

    print(f"--- Starting processing from Excel: {excel_path}, Sheet: {sheet_name} ---")

    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except FileNotFoundError:
        print(f"Error: Excel file not found at {excel_path}")
        return
    except KeyError:
        print(f"Error: Sheet '{sheet_name}' not found in {excel_path}")
        return
    except Exception as e:
        print(f"Error reading Excel file {excel_path}: {e}")
        return

    # Initialize separate lists for positive and negative spectra
    all_spectra_data_p = []
    all_spectra_data_n = []

    processed_count_p = 0
    processed_count_n = 0

    # Dictionary to store opened MSFileReader objects to avoid reopening files
    raw_file_readers = {}

    # Helper function to extract scan number from string like "P-2614" or "N-2574"
    def extract_scan_number(scan_str_with_prefix):
        try:
            # Remove prefix (e.g., "P-" or "N-") and convert to integer
            return int(str(scan_str_with_prefix).split('-')[1])
        except (ValueError, AttributeError, IndexError):
            print(f"Warning: Could not extract scan number from '{scan_str_with_prefix}'. Skipping.")
            return None

    # Group the DataFrame by '原始数据名称' to process each raw file's scans efficiently
    for raw_filename_in_excel, group_df in df.groupby('原始数据名称'):
        # Construct the full path to the .raw file by appending '.raw'
        full_raw_file_path = os.path.join(raw_base_path, raw_filename_in_excel) + '.raw'

        if not os.path.exists(full_raw_file_path):
            print(f"Warning: Raw file not found at {full_raw_file_path}. Skipping data from this file.")
            continue

        # Open the MSFileReader for the current raw file if not already open
        if full_raw_file_path not in raw_file_readers:
            try:
                raw_file_readers[full_raw_file_path] = MSFileReader(full_raw_file_path)
                print(f"Processing spectra from raw file: {raw_filename_in_excel}...")
            except Exception as e:
                print(f"Error opening raw file {full_raw_file_path}: {e}. Skipping spectra from this file.")
                continue

        current_raw_reader = raw_file_readers[full_raw_file_path]

        for index, row in group_df.iterrows():
            scan_num_str = row.get('二级谱图编号') # Column containing scan numbers like "P-2614"
            scan_num = extract_scan_number(scan_num_str)

            if scan_num is None:
                continue

            # Extract classification labels directly from the DataFrame row
            fixed_class = row.get('大类').lower()
            fixed_sub_class = row.get('小类').lower()
            print(fixed_class, fixed_sub_class)
            # Ensure classification labels are not missing
            if pd.isna(fixed_class) or pd.isna(fixed_sub_class):
                print(f"Warning: Missing '大类' or '小类' for scan {scan_num_str} in {raw_filename_in_excel}. Skipping.")
                continue

            try:
                label_data = current_raw_reader.GetLabelData(scan_num)
                collision_energy = current_raw_reader.GetCollisionEnergyForScanNum(scan_num, 2)
                scan_start_time = current_raw_reader.GetScanHeaderInfoForScanNum(scan_num)['StartTime']
                tic = current_raw_reader.GetScanHeaderInfoForScanNum(scan_num)['TIC']
                if label_data and len(label_data) > 0:
                    mz_values = label_data[0][0]
                    intensity_values = label_data[0][1]

                    # Determine ion mode from the scan number string prefix
                    ion_mode = 'positive' if scan_num_str.startswith('P-') else 'negative'
                    name = row['英文名']
                    current_spectrum_data = {
                        'peaks': {
                            'mz': np.array(mz_values, dtype=np.float32),
                            'intensities': np.array(intensity_values, dtype=np.float32)
                        },
                        'metaData': [
                            {'name': 'scanNum', 'value': scan_num},
                            {'name': 'raw filename', 'value': raw_filename_in_excel},
                            {'name': 'ion mode', 'value': ion_mode},
                            {'name': 'tic', 'value': tic},
                            {'name': 'name', 'value': name},
                            {'name': 'collision_energy', 'value': collision_energy},
                            {'name': 'scan_start_time', 'value': scan_start_time},
                            {'name': 'precursor m/z', 'value': 1},
                        ],
                        'classification_L:': [
                            {'name': 'class', 'value': fixed_class},
                            {'name': 'sub_class', 'value': fixed_sub_class}
                        ]
                    }

                    # Append to the correct list based on ion mode
                    if ion_mode == 'positive':
                        all_spectra_data_p.append(current_spectrum_data)
                        processed_count_p += 1
                    else:
                        all_spectra_data_n.append(current_spectrum_data)
                        processed_count_n += 1
                else:
                    print(f"No label data found for scan number {scan_num} in {raw_filename_in_excel}. Skipping.")
            except Exception as e:
                print(f"Error processing scan {scan_num} from {raw_filename_in_excel}: {e}. Skipping.")

    # Close all opened MSFileReader objects to release resources
    for reader in raw_file_readers.values():
        if reader:
            reader.Close()

    # Save positive ion data to its own pickle file
    try:
        output_dir_p = os.path.dirname(output_pickle_file_p)
        if output_dir_p and not os.path.exists(output_dir_p):
            os.makedirs(output_dir_p, exist_ok=True) # Create directory if it doesn't exist
        with open(output_pickle_file_p, 'wb') as f:
            pickle.dump(all_spectra_data_p, f)
        print(f"\n--- Successfully saved {processed_count_p} positive spectra data to {output_pickle_file_p} ---")
    except Exception as e:
        print(f"Error saving positive pickle file {output_pickle_file_p}: {e}")

    # Save negative ion data to its own pickle file
    try:
        output_dir_n = os.path.dirname(output_pickle_file_n)
        if output_dir_n and not os.path.exists(output_dir_n):
            os.makedirs(output_dir_n, exist_ok=True) # Create directory if it doesn't exist
        with open(output_pickle_file_n, 'wb') as f:
            pickle.dump(all_spectra_data_n, f)
        print(f"--- Successfully saved {processed_count_n} negative spectra data to {output_pickle_file_n} ---")
    except Exception as e:
        print(f"Error saving negative pickle file {output_pickle_file_n}: {e}")


if __name__ == "__main__":
    # load_corydalis_decumbens()
    # process_iris_domestica_1_data()
    # process_iris_domestica_2_data()
    process_glycyrrhiza_uralensis_data()
