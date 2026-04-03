# About MS-SMART 2.0

## What is MS-SMART 2.0?

MS-SMART 2.0 can effectively simplify the identification process of compounds in complex samples and enhance the efficiency of non-targeted lipidomics data analysis. Compared with version 1.0, MS-SMART 2.0 incorporates neutral loss filtering as a supplementary or alternative approach. Firstly, diagnostic ions or characteristic neutral loss information is obtained based on existing references. Then, the collected LC-MS/MS data is automatically filtered using these characteristic information. Finally, the filtering results are merged by calculating the mass (m/z) differences between different adduct ions of phospholipid compounds in positive and negative ion modes.

## Repository Purpose

This repository houses a suite of analytical tools, including programs to:

- Filter targeted compounds through diagnostic ions or neutral loss
- Recommendation of the classified scheme for filtered substances

## System Requirements

- python-3.9.19
- ipython-8.5.0
- comtypes-1.1.14
- pandas-1.4.4
- openpyxl-3.0.10
- matplotlib-3.5.2
- scipy-1.9.1
- plotly-5.15.0
- scilit-learn-1.2.2

## Features

### find_common_ions

find_common_ions program extracts common ions from reference datasets.

#### How to Use find_common_ions

Execute the following command:

```bash
python find_common_ions.py
```

The program will process reference data from `original_data/example_data.csv`. Output will be stored in the `results/` directory.

Parameters can be modified in the `find_common_ions.conf` file, including:

- `ppm` (mass tolerance for peak merging across spectra): default is 5.
- `min_spectrums_num` (minimum number of spectra required): default is 3.
- `top_n` (number of most common ions selected for diagnostics): default is 5.

### filter_compounds_and_find_substituents

filter_compounds_and_find_substituents efficiently filter targeted compounds from LC-HRMS datasets.

#### How to Use filter_compounds_and_find_substituents

Execute the following command:

```bash
python filter_compounds_and_find_substituents.py
```

This tool filters compounds from LC-HRMS data (e.g., `original_data/20251025-SP-YY-N.raw`) based on diagnostic ions recorded found in the literature. All results are preserved in the `results/` directory.

Modify the `filter_compounds_and_find_substituents.conf` file to adjust settings as necessary:

- `process_file`: The LC-HRMS data file for extracting targeted components.
- `keyion`: A item for extracting diagnostic ions; set it to `301.2173, 327.2330` if the EPA/DHA ions are extracted.
- `time_delta_threshold`: The difference between the optimal response spectrum time and the output time; set to `0.2` if keep the difference range within 0.2 minutes.

#### How to Use Neutral Loss

Execute the following command:

```bash
python filter_compounds_and_find_substituents.py
```

This tool filters compounds from LC-HRMS data (e.g., `original_data/20251025-SP-YY-P.raw`) based on the neutral loss characteristics recorded found in the literature. All results are preserved in the `results/` directory.

Modify the `filter_compounds_and_find_substituents.conf` file to adjust settings as necessary:

- `process_file`: The LC-HRMS data file for extracting targeted components.
- `keyion`: A item for extracting diagnostic ions; set it to `141.0191` if the characteristic neutral loss of PE is to be extracted.
- `time_delta_threshold`: The difference between the optimal response spectrum time and the output time; set to `0.2` if want to keep the difference range within 0.2 minutes.

### python merge

Match the filtered results and use python merge to recommend possible phospholipid categories.

#### How to Use the python merge

Execute the following command:

```bash
python_merge.py
```

This tool merges compounds from 2 filter table (e.g., `Extracted EPA,DHA.csv`; `Extracted PC.csv`) based on the mass differences recorded found in the literature. All results are preserved in the `merge_results/` directory.

Modify the `python_merge.conf` file to adjust settings as necessary:

- `file_path_a`: The table to be merged: Table 1.
- `file_path_b`: The table to be merged: Table 2. Table 2 requires naming with `"PC"` if the two tables are merged based on the PC characteristics.
- `results_path`: Name of the merged results table (customized).
- `ARRAGNE`: Matching results are sorted in ascending order by the specified column in Table 1. Set it to `RT` if the matching results are sorted by the RT column in Table 1.

For testing purposes, we offer an example `.raw` file, which you can download from: `test data`

To parse this file, we recommend Thermo's Xcalibur software which is available at `thermo.com`. Additionally, this package depends on the MSFileReader Python bindings from `pymsfilereader`.