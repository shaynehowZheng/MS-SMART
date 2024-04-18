# About MS-SMART
## What is MS-SMART?
MS-SMART is an innovative tool that streamlines the identification and analysis of compounds within complex mixtures across a spectrum of families and genera. Utilizing existing reference data, the system commences by isolating diagnostic ions before skillfully sifting through raw liquid chromatography tandem high-resolution mass spectrometry (LC-HRMS) data to filter out pertinent compounds. In its final act, MS-SMART proposes probable structures for these compounds by employing building blocks identified in various plant-based compounds.

## Repository Purpose
This repository houses a suite of analytical tools, including programs to:
* Detect common ions
* Filter targeted compounds
* Recommend structural compositions of the identified substances.

## System Requirements
* python-3.9.19
* ipython-8.5.0
* comtypes-1.1.14
* pandas-1.4.4
* openpyxl-3.0.10
* matplotlib-3.5.2
* scipy-1.9.1
* plotly-5.15.0
* scilit-learn-1.2.2

## Features
### find_common_ions
The **find_common_ions** program meticulously extracts common ions from a plethora of reference datasets.
### How to Use find_common_ions
Execute the following command: ```python find_common_ions.py```

The program will process reference data from ```original_data/example_data.csv```. Any output will be stored in the ```results/``` directory. For fine-tuned performance, parameters can be modified in the ```find_common_ions.conf``` file, including:
* ppm (mass tolerance for peak merging across spectra): default is 5.
* min_spectrums_num (minimum number of spectra required): default is 3.
* top_n (number of most common ions selected for diagnostics): default is 5.
### filter_compounds_and_find_substituents
**filter_compounds_and_find_substituents** efficiently isolates targeted compounds from LC-HRMS datasets and delivers structural recommendations.
### How to Use filter_compounds_and_find_substituents
Execute the following command: ```python filter_compounds_and_find_substituents.py```.

This tool filters compounds from LC-HRMS data (e.g., ```original_data/20230420-XTW-HCD110.raw```) based on diagnostic ions recorded in ```results/common_ions.csv```. It then suggests structures using building blocks of compounds such as berberines, as documented in ```original_data/substituents.xlsx```. All results are preserved in the ```results/``` directory. Modify the ```filter_compounds_and_find_substituents.conf``` file to adjust settings as necessary:
* process_file: The LC-HRMS data file for extracting targeted components.
* if_substituent: A toggle for structure recommendation of components: default is 1.
* plot_file: The LC-HRMS data file for plotting comparative spectra; set to 0 if plotting comparative spectra is unnecessary.

For testing purposes, we offer an example .raw file, which you can download from:
[test data](https://whueducn-my.sharepoint.com/personal/2011301120022_whu_edu_cn/_layouts/15/onedrive.aspx?id=%2Fpersonal%2F2011301120022%5Fwhu%5Fedu%5Fcn%2FDocuments%2FXTW&ga=1)

To parse this file, we recommend installing Thermo's Xcalibur software available at [thermo.com](www.thermofisher.com). Additionally, this package depends on the MSFileReader Python bindings from [pymsfilereader](https://github.com/frallain/pymsfilereader).
