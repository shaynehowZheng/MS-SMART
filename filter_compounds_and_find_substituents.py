import math
import shutil
from public_fn import drawAnnotations
import numpy as np
import pandas as pd
from itertools import product
import re
import datetime
import matplotlib.pyplot as plt
from pykif.getData import keyion_processor, scatter3d, scatter3d_s, spScatter, scatter_group
import os
from pykif.MSFileReader import MSFileReader
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
patterns = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']
rx = fr'({"|".join(sorted(patterns, key=len,reverse=True))})(\d+(?:\.\d+)?)?'


warnings.filterwarnings('ignore', message="kurtosistest only valid for n>=20")
pd.options.mode.chained_assignment = None

def plot_from_rawfile(fileplot, precursor_min, precursor_max, start_time, end_time):
    """
    Creates a plot from a raw mass spectrometry file, focusing on a specified m/z (mass-to-charge) range and retention time window.

    Parameters:
    - fileplot: An object representing the mass spectrometry raw file, with methods for accessing scan data.
    - precursor_min: The minimum m/z value for precursor ion filtering.
    - precursor_max: The maximum m/z value for precursor ion filtering.
    - start_time: The start of the retention time window for filtering scans.
    - end_time: The end of the retention time window for filtering scans.

    Returns:
    - A DataFrame containing the scan numbers that fall within the specified m/z range and retention time window. The DataFrame includes columns for scan number, MS order, retention time, and precursor m/z values. If a scan within this filtered DataFrame has the highest Total Ion Current (TIC), only this scan will be returned.
    """
    df_plot = pd.DataFrame(list(range(fileplot.FirstSpectrumNumber, fileplot.LastSpectrumNumber + 1)),
                           columns=['scanNum'])
    df_plot['msOrder'] = df_plot.apply(lambda z: fileplot.GetMSOrderForScanNum(z), axis=1)
    df_plot = df_plot[df_plot['msOrder'] == 2]
    df_plot['RT'] = df_plot['scanNum'].apply(lambda z: fileplot.RTFromScanNum(z))
    df_plot = df_plot.loc[(df_plot['RT'] > start_time) & (
            df_plot['RT'] < end_time), :]
    df_plot['precursor'] = df_plot['scanNum'].apply(lambda z: fileplot.GetFullMSOrderPrecursorDataFromScanNum(
        z, fileplot.GetMSOrderForScanNum(z) - 2)[0])
    df_plot = df_plot.loc[(df_plot['precursor'] > precursor_min) & (
            df_plot['precursor'] < precursor_max), :]
    if len(df_plot) > 0:
        nl_list = []
        for scannum in df_plot['scanNum']:
            TIC = fileplot.GetScanHeaderInfoForScanNum(scannum)['TIC']
            nl_list.append(TIC)
        nl_index = np.argmax(nl_list)
        df_plot = df_plot.iloc[nl_index, :]
    return df_plot


class Substituent:
    def __init__(self, filepath, filename, df, df_sub, df_mo, mass_tolerance, cal_formula, path, plot_ppm, fileplot=False):
        """
         Check if precursor ions appear in the spectral fragments for data. Return information related to the spectrum if the mass-to-charge ratio of the precursor ion and the minimum relative mass error of fragments exceed certain ppm.
        """
        self.filename = filename
        self.df = df
        self.df_columns = self.df.columns
        self.df_sub_all = df_sub
        self.df_mo = df_mo
        self.mass_tolerance = mass_tolerance
        self.mass_tolerance_nl = 30
        self.df_product = pd.DataFrame()
        self.df_sub = self.df_sub_all
        self.column_sub = self.df_sub['substituents'].values
        self.column_ele = ['C', 'H', 'N', 'O']
        self.mass_ele = [12, 1.007825, 14.003074, 15.994914]
        self.cal_formula = cal_formula
        self.path = path
        self.fileplot = fileplot
        self.plot_ppm = plot_ppm
        self.reader = MSFileReader(filepath)

    def cal_sum_formula(self, x):
        """
        Calculates the sum formula for a given row of substituents in a DataFrame.

        Parameters:
        - self: The instance of the class that contains this method. It should have the following attributes:
            - column_ele: A list of string elements representing columns in a DataFrame that are to be used for calculation.
            - column_sub: A list of string substituents representing columns in 'x' DataFrame used for multiplying factors.
            - df_sub: A DataFrame that maps substituents to their respective contributions to the sum formula.
        - x: A row of a DataFrame that contains numeric values for each substituent how many times it is contained in the molecule.

        Returns:
        - A pandas DataFrame with a single row that is the sum formula, which is the numerical summation of all the elements based on the substituents and their counts found in 'x'.

        The method executes the following steps:
        - Initializes an empty DataFrame 'df' with columns corresponding to the elements in 'self.column_ele'.
        - Iterates through the substituent columns 'self.column_sub' in 'x':
            - Multiplies the count of each substituent from 'x' with its respective element contribution from 'self.df_sub'.
            - Adds the result to the 'df' DataFrame, aligning on columns and filling missing values with zero.
        - Sums up the values in 'df' for all elements and transposes the result so that there is one row with a total count for each element.
        - The resulting DataFrame is a one-row DataFrame representing the sum formula for the molecule described by substituents in 'x'.
        """
        df = pd.DataFrame(columns=self.column_ele)
        for column in self.column_sub:
            b = x[column]
            c = self.df_sub.loc[self.df_sub['substituents'] == column, self.column_ele].astype('int')
            d = b * c
            df = df.add(d, fill_value=0).astype('int')
        df = pd.DataFrame(df.sum().astype('int')).T.reset_index(drop=True)
        return df

    def get_formula(self, df_product):
        """
        Constructs the chemical formula as a string from element counts provided in a DataFrame row.

        Parameters:
        - self: The instance of the class that contains this method. It should have the following attribute:
            - column_ele: A list of element symbols (strings), ordered as they should appear in the chemical formula.
        - df_product: A pandas Series or a row from a DataFrame containing the element counts corresponding to the element symbols in 'self.column_ele'.

        Returns:
        - A string representing the chemical formula. Each element symbol is followed by its count if the count is greater than 1, or just the symbol if the count is 1. Elements with a count of 0 are omitted.

        This method proceeds as follows:
        - Loops over each element and its corresponding count taken from 'df_product':
            - If the element's count 'b' is greater than 1, the element's symbol 'a' and its count are appended to the formula string.
            - If the count 'b' is exactly 1, only the element's symbol 'a' is appended.
            - Elements with a zero count are ignored and not included in the formula.
        """
        df_product_ele = df_product[:len(self.column_ele)]
        m = ''
        for a, b in zip(self.column_ele, df_product_ele):
            if b > 1:
                m = f'{m}{a}{int(b)}'
            elif b == 1:
                m = f'{m}{a}'
        return m

    def cal_precursor_formula(self, precursor, h_num):
        """
        Calculates the potential chemical formulas for a precursor ion based on its mass and a specified number of hydrogen atoms to subtract.

        Parameters:
        - self: The instance of the class that contains this method. It should have the following attributes:
            - mass_ele: A list of atomic masses for elements considered in the calculation.
            - column_ele: A list of element symbols corresponding to 'mass_ele'.
            - mass_tolerance: The mass tolerance for matching theoretical and actual mass, in ppm.
        - precursor: The measured mass of the precursor ion.
        - h_num: The number of hydrogen atoms to subtract from the calculated formulas.

        Returns:
        - A pandas DataFrame containing calculated formulas that match the precursor mass within a given tolerance, including adjustments for unsaturation and the subtraction of 'h_num' hydrogen atoms.

        The method undertakes the following steps:
        - Initializes a DataFrame with element masses and symbols.
        - Estimates the number of carbon atoms that could correspond to the precursor mass and creates ranges for possible counts of C, H, N, and O atoms.
        - Generates all combinations of these atom counts and calculates the theoretical mass for each combination.
        - Filters combinations where the difference between the theoretical and actual mass is within a specified tolerance.
        - Calculates the degree of unsaturation for each formula and filters out those with unrealistic values.
        - Subtracts the specified number of hydrogen atoms and generates the chemical formula strings.
        - Measures computation time from start to finish and prints it along with the initial number of possible formulas.

        This function is useful for generating potential chemical formulas for a given precursor ion mass, which can be further analyzed or used in mass spectrometry data processing.
        """
        start = datetime.datetime.now()
        df = pd.DataFrame({'mass': self.mass_ele,
                           'element': self.column_ele})
        num = int(precursor / df['mass'].iloc[0] * 0.85)
        c_list = np.arange(int(num * 0.5), num)
        h_list = np.arange(int(num * 0.3), num * 2)
        n_list = np.arange(0, int(num / 4))
        o_list = np.arange(1, int(num / 2))
        df_product = pd.DataFrame(product(c_list, h_list, n_list, o_list), columns=df['element'].values)
        print('Number of initial formulas:', len(df_product))
        df_product['theoretical_mass'] = round(df_product.dot(df['mass'].values), 6)
        df_product['observe_mass'] = precursor
        df_product['mass_difference'] = df_product['observe_mass'] - df_product['theoretical_mass']
        df_product['ppm_difference'] = round(abs(df_product['mass_difference']) / df_product['observe_mass'] * 1e6, 5)
        df_product = df_product[abs(df_product['ppm_difference']) <= self.mass_tolerance * 1.5]

        if len(df_product) > 0:
            df_product['unsaturation'] = df_product.apply(lambda x: round(x[0] + 1 - (x[1] - x[2]) / 2, 1), axis=1)
            df_product = df_product[(df_product['unsaturation'] >= 0) & (df_product['unsaturation'] <= 100)]

            df_product['H'] -= h_num
            df_product['formula'] = df_product.apply(lambda x: self.get_formula(x), axis=1)
            end = datetime.datetime.now()
            print(end, 'Time taken to calculate precursor ion formulas:', end - start)
        else:
            print('No formulas found')
        return df_product

    def process(self, index, mass_filter):
        """
        Internal function that calculates possible substituent combinations based on a given precursor ion index and a mass filter value pertaining to -glc (glucosyl) group.

        Parameters
        ----------
        index : int
            The index of the precursor ion of interest.
        mass_filter : float
            The mass threshold value for summing ion hydrogen counts.

        Examples
        --------
        >>> df_res = self.process(0, 192.06338399999999)
        >>> df_res
        cname  startTime    endTime  ... predicted_molecular_formula_correct   -O total_number_of_substituent_combinations
        97      0  10.973974  11.236024  ...            1  NaN           5
        8       0  10.973974  11.236024  ...            4  2.0           5
        21      0  10.973974  11.236024  ...            4  3.0           5
        57      0  10.973974  11.236024  ...            4  3.0           5
        975     0  10.973974  11.236024  ...            4  3.0           5

        This method performs the following tasks:
        - Retrieves the measured precursor mass from the internal DataFrame using the given index.
        - Initializes a new DataFrame for storing products of calculation.
        - Determines the set of potential precursor formulas based on the measured mass and whether to calculate or use a predefined formula.
        - Iterates through potential parent nucleus categories and calculates possible combinations of substituents that could add up to the measured precursor mass, within certain constraints.
        - Filters out substituent combinations that exceed the mass filter and match within a specified mass tolerance.
        - Calculates and appends additional data, such as predicted formulas, ppm differences, and substituent compositions expressed as text.
        - Tracks and reports calculation time.
        - Returns a DataFrame containing all potential substituent combinations that fit the criteria for the given precursor ion and mass filter.
        
        Note: The function requires access to several attributes and methods within the containing class to carry out its computations, including lists of possible substituents, mass elements, and class-defined methods for formula calculation.
        """
        precursor = self.df.loc[index, 'precursor']
        df_product_s = pd.DataFrame()
        h_num = 0

        if self.cal_formula == 1:
            df_precursor = self.cal_precursor_formula(precursor, h_num)
        else:
            df_precursor = pd.DataFrame(re.findall(rx, self.df.loc[index, 'formula'])).replace(
                '', 1).set_index(0).T.astype('int').reset_index(drop=True)
            df_precursor['formula'] = df_precursor.apply(lambda x: self.get_formula(x), axis=1)

        start = datetime.datetime.now()
        for index_mo in self.df_mo.index:
            if self.df_mo.loc[index_mo, 'sub_category'] != 'dihydroberberines':
                self.df_sub = self.df_sub_all.loc[self.df_sub_all['substituents'] != '-O']
            else:
                self.df_sub = self.df_sub_all
            self.column_sub = self.df_sub['substituents'].values
            sub_mass = self.df_sub['mass'].values
            mo_mass = self.df_mo.loc[index_mo, 'mass']
            mass_delta = precursor - mo_mass

            if mass_delta > mass_filter:
                self.df_sub['num'] = self.df_sub['mass'].apply(lambda x: round((mass_delta - mass_filter) / x), 0)
                self.df_sub.loc[self.df_sub['substituents'] == '-Glc', 'num'] += 1
                num_list = []
                for index_sub, row in self.df_sub.iterrows():
                    if row['substituents'] == '-Glc':
                        max_num = min(14, (row['num'] + 2))
                        num = np.arange(1, max_num)
                    else:
                        max_num = min(14, (row['num'] + 1))
                        num = np.arange(0, max_num)
                    num_list.append(num)
                df_product = pd.DataFrame(product(*num_list), columns=self.df_sub['substituents'].values)

                print('Initial number of substituent combinations:', len(df_product))
            else:
                self.df_sub.loc[:, 'num'] = self.df_sub['mass'].apply(lambda x: round(mass_delta / x), 0)
                num_list = [np.arange(0, min(14, x + 1)) if x > 0 else [0] for x in self.df_sub['num']]
                df_product = pd.DataFrame(product(*num_list), columns=self.df_sub['substituents'].values)
                print('Initial number of substituent combinations:', len(df_product))

            sum_product = df_product.sum(axis=1)
            df_product = df_product[sum_product <= 14]

            df_product['theoretical_sum_mass_substituents'] = round(df_product.dot(sub_mass), 6)
            df_product['delta_mass_predicted_m/z'] = - mass_delta + df_product['theoretical_sum_mass_substituents']
            df_product['ppm_difference'] = round(abs(df_product['delta_mass_predicted_m/z']) / precursor * 1e6, 5)
            df_product = df_product[abs(df_product['ppm_difference']) <= self.mass_tolerance]
            df_product['parent_nucleus_category'] = mo_mass
            df_product['core_formula'] = self.df_mo.loc[index_mo, 'formula']
            df_product['delta_mass_m/z_parent_nucleus'] = mass_delta

            df_core = pd.DataFrame(re.findall(rx, self.df_mo.loc[index_mo, 'formula'])).replace(
                '', 1).set_index(0).T.astype('int').reset_index(drop=True)

            if len(df_product) == 0:
                print('No substituents found', index)
                continue

            for i in df_product.index:
                df_sum = self.cal_sum_formula(df_product.loc[i, :])
                df_plus = df_sum.add(df_core, fill_value=0).astype('int')
                df_plus = df_plus[self.column_ele].reset_index(drop=True)
                df_plus['H'] -= h_num
                df_precursor_0 = df_precursor.loc[(df_precursor[self.column_ele] == df_plus.values).all(1), :]

                predict_formula = df_precursor_0['formula'].iloc[0]
                df_product.loc[i, 'predicted_molecular_formula'] = predict_formula
                df_product.loc[i, 'predicted_molecular_mass'] = df_precursor_0.iloc[0, :4].dot(self.mass_ele)
                df_product.loc[i, 'delta_mass_predicted_formula_m/z'] = df_product.loc[i, 'predicted_molecular_mass'] - precursor
                df_product.loc[i, 'delta_mass_predicted_formula_m/z_ppm'] = round(
                    df_product.loc[i, 'delta_mass_predicted_formula_m/z'] / precursor * 1e6, 5)

                sub_formula = self.get_formula(df_plus.iloc[0, :])
                df_product.loc[i, 'substituents+formula'] = sub_formula
                df_product.loc[i, 'substituents+parent_nucleus_mass'] = df_product.loc[i, 'theoretical_sum_mass_substituents'] + mo_mass

            df_product['substituent_combinations_text'] = df_product.apply(lambda x: ''.join([str(int(x[j])) + str(
                df_product.columns[j] + '+') for j in range(len(sub_mass)) if x[j] > 0])[:-1], axis=1)
            product_columns = np.concatenate((self.df_columns, df_product.columns))
            df_product[self.df_columns] = [self.df.loc[index, self.df_columns]] * len(df_product)
            df_product = df_product[product_columns]
            if self.cal_formula == 0:
                df_product['predicted_molecular_formula_correct'] = df_product.apply(
                    lambda x: x['formula'].strip('+') == x['predicted_molecular_formula'], axis=1)
            else:
                df_product['predicted_molecular_formula_correct'] = np.nan
            df_product['combine_num'] = int(len(df_product))
            df_product_s = pd.concat([df_product_s, df_product])
        end = datetime.datetime.now()
        print(end, 'Substituent calculation time:', end - start, index)

        df_product_s['total_number_of_substituent_combinations'] = len(df_product_s)
        df_product_s['link'] = ''
        df_product_s['top_20_abundance_fragments_HCD110'] = ''
        df_product_s['top_20_abundance_fragments_HCD50'] = ''
        df_product_s['CH3/CH4_neutral_loss_judgement'] = ''
        df_product_s['minimum_ppm_CH3/CH4_neutral_loss'] = ''
        if isinstance(self.fileplot, str):
            try:
                path_html, mass_top_1, mass_top_2 = self.create_mass_spectra_plot(index)
                if (len(df_product_s) > 0):
                    if len(mass_top_2) > 0:
                        numbers = mass_top_2['x'].iloc[:5].tolist()
                    else:
                        numbers = []
                    new_num = float(df_product_s['precursor'].iloc[0])
                    if new_num not in numbers:
                        numbers.append(new_num)

                    differences = []
                    for i in range(len(numbers)):
                        for j in range(i + 1, len(numbers)):
                            difference = abs(numbers[j] - numbers[i])
                            differences.append(difference)

                    mass_ch3 = sum(i * j for i, j in zip(self.mass_ele, [1, 3, 0, 0]))
                    mass_ch4 = sum(i * j for i, j in zip(self.mass_ele, [1, 4, 0, 0]))
                    count = 0
                    diff_min = 100
                    for num in differences:
                        diff1 = abs(num - mass_ch3) / mass_ch3 * 1e6
                        diff2 = abs(num - mass_ch4) / mass_ch4 * 1e6
                        if (diff1 <= self.mass_tolerance_nl) or (diff2 <= self.mass_tolerance_nl):
                            filter_neutral = df_product_s['-OCH3'] == 0
                            df_product_s.loc[filter_neutral, 'CH3/CH4_neutral_loss_judgement'] = 0

                            diff_min = min(diff_min, diff1, diff2)
                            df_product_s.loc[filter_neutral, 'minimum_ppm_CH3/CH4_neutral_loss'] = diff_min
                            break
                        else:
                            diff_min = min(diff_min, diff1, diff2)
                            count += 1
                    if count == len(differences):
                        filter_neutral = df_product_s['-OCH3'] > 0
                        df_product_s.loc[filter_neutral, 'CH3/CH4_neutral_loss_judgement'] = 0
                        df_product_s.loc[filter_neutral, 'minimum_ppm_CH3/CH4_neutral_loss'] = diff_min

                    df_product_s['link'].iloc[0] = 'file:\\\\\\' + path_html
                    mass_string = ''
                    for index_x, row in mass_top_1.iterrows():
                        mass_string += str(round(row['x'], 5)) + ' (' + str(round(row['rel_abundance'], 1)) + ')' + ' / '
                    df_product_s['top_20_abundance_fragments_HCD110'].iloc[0] = mass_string
                    self.df.loc[index, 'top_20_abundance_fragments_HCD110'] = mass_string

                    if len(mass_top_2) > 0:
                        mass_string = ''
                        for index_x, row in mass_top_2.iterrows():
                            mass_string += str(round(row['x'], 5)) + ' (' + str(round(row['rel_abundance'], 1)) + ')' + ' / '
                        df_product_s['top_20_abundance_fragments_HCD50'].iloc[0] = mass_string
                        self.df.loc[index, 'top_20_abundance_fragments_HCD50'] = mass_string
                    else:
                        self.df = self.df.drop(index)
                        return None
            except Exception as error:
                print('Error in the plotting process:', error)
        if len(df_product_s) == 0:
            df_product_s = pd.DataFrame(self.df.loc[index, self.df_columns]).T
            df_product_s['predicted_molecular_formula_correct'] = 0
            df_product_s['total_number_of_substituent_combinations'] = 0
            df_product_s
        return df_product_s

    def run(self):
        """
        Executes the main process flow of the class.

        This method orchestrates various steps to compute possible substituent combinations for precursor ions listed in a DataFrame and subsequently generates two DataFrames: one with the results of those calculations and another with statistics derived from the results.

        Steps:
        1. Initialization of variables and mass calculations for specific substituents.
        2. Iteration over all precursor ions in the internal DataFrame, executing the process function for each.
        3. Accumulation of the calculation results into a final DataFrame.
        4. Computation of various substituent statistics, such as the counts of specific substituent occurrences and combinations.
        5. Removal of temporary columns and rearranging of final result columns for clarity.
        6. Renaming columns from Chinese to English for the resulting DataFrame.
        7. Calculation of total runtime for the entire process.
        8. Return of the results and statistics DataFrames.

        Returns
        -------
        df_res : DataFrame
            Contains detailed information about each precursor ion, the predicted molecular formulas, substituent combinations, mass differences, and other metrics computed during the process.
        
        df_stat : DataFrame
            Contains statistics on the presence or absence of specific substituents across all analyzed precursor ions.
        
        Examples
        --------
        # Assuming this method is part of a class with a properly initialized DataFrame 'df':
        >>> df_results, df_statistics = self.run()
        >>> df_results.head()  # Displays the first few rows of the results DataFrame.
        >>> df_statistics       # Displays the statistics DataFrame.

        Note: This method relies on the 'process' method defined within the same class and various class attributes that are expected to be present and correctly initialized.
        """
        start = datetime.datetime.now()
        self.df_sub_all[self.column_ele] = 0
        mass_glc = self.df_sub_all.loc[self.df_sub_all['substituents'] == '-Glc', 'mass'].iloc[0]
        mass_filter = (self.df_sub_all.loc[self.df_sub_all['substituents'] == '-CH3', 'mass'].iloc[0] + mass_glc)
        for i in self.df_sub_all.index:
            parse = pd.DataFrame(re.findall(rx, self.df_sub_all.loc[i, 'formula'])).replace('', 1).set_index(0).T
            columns_ele_t = parse.columns
            self.df_sub_all.loc[i, columns_ele_t] = parse.loc[1, columns_ele_t]
        self.df_sub_all.loc[self.df_sub_all['substituents'] == '-O', 'H'] = -2

        df_res = pd.DataFrame()
        num_list = []
        for index in self.df.index:
            res = self.process(index, mass_filter)
            if res is not None:
                df_res = pd.concat([df_res, res])
                num_list.append(len(res))
        self.df['substituent_possibility'] = num_list

        df_res = df_res.drop(self.column_sub, axis=1).drop('theoretical_sum_mass_substituents', axis=1)
        column = ['cname', 'RT', 'precursor', 'formula', 'predicted_molecular_formula', 'delta_mass_predicted_formula_m/z', 'delta_mass_predicted_formula_m/z_ppm',
                  'parent_nucleus_category', 'delta_mass_m/z_parent_nucleus', 'substituent_combinations_text', 'substituents+formula', 'substituents+parent_nucleus_mass', 'delta_mass_predicted_m/z',
                  'ppm_difference', 'predicted_molecular_formula_correct', 'predicted_molecular_formula_correct', 'total_number_of_substituent_combinations', 'link',
                  'top_20_abundance_fragments_HCD110', 'top_20_abundance_fragments_HCD50', 'CH3/CH4_neutral_loss_judgement', 'minimum_ppm_CH3/CH4_neutral_loss']
        df_res = df_res[column]
        df_res = df_res.rename({'cname': 'ID', 'precursor': 'mass'})
        end = datetime.datetime.now()
        print('total time is ', end - start)
        return df_res

    def create_mass_spectra_plot(self, index):
        """
        Generates mass spectra plots for a given index from the mass spectrometer data file.

        This method performs the following tasks:
        - Retrieves the precursor, retention time, and other relevant data from the DataFrame for the specified index.
        - Sets up the precursor mass range for the mass spectrum plot based on ppm tolerances.
        - Retrieves and processes mass spectrum data from the reader object.
        - Sorts the top mass fragments by relative abundance and prepares them for plotting.
        - Constructs mass spectrum plots using the Plotly library and adds annotations.
        - Saves the plot to disk as both an HTML and an SVG file. Optionally, it may be opened in a web browser.
        - Returns the path to the HTML plot and lists of the mass fragments for further analyses.

        Parameters
        ----------
        index : int
            The index of the target entry in the DataFrame, which corresponds to a specific scan in the mass spectroscopy data.

        Returns
        -------
        path_html : str
            The file path to the saved HTML plot.
        mass_top_list_1 : DataFrame
            A DataFrame of the top mass fragments by relative abundance for one plot.
        mass_top_list_2 : DataFrame
            A DataFrame of the top mass fragments by relative abundance for another plot if applicable.
        
        Examples
        --------
        # Assuming this method is part of a class with a properly initialized DataFrame 'df':
        >>> path, top_list_1, top_list_2 = self.create_mass_spectra_plot(171)
        
        Note: The actual implementation of methods like `drawAnnotations` and `plot_from_rawfile` are assumed to exist with meaningful implementations.
        """
        print(index)
        precursor = self.df.loc[index, 'precursor']
        rt = self.df.loc[index, 'RT']
        precursor_max = precursor * (1 + 10e-6)
        precursor_min = precursor * (1 - 10e-6)
        start_time = self.df.loc[index, 'startTime']
        end_time = self.df.loc[index, 'endTime']

        fig = make_subplots(rows=2, cols=1, vertical_spacing=0.06, horizontal_spacing=0.06)
        anns = []
        scannum = self.df.loc[index, 'scanNum']
        precurs = round(self.df.loc[index, 'precursor'], 5)

        data = self.reader.GetMassPrecisionEstimate(scannum)
        x = np.array(data[1])
        y = np.array(data[0])
        df_mass = pd.DataFrame({'x': x, 'y': y}).sort_values(by='y', ascending=False)
        max_intensity = df_mass['y'].iloc[0]
        df_mass['rel_abundance'] = df_mass['y'] / max_intensity * 100
        mass_top_list_1 = df_mass.iloc[:20]

        fig.add_trace(go.Bar(
            x=df_mass['x'],
            y=df_mass['y'],
            width=0.3,
            name="HCD110_{}_{}_{}".format(round(rt, 2), scannum, precurs),
            hovertemplate='%{x:.5f}, %{y:.0f}'
        ))
        customdata = pd.DataFrame({'x': x, 'y': y}).to_dict('records')
        ann = drawAnnotations(customdata, 1, 14.5, 10000000, [])
        for item in ann:
            del item['overlap']
        anns += ann

        fileplot = MSFileReader(self.fileplot)
        df_plot = plot_from_rawfile(fileplot, precursor_min, precursor_max, start_time, end_time)
        mass_top_list_2 = []
        if len(df_plot) > 0:
            precurs = round(df_plot['precursor'], 5)
            scannum = int(df_plot['scanNum'])
            data = fileplot.GetMassPrecisionEstimate(scannum)
            x = data[1]
            y = data[0]

            df_mass = pd.DataFrame({'x': x, 'y': y}).sort_values(by='y', ascending=False)
            max_intensity = df_mass['y'].iloc[0]
            df_mass['rel_abundance'] = df_mass['y'] / max_intensity * 100
            mass_top_list_2 = df_mass.iloc[:20]


            fig.add_trace(go.Bar(
                x=df_mass['x'],
                y=df_mass['y'],
                width=0.5,
                name="HCD305070_{}_{}_{}".format(round(rt, 2), scannum, precurs),
                hovertemplate='%{x:.5f}, %{y:.0f}',
                xaxis='x2', yaxis='y2'))
            customdata = pd.DataFrame({'x': x, 'y': y}).to_dict('records')
            ann = drawAnnotations(customdata, 2, 14.5, 10000000, [])
            for item in ann:
                del item['overlap']
            anns += ann
        else:
            print('No corresponding spectrum found in the 305070 data')

        fig.update_layout(
            annotations=anns,
            template='simple_white',
            font={'color': '#000', 'family': 'Times New Roman'},
            hoverlabel=dict(
                font_size=16,
                font_family='Times New Roman'
            ),
            # showlegend=True,
            xaxis=dict(title='mass', fixedrange=True),
            yaxis=dict(title='intensity'),
            legend=dict(yanchor='top', y=0.94),
            margin=dict(t=5, b=0, l=10, r=10),
        )
        path_html = f'{self.filename}\_ID{index}_ppm4={self.plot_ppm}_{precursor}.html'

        path_jpg = f'{self.filename}\_ID{index}_ppm4={self.plot_ppm}_{precursor}.svg'
        # fig.write_image(path_jpg, width=1200, height=600)
        return path_html, mass_top_list_1, mass_top_list_2


def find_peaks_substituent_plot(process_file, ppm4=5, ppm5=5, ppm6=5, integrate=3, if_substituent=1, plot_file=0, top=5, read_ref=0):
    """
    Processes a given LC-MS/MS file to find peaks, calculate key ion statistics, and plot mass spectra.

    This function will:
    - Locate the directory of the process file and prepare the output directory for the results.
    - Extract the collision energy from the process file name.
    - Read the common ions CSV file filtered by the energy value and abundance rank.
    - Calculate theoretical m/z values for key ions.
    - Process raw data using a defined key-ions processor function.
    - Optionally read a reference file to filter out impurities and poor fittings.
    - Generate and save Total Ion Chromatogram (TIC) for MS2 spectra and Key Ions Filtered MS2 spectra.
    - Variable 's' represents an instance of a Substituent object (presumably defined elsewhere), which is used for further processing.
    - Optionally run a substituents analysis method on the 's' object if specified.
    - Save various outputs and statistics to the designated outputs directory.
    - Produce mass spectra plots.
    
    Parameters
    ----------
    process_file : str
        The name of the LC-MS/MS file to be processed.
    ppm4 : int, optional
        Parts-per-million tolerance for peak detection (default is 5).
    ppm5 : int, optional
        Parts-per-million tolerance used in a secondary processing step (default is 5).
    ppm6 : int, optional
        Parts-per-million tolerance used in a third processing step (default is 5).
    integrate : int, optional
        The integration setting used in data processing (default is 3).
    if_substituent : int, optional
        If set to 1, runs substituent analysis (default is 1).
    plot_file : str, optional
        The name of the file used for plotting, if not the same as the process file (default is 0).
    top : int, optional
        The number of top ions to consider in analysis (default is 5).
    read_ref : int, optional
        If set to 1, will read in a reference file for filtering data (default is 0).
    
    Returns
    -------
    count_info : DataFrame
        A DataFrame with the name of the processed file and count statistics of filtered spectra.
    
    Notes
    -----
    This function assumes the existence of many internal functions (like keyion_processor)
    and classes (like Substituent), which should be defined elsewhere in the codebase.

    Examples
    --------
    >>> count_information = find_peaks_substituent_plot('sample_data_file', ppm4=10, integrate=2, if_substituent=0)
    """
    path = os.getcwd()
    path_data = os.path.join(path, 'original_data')
    process_file_path = os.path.join(path_data, process_file)
    print(process_file)
    energy = process_file.split('HCD')[1]
    file_path = os.path.join(path_data, process_file)
    out_path = os.path.join(path, 'results', process_file)
    file_plot = False
    
    if plot_file != 0:
        file_plot = os.path.join(path_data, plot_file)

    if os.path.exists(out_path):
        pass
    else:
        os.makedirs(out_path)
        print(f"Directory '{process_file}' created successfully!")

    df_keyion = pd.read_csv(f'{path}/results/common_ions_ppm3=5.csv')
    df_keyion = df_keyion[df_keyion['energy'] == int(energy)]
    if len(df_keyion) == 0:
        print('No common ions found for this energy level.')
        return 
    df_keyion = df_keyion[df_keyion['rel_abundance_rank'] <= top]
    keyion = df_keyion['m/z'].tolist()
    keyion_theory = cal_theoretical_mz(keyion)

    # KIF
    df_f_all, df_raw, chro_data, df_ms2, df_ms2_kif = keyion_processor(process_file_path, keyion_theory, ppm4=ppm4, ppm5=ppm5, ppm6=ppm6, integrate=integrate)

    # df_ms2.to_csv(f'{out_path}/df_ms2.csv', encoding='utf_8_sig', index=False)
    # df_ms2_kif.to_csv(f'{out_path}/df_ms2_kif_ppm4={ppm4}_top{top}.csv', encoding='utf_8_sig', index=False)
    # df_raw.to_csv(f'{out_path}/df_raw_ppm4={ppm4}_top{top}.csv', encoding='utf_8_sig', index=False)
    # chro_data.to_csv(f'{out_path}/chro_data.csv', encoding='utf_8_sig', index=False)

    # chro_data = pd.read_csv(f'{out_path}/chro_data.csv')
    # df_ms2 = pd.read_csv(f'{out_path}/df_ms2.csv')
    # df_ms2_kif = pd.read_csv(f'{out_path}/df_ms2_kif_ppm4={ppm4}_top{top}.csv')
    # df_raw['y'] = df_raw['y'].apply(lambda x: np.fromstring(' '.join(x.strip('[]').split()), sep=' '))  #
    # df_raw['y'] = df_raw['y'].apply(lambda x: [float(i) for i in x])
    # df_raw['x'] = df_raw['x'].apply(lambda x: np.fromstring(' '.join(x.strip('[]').split()), sep=' '))  #
    # df_raw['x'] = df_raw['x'].apply(lambda x: [float(i) for i in x])

    if read_ref == 1:
        df_ref = pd.read_csv(f'{out_path}/peaks_energy={energy}_ppm4={ppm4}_ppm5={ppm5}_ppm6={ppm6}_top={top}.csv')
        df_ref = df_ref[df_ref['impurity'] == '0']
        df_ref = df_ref[df_ref['fitting'] == 1]
        df_ref = df_ref.set_index('scanNum', drop=False)
        df_raw = df_raw.set_index('scanNum', drop=False)
        df_raw = df_raw.loc[df_ref.index, :]
        df_raw['precursor'] = df_ref['precursor']
        df_raw = df_raw.reset_index(drop=True)
        df_raw.index += 1
        df_raw['cname'] = df_raw.index
    
    count_info = pd.DataFrame({'file_name': [process_file], 'count_ms2': [len(df_ms2)], 'count_ms2_kif': [len(df_ms2_kif)], 'top': [top]})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.vlines(df_ms2['StartTime'], ymin=0, ymax=df_ms2['TIC'], color='black', linewidth=1)
    ax1.set_ylim(bottom=0)
    ax1.set_ylabel('TIC')
    ax2.vlines(df_ms2_kif['StartTime'], ymin=0, ymax=df_ms2_kif['TIC'], color='black', linewidth=1)
    ax2.set_xlabel('RT(min)')
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel('TIC')
    plt.savefig(f'{out_path}/TIC_ppm4={ppm4}_top{top}.svg')
    plt.close()
    scatter3d_s(df_f_all, df_raw, out_path, ppm4)
    
    if process_file == '20221017-DBR-HCD110':
        x_min = 6
    else:
        x_min = 0
    
    fig = scatter_group(out_path, df_raw, chro_data, x_min=x_min, showlegend=True, add_annotation=1)
    # fig.write_image(f'{out_path}\ppm4={ppm4}_kif_{integrate}_fitting.svg')
    df_substituent = pd.read_excel('{}\substituents.xlsx'.format(path_data), sheet_name='substituents', index_col=False)
    df_mother = pd.read_excel('{}\substituents.xlsx'.format(path_data), sheet_name='core', index_col=False)
    s = Substituent(file_path, out_path, df_raw, df_substituent, df_mother, 5, 1, path, ppm4, file_plot)
    
    if if_substituent == 1:
        df_result= s.run()
        len_cname = len(df_result['cname'].unique())
        df_filter = s.df[s.df['substituent_possibility'] > 0]
        print(f'{datetime.datetime.now()}: plot subs, {len(df_filter)}个峰')
        df_result.to_csv(f'{out_path}\substituent_energy={energy}_ppm4={ppm4}_ppm5={ppm5}_ppm6={ppm6}_top={top}_num={len_cname}.csv',
                         encoding='utf_8_sig', index=False)
    
    s.df.to_csv(f'{out_path}\peaks_energy={energy}_ppm4={ppm4}_ppm5={ppm5}_ppm6={ppm6}_top={top}.csv', encoding='utf_8_sig',
                index=False)
    print(f'{datetime.datetime.now()}: finished task %s' % (os.getpid()))
    return count_info

def cal_theoretical_mz(keyion_observed):
    """
    Calculates the theoretical m/z (mass-to-charge ratio) values for a list of observed key ions.

    This function takes a list of observed key ion m/z values and matches them to a predefined list of
    theoretical m/z values. The match is based on the relative tolerance of 50 ppm (parts per million). 
    If an observed m/z value is close enough to a theoretical value within the specified tolerance, 
    the observed value is replaced by the theoretical value.

    Parameters
    ----------
    keyion_observed : list of float
        A list of observed m/z values for key ions.

    Returns
    -------
    keyion_theory : list of float
        A list of m/z values where observed values close to theoretical values have been replaced 
        by those theoretical values.

    Notes
    -----
    The purpose of this function is to correct the observed m/z values to known standard theoretical 
    values. This can help in improving the accuracy of mass spectrometry analysis by ensuring that 
    identified peaks correspond to specific known ions.
    
    The predefined list of theoretical m/z values represents common key ions or fragments observed 
    in mass spectrometry experiments.

    Examples
    --------
    >>> observed_mz = [191.073, 192.081, 190.065]
    >>> theoretical_mz = cal_theoretical_mz(observed_mz)
    >>> print(theoretical_mz)
    [191.07295, 192.08078, 190.06513]
    """
    keyion_theoretical = [191.07295, 192.08078, 190.06513, 204.08078, 180.08078, 208.07569, 177.05730, 178.06513,
                          203.07295, 176.06205]
    keyion_theory = keyion_observed.copy()
    for i in range(len(keyion_theory)):
        for j in range(len(keyion_theoretical)):
            if math.isclose(keyion_theory[i], keyion_theoretical[j], rel_tol=50e-6):
                keyion_theory[i] = keyion_theoretical[j]
                break
    return keyion_theory

if __name__ == "__main__":
    f = open('filter_compounds_and_find_substituents.conf', 'r')  
    config = f.read()
    f.close()
    dbcondict = {}
    for params in config.split():
        item = params.split('=')
        dbcondict[item[0]] = item[1]
    
    find_peaks_substituent_plot(process_file=dbcondict['process_file'], 
                                ppm4=int(dbcondict['ppm4']), 
                                ppm5=int(dbcondict['ppm5']), 
                                ppm6=int(dbcondict['ppm6']), 
                                integrate=int(dbcondict['integrate']), 
                                if_substituent=int(dbcondict['if_substituent']), 
                                plot_file=dbcondict['plot_file'], 
                                top=int(dbcondict['top']))
