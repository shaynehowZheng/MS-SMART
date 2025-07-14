import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import plotly.express as px
import os
import re
import json
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False


def find_common_ion(df, totalFileNum):
    """
    Identifies common ions across multiple files in a DataFrame based on relative abundance and cname counts.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing the mass spectrometry data. It should have at least two columns:
        'relAbundance' representing the relative abundance of each ion, and 'cname' which includes
        a list of file identifiers where the ion was found.

    totalFileNum : int
        The number of different files (or total experiments) from which an ion must be present to be 
        considered a common ion. It's a threshold that determines how common an ion must be across 
        different files/experiments to be included in the output list.

    Returns:
    common_ion : list
        A list of mass values corresponding to the ions that are common across the given number of files.
        These ions are sorted by their relative abundance in descending order, with the most abundant ions first.

    The function works by first sorting the DataFrame `df` by relative abundance in descending order,
    ensuring the most abundant ions are considered first. Then, it uses the `apply` function on the 'cname'
    column to filter out only those ions which exist across a number of files equal to or greater than
    `totalFileNum`. The mass values of the common ions are then returned as a list.
    """
    df = df.sort_values(by='relAbundance', ascending=False)
    substract = df['cname'].apply(lambda x: len(set(x)) >= totalFileNum)
    common_ion = df.loc[substract, 'mass'].to_list()
    return common_ion


def check_common_ions_core(df_cname_rid, mass, ppm3, n, df_row, min_abundance=0):
    """
    This function checks for the closest ion mass in a DataFrame to a given target mass within a specified ppm tolerance and updates the count if the closest mass is within that tolerance. It appends the closest mass, its ppm difference, and relative abundance to a given list.

    Parameters:
    df_cname_rid : pandas.DataFrame
        The DataFrame containing mass spectrometry ions data. Expected to have the columns 'mass' and 'relAbundance' for mass-to-charge ratio and relative abundance, respectively.
        
    mass : float
        The target mass/charge (mass) ratio for which the closest ion in the DataFrame is to be found.
        
    ppm3 : float
        The parts per million (ppm) tolerance within which the closest ion's mass must fall to be considered a match.
        
    n : int
        Count of how many times a mass within the specified ppm tolerance has been found. It gets updated if a new closest mass falls within the tolerance.
        
    df_row : list
        The list to which the closest mass, its ppm difference from the target mass, and its relative abundance are appended.

    Returns:
    n : int
        Updated count after checking the given mass against the DataFrame's 'mass' values.
        
    df_row : list
        Updated list with appended closest mass, its mass ppm difference from the target mass, and its relative abundance.

    The function first identifies the ion with the closest mass to the target mass by computing the absolute difference and sorting these differences, selecting the ion with the smallest difference. Then it calculates the ppm difference of the closest found mass from the target mass. If this ppm difference is within the specified tolerance, it increments the count. The closest mass, its ppm difference, and the relative abundance are appended to the `df_row` list.
    """
    closet = (df_cname_rid['mass'] - mass).abs().argsort()[:1]
    mass_cloeset = df_cname_rid['mass'].iloc[closet].values[0]
    relAbundance = df_cname_rid['relAbundance'].iloc[closet].values[0]
    mass_ppm = round((mass_cloeset - mass) / mass * 1e6, 1)
    if (abs(mass_ppm) <= ppm3) & (relAbundance > min_abundance):
        n += 1
    df_row.append(mass_cloeset)
    df_row.append(mass_ppm)
    df_row.append(relAbundance)
    return n, df_row


def check_common_ions(df, mass_check_list, energy, category, ppm2=5, ppm3=5, any_ion=False, min_abundance=0):
    """
    This function processes a DataFrame containing mass spectrometry data to analyze and classify true positives, false negatives, and false positives based on the presence of specified ion masses, their abundance, and given ppm tolerances. It outputs these classifications along with statistics on detection performance.

    Parameters:
    - df : pandas.DataFrame - The DataFrame containing columns for mass spectrometry data, including 'cname' for naming and 'category' for classification.
    - mass_check_list : list - A list of ion masses to check within the provided DataFrame.
    - category, energy, ppm3 : Various parameters used for filtering the DataFrame and defining the ppm tolerance level for ion mass matching.
    - FilePath : str - The path where the resulting Excel file containing the analysis results will be stored.

    Returns:
    - result : list - A summary of the analysis, containing the category, ppm3, mass_check_list, number of false negatives, count of rid(s), ratio of false negatives, number of false positives, uncounted rid(s), and ratio of false positives.
    - false_negative : pandas.DataFrame - A DataFrame listing the false negatives identified during the analysis.
    - false_positive : pandas.DataFrame - A DataFrame listing the false positives identified during the analysis.
    - true_positive : pandas.DataFrame - A DataFrame listing the true positives identified during the analysis.

    The function first checks if the mass_check_list is empty and prints a message indicating whether common ions were found based on the input parameters. It then iterates over unique 'cname' values in the DataFrame to calculate and classify true positives, false negatives, and false positives, creating respective DataFrames for each category. These results, along with the calculated false negative and false positive ratios, are then saved to an Excel file for further analysis.
    """
    columns = ['name', 'category', 'cname']
    for i in range(len(mass_check_list)):
        columns.append(str(mass_check_list[i]))
        columns.append('ppm')
        columns.append('relAbundance')
    false_positive = pd.DataFrame(columns=columns)
    false_negative = false_positive.copy()
    true_positive = false_positive.copy()

    if len(mass_check_list) == 0:
        text = 'category={}, energy={}, ppm3={},not find common ion'.format(category, energy, ppm3)
        print(text)
        result = [category, ppm3, mass_check_list, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        return result
    else:
        text = f'category={category}, energy={energy}, ppm3={ppm3},find common ion{mass_check_list}'
        print(text)
    list_rid = df['cname'].unique()
    list_true_rid = df.loc[df['category'] == category, 'cname'].unique()

    rid_count = 0
    for cname in list_rid:
        df_cname_rid = df.loc[(df['cname'] == cname) & (df['energy'] == energy), :]
        n = 0
        df_row = [df_cname_rid['cname'].values[0], df_cname_rid['category'].values[0], cname]
        for mass in mass_check_list:
            n, df_row = check_common_ions_core(df_cname_rid, mass, ppm3, n, df_row, min_abundance=min_abundance)
        if any_ion:
            judge = n > 0
        else:
            judge = n == len(mass_check_list)
        if (judge):
            if cname in list_true_rid:
                rid_count += 1
                true_positive.loc[str(cname)] = df_row
            else:
                false_positive.loc[str(cname)] = df_row
        else:
            if cname in list_true_rid:
                rid_count += 1
                false_negative.loc[str(cname)] = df_row
            else:
                pass

    rid_total = len(df['cname'].unique())
    len_false_negative = len(false_negative)
    ratio_false_negative = round(len_false_negative / rid_count * 100, 2)

    rid_un_count = rid_total - rid_count
    if rid_un_count == 0:
        rid_un_count = 1
    len_false_positive = len(false_positive)
    ratio_len_false_positive = round(len_false_positive / rid_un_count * 100, 2)
    result = 'check_common_ions_category={}, energy={}, ppm2={}, ppm3={}, false_neg_perc={}%{}={}%, false_pos_perc={}%{}={}%'.format(
        category, energy, ppm2, ppm3, len_false_negative, rid_count, ratio_false_negative, len_false_positive,
        rid_un_count, ratio_len_false_positive)

    false_negative = false_negative.set_index('cname', drop=True)
    false_positive = false_positive.set_index('cname', drop=True)
    true_positive = true_positive.set_index('cname', drop=True)
    with pd.ExcelWriter('{}/results/check_common_ions/{}.xlsx'.format(FilePath, result)) as writer:
        false_negative.to_excel(writer, sheet_name='false_negative')
        false_positive.to_excel(writer, sheet_name='false_positive')
        true_positive.to_excel(writer, sheet_name='true_positive')
    result = [category, ppm3, mass_check_list, len_false_negative, rid_count, ratio_false_negative,
              len_false_positive, rid_un_count, ratio_len_false_positive]
    return result

def ppm_window(mass, ppm=5):
    return mass * (ppm / 1e6)

class ScriptIon:
    def __init__(self, ppm1, ppm2, category, min_number):
        """
        This function processes raw mass spectrometry data, particularly focusing on characteristic ions and analyzes them in a blind sample test. Using the categories present in the dataset, it finds common ions, and checks the instances of false positives and negatives at different energy levels defined by High Collision Dissociation (HCD).

        Parameters:
        - category : str/list - The category or categories to be analyzed. The function handles 'all' categories or a specific list of categories.
        - in_number : int - The minimum number of unique 'cname' values per category required to perform the analysis.
        - ppm1 : float - The primary ppm value used for tolerance in ion mass comparisons and possibly as ppm3.
        - ppm2 : float - The secondary ppm value used for tolerance in subsequent analysis steps.

        The function updates the DataFrame 'df_cate_result' with results including ppm values, energy, categories, and the list of masses checked, alongside the number and percentage of false negatives and positives. It then writes this analysis to an Excel file named 'check_common_ions_result.xlsx' in the specified results directory.
        """
        self.hcd_list = []
        self.ppm1 = ppm1
        self.ppm2 = ppm2
        self.category = category
        self.min_number = min_number
        self.substract_list = []
        self.intensity = []

    def merge_ions(self, df, energy, categories, ppm1=10, ppm2=5, merge_col='cname'):
        """
        This function takes a DataFrame of mass spectrometry data and merges ions based on a permissible error window defined by the parts per million (ppm) parameter values, aiming to create an averaged spectrum. It filters the data by specific categories and energy levels before proceeding with the merge and averaging process.

        Parameters:
        self : object
            The instance of the class that this method is a part of. This parameter is not explicitly used in the function body but is necessary for methods within classes.
        df : pandas.DataFrame
            The DataFrame containing the mass spectrometry data. It is expected to have columns for 'category', the merge column (e.g., 'cname'), 'energy', 'mass', and 'intensity'.
        energy : int
            The specific energy level to filter the DataFrame by before proceeding with the merging process.
        categories : list
            A list of categories to filter the DataFrame by before proceeding with the merging process.
        ppm1 : int, optional
            The first error window for mass/charge (mass) ratio tolerance, specified in parts per million (ppm). This is used initially to find matching ions. The default value is 10 ppm.
        ppm2 : int, optional
            The second, usually smaller, error window for mass ratio tolerance, used for merging ions in the averaged DataFrame. The default value is 5 ppm.
        merge_col : str, optional
            The column name to be used for identifying unique items that should be merged. The default column is 'cname'.

        Returns:
        df_average : pandas.DataFrame
            A new DataFrame containing the merged and averaged ions' data. It includes columns for 'mass', 'intensity', 'max_mass', 'min_mass', list of identifiers from merge_col (e.g., 'cname'), list of 'ionIDs', and the number of items merged ('ridNum'), along with their calculated relative abundance.

        The function performs several key operations:
        1. Filters the input DataFrame by specified categories and energy level.
        2. Iterates over the filtered DataFrame, identifying ions to be merged based on the error windows (ppm1 and ppm2) and calculates averaged values for mass/charge ratio, intensity, etc.
        3. Merges ions that are within the specified error tolerance and recalculates averaged values if matches are found within the newly averaged DataFrame.
        4. Compiles a new DataFrame of averaged ion data, including recalculated relative abundances.
        5. Saves this new DataFrame to an Excel file, naming it according to the specified parameters.
        6. Returns the new DataFrame containing the averaged mass spectrometry data.
        """
        # Assuming df is your DataFrame and categories is a list of categories to filter by
        # Generate a list of unique 'cname' values that match the given categories
        true_rid_list = df.loc[df['category'].isin(categories), 'cname'].unique()

        # Further filter the DataFrame to include only those rows where 'cname' is in true_rid_list and 'energy' matches the given energy
        df = df[(df['cname'].isin(true_rid_list)) & (df['energy'] == energy)]

        # Initialize a new DataFrame to store average values with specified columns
        df_average = pd.DataFrame(columns=['mass', 'intensity', 'max_mass', 'min_mass', 'rid_list', 'ionIDs', 'ridNum'])

        # Iterate through the index of the filtered DataFrame
        for index in df.index:
            # Check if the index exists in the DataFrame index, otherwise skip to next iteration
            if index not in df.index:
                continue

            # Retrieve the mass/charge ratio (mass) for the current index
            mass = df.loc[index, 'mass']
            
            # Calculate the error window for the mass using a custom ppm_window function
            delta = ppm_window(mass, ppm1)

            # Find all the fragments in the spectra that are within the error range of the mass
            between_index = (df['mass'] >= mass - delta) & (df['mass'] <= mass + delta)
            matches = df[between_index]

            # Calculate the average intensity of the matched fragments
            avg_intensity = matches['intensity'].mean()
            # Find the maximum and minimum mass values of the matched fragments
            max_mass = matches['mass'].max()
            min_mass = matches['mass'].min()
            
            # Calculate the new average mass using max_mass and min_mass
            mass = (max_mass + min_mass) / 2
            # Create lists of rid and peakID from the matched fragments
            rid = matches[merge_col].to_list()
            ridNum = len(rid)
            peakID = matches['_id'].to_list()
            
            # Check for existing mass within the new error window in the averaged DataFrame
            delta = ppm_window(mass, ppm2)
            filter = (df_average['mass'] >= mass - delta) & (df_average['mass'] <= mass + delta)
            
            # If there's an existing mass that falls within the error window, merge the data
            if filter.any():
                df_f = df_average[filter]
                for index_f in df_f.index:
                    df_filter = df_f.loc[index_f, :]
                    max_mass = max([df_filter['max_mass'], max_mass])
                    min_mass = min([df_filter['min_mass'], min_mass])
                    mass = (max_mass + min_mass) / 2
                    length_1 = len(df_filter[merge_col])
                    length_2 = len(rid)
                    avg_intensity = (df_filter['intensity'] * length_1 + avg_intensity * length_2) / (length_1 + length_2)
                    rid = df_average.loc[filter, merge_col].iloc[0] + rid
                    ridNum = len(rid)
                    peakID = df_average.loc[filter, 'ionIDs'].iloc[0] + peakID
                    if ridNum > 16:
                        print(f'{mass}')
                df_average = df_average[~filter]

            # Add the new row to df_average, ignoring the existing indices in filter
            new_row = pd.DataFrame({'mass': [mass], 'intensity': [avg_intensity], 
                                    'max_mass': [max_mass], 'min_mass': [min_mass],
                                    f'{merge_col}': [rid], 'ionIDs': [peakID], 'ridNum': [ridNum]})
            df_average = pd.concat([df_average, new_row], ignore_index=True)

            # Remove the matched indices from the original DataFrame
            df = df[~between_index]

        # Calculate the relative abundance as a percentage of max intensity
        df_average['relAbundance'] = df_average['intensity'] / df_average['intensity'].max() * 100

        # Write the processed DataFrame to an Excel file with a formatted name indicating the parameters used
        with pd.ExcelWriter('{}/processed_data/merged_peaks_{}_ppm1={}_ppm2={}.xlsx'.format(FilePath, energy, ppm1, ppm2)) as writer:
            df_average.to_excel(writer, sheet_name='merged_peaks')

        # Return the processed average DataFrame
        return df_average
    
    def main(self, df_raw, common_ion, any_ion, min_abundance=0):
        """
        This function processes raw mass spectrometry data, particularly focusing on characteristic ions and analyzes them in a blind sample test. Using the categories present in the dataset, it finds common ions, and checks the instances of false positives and negatives at different energy levels defined by High Collision Dissociation (HCD).

        Parameters:
        - df_raw : pandas.DataFrame - The raw DataFrame containing mass spectrometry data with 'category', 'MSOrder', 'activeType', 'energy', and 'cname' columns.
        - FilePath : str - The directory path where results will be saved.
        - self.category : str/list - The category or categories to be analyzed. The function handles 'all' categories or a specific list of categories.
        - self.min_number : int - The minimum number of unique 'cname' values per category required to perform the analysis.
        - self.ppm1 : float - The primary ppm value used for tolerance in ion mass comparisons and possibly as ppm3.
        - self.ppm2 : float - The secondary ppm value used for tolerance in subsequent analysis steps.
        - common_ion : list/bool - A flag that indicates whether common ions are to be computed on the fly or already provided.

        The function updates the DataFrame 'df_cate_result' with results including ppm values, energy, categories, and the list of masses checked, alongside the number and percentage of false negatives and positives. It then writes this analysis to an Excel file named 'check_common_ions_result.xlsx' in the specified results directory.
        """
        totalFileNum = len(df_raw['cname'].unique())
        df_raw = df_raw[(df_raw['MSOrder'] == 2) & (df_raw['activeType'] == 'HCD')]
        self.hcd_list =  np.sort(df_raw['energy'].unique())
        df_cate_result = pd.DataFrame(
            columns=['ppm1', 'ppm2', 'energy', 'category', 'ppm3', 'mass_check_list', 'false_neg_num',
                     'type_A_num', 'false_nag_perc', 'false_pos_num', 'type_B_num', 'false_pos_perc'])

        categories_list = self.category

        ppm3 = self.ppm1
        for energy in self.hcd_list:  #
            # Finding common ions
            for category in categories_list:
                df_label = df_raw[df_raw['category'] == category]
                if len(df_label['cname'].unique()) >= self.min_number:
                    if not common_ion:
                        df_average = self.merge_ions(
                            df_label, energy, categories=[category], ppm1=self.ppm1, ppm2=self.ppm2)
                        common_ion_t = find_common_ion(df_average, totalFileNum)
                    else:
                        common_ion_t = common_ion
                    result = check_common_ions(df_raw, common_ion_t, energy, category=category, ppm2=self.ppm2, ppm3=ppm3, any_ion=any_ion, min_abundance=min_abundance)
                    row = [self.ppm1, self.ppm2, energy]
                    row.extend(result)
                    row = pd.DataFrame([row], columns=df_cate_result.columns)
                    df_cate_result = pd.concat([df_cate_result, row])

        with pd.ExcelWriter('{}/results/check_common_ions_result.xlsx'.format(FilePath)) as writer:
            df_cate_result.to_excel(writer, sheet_name='ion_cate_result')

def get_ions_by_abundance(df_energy, n):
    """
    Identifies and ranks the most abundant ions from a DataFrame within a specified mass range.

    This function calculates the sum of relative abundance for each unique mass value in the provided DataFrame. It then filters the results to include only those ions within the mass range of 170 to 210, sorting them in descending order based on their summed relative abundance. The function further narrows down the list to the top 'n' most abundant ions, assigns a rank based on their relative abundance, and returns this information in a new DataFrame.

    Parameters:
    - df_energy (DataFrame): The DataFrame containing ion data with columns for 'mass' and 'relAbundance'.
    - n (int): The number of top ions to return based on their relative abundance sum.

    Returns:
    - DataFrame: A DataFrame containing the top 'n' ions, sorted by their summed relative abundance within the specified mass range, with duplicate mass values removed, and a new column 'rel_abundance_rank' indicating the rank of each ion.
    """
    df_sums = df_energy.groupby('mass')['relAbundance'].sum().reset_index()
    df_sums['m/z_int'] = df_sums['mass'].astype(int)
    df_sums = df_sums.sort_values(by=['relAbundance'], ascending=False)
    df_sums = df_sums.loc[(df_sums['m/z_int'] <= 210) & (df_sums['m/z_int'] >= 170), :]
    df_sums = df_sums.drop_duplicates(subset=['m/z_int']).iloc[:n, :]
    df_sums = df_sums.reset_index(drop=True)
    df_sums['rel_abundance_rank'] = df_sums.index + 1
    return df_sums


def find_common_ions(ppm1=5, ppm2=5, min_spectrums_num=3, common_ion=[], any_ion=False, min_abundance=0):
    """
    Finds common ions in spectral data using predefined ppm thresholds and a minimum number of spectrums.

    The function ingests spectral data from a CSV file, identifies the unique categories, and standardizes the filename format. It then triggers a common ion search by utilizing a ScriptIon instance with parameters set as per user preferences. The common ion search is performed within the ScriptIon's main method using the raw DataFrame and the list of common ions.

    Parameters:
    - ppm1 (int): The first ppm (parts per million) threshold for the ion search (default is 5).
    - ppm2 (int): The second ppm threshold for the ion search (default is 5).
    - min_spectrums_num (int): The minimum number of spectrums required for an ion to be considered common (default is 3).
    - common_ion (list): An optional list of pre-identified common ions to consider in the search (default is an empty list).If it is provided, it is used as the common ion list for checking false postive rate and false negative rate.

    This function doesn't have a return statement; results are presumably stored or manipulated within the 'ScriptIon' instance.
    """
    df_raw = pd.read_csv('{}\original_data\example_data.csv'.format(FilePath), index_col=False)
    cate_name = df_raw['category'].unique()
    df_raw['filename_s'] = df_raw['cname'].apply(lambda x: x[:15])
    # 191.07295, 192.08078, 204.08078, 220.07569
    s = ScriptIon(ppm1, ppm2, cate_name, min_spectrums_num)
    s.main(df_raw, common_ion, any_ion, min_abundance=min_abundance)


def transform_ion_data(file, df_s, ridNum):
    """
    Transforms ion data by extracting energy level and ppm threshold from the file name and filtering the DataFrame.

    This function searches for matches in the provided file name using regular expressions to find the energy level and ppm1 value. It then assigns these values to new columns 'energy' and 'ppm' within the DataFrame. The DataFrame is further filtered to include rows where 'ridNum' values are greater than or equal to the specified threshold.

    Parameters:
    - file (str): The file name which contains encoded information about energy levels and ppm values in the format "merged_peaks_{energy}" and "ppm1={value}".
    - df_s (DataFrame): The pandas DataFrame to be transformed.
    - ridNum (int): The threshold number of 'ridNum' below which rows will be filtered out.

    Returns:
    - DataFrame: The transformed DataFrame with additional columns for 'energy' and 'ppm', and filtered according to 'ridNum'.
    """
    energy_match = re.search(r"merged_peaks_(\d+)", file)
    ppm1_match = re.search(r"ppm1=(\d+)", file)
    df_s['energy'] = int(energy_match.group(1))
    df_s['ppm'] = int(ppm1_match.group(1))
    df_s = df_s[df_s['ridNum'] >= ridNum]
    return df_s
    
def script_plot_hot(n=5, ridNum=16):
    """
    Generates hot-spot plots for common ions across different collision energies and ppm levels.

    The function reads all Excel files with spectral data, consolidates the information, and generates scatterplots. It sorts and aggregates data according to the mass values and collision energies, sums the relative abundance of ions, matches mass values to known targets, and saves the results as CSV files. Scatter plots are created to visualize the summed abundance across different mass and collision energies, as well as the top n most abundant ions.

    Parameters:
    - n (int): Number of top ions to be plotted based on abundance.
    - ridNum (int): Identifier number that is used during the transformation of ion data (default is 16).

    Process:
    1. Reads all Excel files in the specified 'FilePath/processed_data' directory.
    2. Transforms and consolidates the ion data into a single DataFrame.
    3. Loops over unique ppm values and generates plots by grouping data by 'mass' and 'energy'.
    4. Assigns known target mass values to integer mass keys and exports this data.
    5. Formats abundance data for various energy levels and exports it.
    6. Generates and exports a scatter plot for the top n abundant ions specified by the user.
    """
    df_con = pd.DataFrame()
    path = os.path.join(FilePath, 'processed_data')
    for file in os.listdir(path):
        if file.endswith('.xlsx'):
            df = pd.read_excel(f'{path}\\{file}', sheet_name='merged_peaks', index_col=0)
            df = transform_ion_data(file, df, ridNum)
            df_con = pd.concat([df_con, df])

    for PPM in df_con['ppm'].unique():
        df_ppm = df_con[df_con['ppm'] == PPM]
        df_group = df_ppm.groupby(['mass', 'energy'])['relAbundance'].sum().reset_index()
        fig = px.scatter(df_group, x="mass", y="energy", color='relAbundance')
        fig.update_layout(
            font=dict(family="Times New Roman", size=11), 
            xaxis=dict(title="<i>mass</i>", titlefont=dict(size=20),
                       tickfont=dict(family="Times New Roman", size=11)),
            yaxis=dict(title="collision energy(%)", titlefont=dict(size=20, family="Arial"),
                       tickfont=dict(family="Times New Roman", size=11)),
        )
        path = f'{FilePath}/results/plot/ppm={PPM}_sum_abundance_all.svg'
        # fig.write_image(path)

        df_sums = pd.DataFrame() 
        for Energy in df_ppm['energy'].unique():
            df_energy = df_ppm[df_ppm['energy'] == Energy]
            df_sum = get_ions_by_abundance(df_energy, n)
            df_sum['energy'] = Energy
            df_sums = pd.concat([df_sums, df_sum])

        df_sums.loc[df_sums['m/z_int'] == 176, 'm/z_t'] = '176.06205'
        df_sums.loc[df_sums['m/z_int'] == 177, 'm/z_t'] = '177.05730'
        df_sums.loc[df_sums['m/z_int'] == 178, 'm/z_t'] = '178.06513'
        df_sums.loc[df_sums['m/z_int'] == 180, 'm/z_t'] = '180.08078'
        df_sums.loc[df_sums['m/z_int'] == 190, 'm/z_t'] = '190.06513'
        df_sums.loc[df_sums['m/z_int'] == 191, 'm/z_t'] = '191.07295'
        df_sums.loc[df_sums['m/z_int'] == 192, 'm/z_t'] = '192.08078'
        df_sums.loc[df_sums['m/z_int'] == 203, 'm/z_t'] = '203.07295'
        df_sums.loc[df_sums['m/z_int'] == 204, 'm/z_t'] = '204.08078'
        df_sums.loc[df_sums['m/z_int'] == 208, 'm/z_t'] = '208.07569'
        df_sums.loc[df_sums['m/z_int'] == 174, 'm/z_t'] = '174.04640'
        df_sums.loc[df_sums['m/z_int'] == 188, 'm/z_t'] = '188.04948'
        df_sums.loc[df_sums['m/z_int'] == 189, 'm/z_t'] = '189.05730'
        df_sums.loc[df_sums['m/z_int'] == 194, 'm/z_t'] = '194.06004'
        df_sums.loc[df_sums['m/z_int'] == 201, 'm/z_t'] = '201.05730'
        df_sums.loc[df_sums['m/z_int'] == 202, 'm/z_t'] = '202.06513'
        df_sums.loc[df_sums['m/z_int'] == 206, 'm/z_t'] = '206.06004'
        df_sums.loc[df_sums['m/z_int'] == 207, 'm/z_t'] = '207.06787'
        df_sums.to_csv(f'{FilePath}/results/common_ions_ppm3={PPM}.csv')
        
        df_format = pd.DataFrame({'energy':df_sums['energy'].unique(), 'formatted': None})
        for index, row in df_format.iterrows():
            df = df_sums[df_sums['energy'] == row['energy']]
            formatted = df.apply(lambda row: f"{float(row['m/z_t'])} ({round(row['relAbundance'], 2)})", axis=1).tolist()
            formatted = ', '.join(formatted).replace("'", '')
            df_format.loc[index, 'formatted'] =formatted
        df_format.to_csv(f'{FilePath}/results/common_ion_and_abundance_for_various_energy_ppm3={PPM}.csv')
        df_sums = df_sums.rename(columns={"relAbundance": "abundance"})
        df_sums = df_sums.sort_values(by='m/z_int')
        fig = px.scatter(df_sums, x="m/z_t", y="energy", color='abundance')
        fig.update_layout(
            font=dict(family="Times New Roman", size=11),  
            xaxis=dict(title="<i>mass</i>", titlefont=dict(size=20),
                       tickfont=dict(family="Times New Roman", size=11)),
            yaxis=dict(title="collision energy(%)", titlefont=dict(size=20, family="Arial"),
                       tickfont=dict(family="Times New Roman", size=11)),
        )
        path = f'{FilePath}/results/plot/ppm={PPM}_sum_abundance_top{n}.svg'
        # fig.write_image(path)


if __name__ == '__main__':
    # get current work path
    f = open('find_common_ions.conf', 'r')  
    config = f.read()
    f.close()
    dbcondict = {}
    for params in config.split():
        item = params.split('=')
        dbcondict[item[0]] = item[1]
        
    FilePath = os.getcwd()
    find_common_ions(ppm1=int(dbcondict['ppm1']), 
                     ppm2=int(dbcondict['ppm2']), 
                     min_spectrums_num=int(dbcondict['min_spectrums_num']), 
                     common_ion=json.loads(dbcondict['common_ion']),
                     any_ion=bool(dbcondict['any_ion']),
                     min_abundance=int(dbcondict['min_abundance']))
    # script_plot_hot(n=int(dbcondict['top_n']), 
    #                 ridNum=int(dbcondict['max_spectrums_num']))
