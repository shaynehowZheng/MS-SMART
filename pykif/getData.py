import datetime
import math
import os
import sys
from math import exp
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from matplotlib import pyplot as plt
plt.rcParams['font.family'] = 'Arial'
from plotly.subplots import make_subplots
from scipy import exp
from scipy.optimize import curve_fit
from PICore.DetectMethod import Trad_Integral as Ti
from PICore.Events import EventTool as Et
from pykif.MSFileReader import MSFileReader
from scipy.stats import normaltest
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from sklearn.cluster import AgglomerativeClustering

os.environ['OMP_NUM_THREADS'] = '1'
sys.path.append('..')
sys.path.extend([os.path.join(root, name) for root, dirs, _ in os.walk("../") for name in dirs])


def keyionfilter(
    raw_file, key_ion_list, scan_num, mass_tolerance, min_abundance=0, any_ion=False,
    neutral_loss_list=None
):
    massarray = np.array(raw_file.GetLabelData(scan_num)[0][0])
    intensity = np.array(raw_file.GetLabelData(scan_num)[0][1])

    if len(massarray) == 0:
        return False, 0, 0, 0

    total_intensity = max(intensity)
    rel_abundance = intensity / total_intensity * 100

    # --- Original Key Ion Filtering ---
    # Filter massarray and rel_abundance for values within a broader range of key_ion for initial processing
    # This helps in reducing the array size before precise matching
    filtered_indices = (massarray > (min(key_ion_list) - 0.5)) & (massarray < (max(key_ion_list) + 0.5))
    filtered_massarray = massarray[filtered_indices]
    filtered_rel_abundance = rel_abundance[filtered_indices]

    num = 0
    delta_list = []
    keyion_observe_list = []
    rel_abundance_list = []
    if len(filtered_massarray) > 0: # Only proceed if there are ions in the filtered range
        for ion in key_ion_list:
            deltaarray = abs((filtered_massarray - ion) / ion * 1e6)
            min_index = np.argmin(deltaarray)
            keyion_observe = filtered_massarray[min_index]
            rel_abun = filtered_rel_abundance[min_index]
            delta = round(deltaarray[min_index], 2)

            if delta < mass_tolerance and rel_abun >= min_abundance:
                delta_list.append(delta)
                keyion_observe_list.append(round(keyion_observe, 5))
                rel_abundance_list.append(round(rel_abun, 5))
                num += 1

        judge = (num > 0) if any_ion else (num == len(key_ion_list))

        if judge:
            print('AA')
            return True, keyion_observe_list, delta_list, rel_abundance_list

    # --- Neutral Loss Filtering (if original conditions not met) ---
    if neutral_loss_list is not None:

        # Get indices of the top 5 most abundant ions
        sorted_indices = np.argsort(intensity)[::-1]  # Sort in descending order
        top_5_indices = sorted_indices[:min(5, len(intensity))] # Take up to 5, or fewer if less than 5 ions

        top_5_masses = massarray[top_5_indices]
        top_5_intensities = intensity[top_5_indices]
        top_5_rel_abundances = rel_abundance[top_5_indices]

        neutral_loss_observed_pairs_mz = []
        neutral_loss_deltas = []
        neutral_loss_observed_pairs_rel_abun = []

        # Check all pairs for neutral loss
        for i in range(len(top_5_masses)):
            for j in range(i + 1, len(top_5_masses)):
                mz1 = top_5_masses[i]
                mz2 = top_5_masses[j]
                rel_abun1 = top_5_rel_abundances[i]
                rel_abun2 = top_5_rel_abundances[j]

                observed_neutral_loss = abs(mz1 - mz2)

                for expected_nl_mass in neutral_loss_list:
                    delta_nl = abs((observed_neutral_loss - expected_nl_mass) / expected_nl_mass * 1e6)

                    if delta_nl < mass_tolerance:
                        neutral_loss_observed_pairs_mz.append((round(mz1, 5), round(mz2, 5)))
                        neutral_loss_deltas.append(round(delta_nl, 2))
                        neutral_loss_observed_pairs_rel_abun.append((round(rel_abun1, 5), round(rel_abun2, 5)))
                        # If we find any neutral loss that matches, we can return True
                        # You might want to adjust this logic if you need to find ALL specified neutral losses
                        # For now, if at least one pair satisfies any neutral loss, it returns True
                        print('BB')
                        return True, neutral_loss_observed_pairs_mz, neutral_loss_deltas, neutral_loss_observed_pairs_rel_abun

    # If neither condition is met
    return False, 0, 0, 0
def getprecursor_prescannum(rawfile, scan_num):
    precursor = rawfile.GetFullMSOrderPrecursorDataFromScanNum(scan_num, rawfile.GetMSOrderForScanNum(scan_num) - 2)[0]
    for i in range(scan_num - 1, scan_num - 4, -1):
        if rawfile.GetMSOrderForScanNum(i) == 1:
            pre_scan_num = i
            break
    return pd.Series([pre_scan_num, precursor])


def initdata(rawfile, keyion, mass_tolerance, min_abundance=0, any_ion=False, neutral_loss_list=None):
    df_ms2 = pd.DataFrame(list(range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1)), columns=['scanNum'])
    df_ms2['msOrder'] = df_ms2.apply(lambda x: rawfile.GetMSOrderForScanNum(x), axis=1)
    print('num of MS1 spectra:', len(df_ms2[df_ms2['msOrder'] == 1]))
    df_ms2 = df_ms2[df_ms2['msOrder'] == 2]
    print('num of MS2 spectra:', len(df_ms2))
    res = df_ms2.apply(lambda x: keyionfilter(rawfile, keyion, x.scanNum, mass_tolerance, min_abundance=min_abundance, any_ion=any_ion, neutral_loss_list=neutral_loss_list), axis=1)
    res = res.apply(lambda x: pd.Series([x[0], x[1], x[2], x[3]]))
    res = res.rename(columns={0: 'keyfilter', 1: 'keyion_observe', 2: 'delta', 3:'rel_abundance'})
    df_ms2[['keyfilter', 'keyion_observe', 'delta', 'rel_abundance']] = res
    df_ms2[['preScanNum', 'precursor']] = df_ms2.apply(lambda x: getprecursor_prescannum(rawfile, x.scanNum), axis=1)
    df_ms2['keyion'] = str(keyion)
    df_ms2['TIC'] = df_ms2['scanNum'].apply(lambda x: rawfile.GetScanHeaderInfoForScanNum(x)['TIC'])
    df_ms2['StartTime'] = df_ms2['scanNum'].apply(lambda x: rawfile.GetScanHeaderInfoForScanNum(x)['StartTime'])
    df_ms2_all = df_ms2
    df_ms2 = df_ms2[df_ms2.keyfilter == True]
    return df_ms2_all, df_ms2


def getstartnum(rawfile, row, thre):
    startnum = int(row.preScanNum)
    for i in range(int(row.preScanNum) - 1, -1, -1):
        if rawfile.GetMSOrderForScanNum(i) == 1:
            x = np.array(rawfile.GetLabelData(i)[0][0])
            y = list(rawfile.GetLabelData(i)[0][1])
            resbool = (abs(x - row.precursor) / row.precursor * 1e6 <= 5)
            if resbool.any():
                if y[np.where(resbool)[0][0]] >= thre:
                    startnum = i
                else:
                    return startnum
            else:
                return startnum
    return startnum


def getendnum(rawfile, row, thre):
    endnum = int(row.preScanNum)
    for i in range(int(row.preScanNum) + 1, rawfile.LastSpectrumNumber + 1):
        if rawfile.GetMSOrderForScanNum(i) == 1:
            x = np.array(rawfile.GetLabelData(i)[0][0])
            y = list(rawfile.GetLabelData(i)[0][1])
            resbool = (abs(x - row.precursor) / row.precursor * 1e6 <= 5)
            if resbool.any():
                if y[np.where(resbool)[0][0]] >= thre:
                    endnum = i
                else:
                    return endnum
            else:
                return endnum
    return endnum


def getrelIntensity(y):
    df_data = pd.DataFrame(y, columns=['y'])
    res = y / df_data.max() * 100

    return df_data['relIntensity'].to_list()


def getChorInfor(rawfile, x):
    scanNum = rawfile.ScanNumFromRT(x)
    msfilter = rawfile.GetFilterForScanNum(scanNum)
    basepeak = dict(rawfile.GetScanHeaderInfoForScanNum(scanNum))['BasePeakMass']
    masslist = rawfile.GetLabelData(scanNum)
    mass_x = list(masslist[0][0])
    mass_y = list(masslist[0][1])
    sort_mass_y = frozenset(mass_y)
    sort_mass_y = sorted(sort_mass_y, reverse=True)
    topmass = []
    for i in range(3):
        topmass.append(mass_x[mass_y.index(sort_mass_y[i])])

    return pd.Series([scanNum, msfilter, basepeak, topmass])


def getAllChorData(rawfile, startTime, endTime, massTolerance=5, units=1, msfilter='', precursormass='50-1000'):
    """
    :Description  :
    :Param  rawfile:MSFileReader
    :Param  msfilter:
    :Returns  x,y,relIntensity,basepeak:
    """
    df_result = pd.DataFrame()
    rawfile.SetMassTolerance(userDefined=True, massTolerance=massTolerance, units=units)
    allChorData = rawfile.GetChroData(startTime, endTime, precursormass, '', msfilter)
    df_result['x'] = allChorData[0][0]
    df_result['y'] = allChorData[0][1]
    df_result['relIntensity'] = df_result.apply(lambda res: res.y / df_result['y'].max() * 100, axis=1)
    df_result[['scanNum', 'msfilter', 'basepeak', 'topmass']] = df_result['x'].apply(lambda x: getChorInfor(rawfile, x))

    return df_result


def add(y):
    sum_y = np.zeros(len(y[1]))
    for index, item in y.items():
        sum_y += np.array(item)
    return sum_y


def getChorDataByMass(rawfile, precursor, startTime, endTime, massTolerance=5, units=1):
    """
    :Description  :
    :Param  rawfile:MSFileReader
    :Param  massTolerance:
    :Param  units:(amu = 2, mmu = 0, or ppm = 1)
    :Param  precursor:
    :Returns  df_chorData:
    """

    df_chorData = pd.DataFrame(columns=['x', 'y', 'relIntensity', 'name'])
    if len(precursor) == 0:
        chorData = rawfile.GetChroData(startTime, endTime, '50-1000', '', 'Full ms')
        x = list(chorData[0][0])
        y = list(np.linspace(0, 0, len(x)))
        df_chorData = df_chorData.append([{'x': x, 'y': y}], ignore_index=True)
        return df_chorData
    # The type of tolerance value (amu = 2, mmu = 0, or ppm = 1)
    rawfile.SetMassTolerance(userDefined=True, massTolerance=massTolerance, units=units)

    for i in range(len(precursor)):
        chorData = rawfile.GetChroData(startTime, endTime, str(precursor[i]), '', 'Full ms')
        x = list(chorData[0][0])
        y = list(chorData[0][1])
        df_chorData = df_chorData.append([{'x': x, 'y': y, 'name': str(precursor[i])}], ignore_index=True)
        # df_chorData['relIntensity'] = df_chorData.apply(lambda res:res.y/df_chorData['y'].max()*100,axis=1)
    # df_chorData = df_chorData.append([{'x':x,'y':add(df_chorData['y']),'name':'sum'}],ignore_index=True)

    return df_chorData


def getSpectrumData(rawfile, time):
    scanNum = rawfile.ScanNumFromRT(time)
    data = rawfile.GetLabelData(scanNum)
    x = data[0][0]
    y = data[0][1]
    charge = data[0][5]
    relIntensity = getrelIntensity(y)

    return x, y, charge, relIntensity


def scatter3d(data1, data2, min_z=10000):
    fig = make_subplots(rows=2, cols=1, specs=[[{'type': 'scene'}], [{'type': 'scene'}]])
    z_max_t = 0
    for n, data in zip([1, 2], [data1, data2]):
        for index, row in data.iterrows():
            x = row['x']
            z = row['y']
            z_max = max(z)
            z_max_t = max([z_max_t, z_max])
            # if max(z) > 3e8:
            #     z = [m / max(z) * 3e8 for m in z]
            x = [x2 for x1, x2 in zip(z, x) if x1 > min_z]
            z = [m for m in z if m > 0]
            fig.add_trace(go.Scatter3d(
                x=x,
                y=[row['precursor']] * len(x),
                z=z,
                mode='lines',
                line=dict(
                    color='blue',
                    width=4
                ),
                name=''
            ), row=n, col=1)

    fig.update_scenes(
        xaxis=dict(title='time(min)', linecolor='black'),
        yaxis=dict(title='precursor', linecolor='black'),
        zaxis=dict(title='Absolute Intensity', linecolor='black'),
        yaxis_range=[250, 450],  # 
        zaxis_range=[0, 3e8],  # 
        camera=dict(
            center=dict(x=0.5, y=0.5, z=0),
            eye=dict(x=1.6, y=2.1, z=0.2)),
        row=1, col=1  # 
    )
    fig.update_scenes(
        xaxis=dict(title='time(min)', linecolor='black'),
        yaxis=dict(title='precursor', linecolor='black'),
        zaxis=dict(title='Absolute Intensity', linecolor='black'),
        yaxis_range=[250, 450],  # 
        zaxis_range=[0, 3e8],
        camera=dict(
            center=dict(x=0.5, y=0.5, z=0),
            eye=dict(x=1.6, y=2.1, z=0.2)),
        row=2, col=1  # 
    )

    fig.update_layout(
        # title='',
        showlegend=False,
        width=1800,
        height=1800,
        margin=dict(l=0, r=0, b=0, t=0),
        template="none"
    )
    return fig


def scatter3d_s(data1, data2=None, file_path=None, ppm4=None, x_min=0, x_max=30, min_z=1e5):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    data1['x'] = data1['x'].apply(lambda x: [float(item) for item in x.strip('(').strip(')').split(',')])
    data1['y'] = data1['y'].apply(lambda x: [float(item) for item in x.strip('(').strip(')').split(',')])
    
    data1['max_x'] = data1['x'].apply(lambda x: max(x))
    data1['min_x'] = data1['x'].apply(lambda x: min(x))
    x_max = data1['max_x'].max()
    x_min = data1['min_x'].min()
    ax.set_xlim(x_min, x_max)

    data1['max_z'] = data1['y'].apply(lambda x: max(x))
    data1['min_z'] = data1['y'].apply(lambda x: min(x))
    z_max = data1['max_z'].max()
    z_min = data1['min_z'].min()

    z_max = 0
    for index, row in data1.iterrows():
        x = row['x']
        z = row['y']

        y = row['precursor']
        x = [x2 for x1, x2 in zip(z, x) if x1 >= min_z]
        z = [m for m in z if m >= min_z]

        # z = [z[i] for i in range(len(z)) if ((x[i] >= x_min) and (x[i] <= x_max))]
        # x = [val for val in x if ((val >= x_min) and (val <= x_max))]
        if len(z) > 0:
            z_max = max(z_max, max(z))

        y = [y] * len(x)
        ax.plot(x, y, z, color='blue', linewidth=0.7)

    if data2 is not None:
        for index, row in data2.iterrows():
            x = row['x']
            z = row['y']
            x = [x2 for x1, x2 in zip(z, x) if x1 > min_z]
            z = [m for m in z if m > min_z]
            y= [row['precursor']] * len(x)
            ax.plot(x, y, z, color='orange', linewidth=0.7)
            
    x_locator = MaxNLocator(5)
    y_locator = MaxNLocator(4)
    ax.w_yaxis.set_major_locator(y_locator)
    ax.w_xaxis.set_major_locator(x_locator)
    ax.set_xlabel('time(min)')
    ax.set_ylabel('m/z', fontstyle='italic')
    ax.set_zlabel('Intensity')
    ax.view_init(elev=15, azim=240) 
    ax.w_xaxis.pane.fill = False
    ax.w_yaxis.pane.fill = False
    ax.w_zaxis.pane.fill = False

    plt.savefig(f'{file_path}/ppm4={ppm4}_3d.svg', format='svg')
    plt.close()


def spScatter(data, chro_data):
    x_min = 6
    x_max = 23
    fig = make_subplots(rows=4, cols=1)
    fig.update_xaxes(range=(x_min, x_max), row=1, col=1)
    fig.update_xaxes(range=(x_min, x_max), row=2, col=1)
    fig.update_xaxes(range=(x_min, x_max), row=3, col=1)
    fig.update_xaxes(range=(x_min, x_max), row=4, col=1)
    y_max_all = 0
    fig.add_trace(
        go.Scatter(
            x=chro_data['x'],
            y=chro_data['y'],
            name='TIF',
            mode='lines',
            # hovertemplate=(
            #                   data['ylabel'] if data['ylabel'] else 'y') + ': %{y:}<br>' + (
            #                   data['xlabel'] if data['xlabel'] else 'x') + ': %{x:.2f}'
        ), row=1, col=1
    )
    for index, row in data.iterrows():
        y_max = max(row.y)
        y_max_all = max([y_max_all, y_max])

    for index, row in data.iterrows():
        y_max = max(row.y)
        if y_max > y_max_all * 0.1:
            r = 2
        elif y_max > y_max_all * 0.01:
            r = 3
        else:
            r = 4

        fig.add_trace(
            go.Scatter(
                x=row.x,
                y=row.y,
                name=str(row.precursor),
                mode='lines',
                # hovertemplate=(
                #                   data['ylabel'] if data['ylabel'] else 'y') + ': %{y:}<br>' + (
                #                   data['xlabel'] if data['xlabel'] else 'x') + ': %{x:.2f}'
            ), row=r, col=1
        )
    fig.update_xaxes(tickfont={'size': 20}, title_text='time(min)', titlefont={"size": 24})
    fig.update_yaxes(tickfont={'size': 20}, title_text='AU', titlefont={"size": 24})

    # y_max_all = math.ceil(y_max_all)  # 
    fig.update_layout(
        height=1000, width=1000,
        template='simple_white',
        font={'color': '#000', 'family': 'Times New Roman'},
        title=dict(text='TIF and filtered chro', x=0.5, font_size=24),
        hoverlabel=dict(
            font_size=14,
            font_family='Times New Roman'
        ),
        showlegend=True,
        margin=dict(t=50, b=0, l=10, r=10),
        legend=dict(font=dict(size=14)),
    )

    # 
    # if 'annotation' in data.keys():
    #     for item in data['annotation'].values:
    #         fig.add_annotation(x=item[0], y=item[1], text=str(item[1]), showarrow=True, arrowhead=2, arrowsize=1,
    #                            arrowwidth=1,
    #                            xref='x', yref='y', font={'size': 12})
    return fig


def sort_row(data, candidate):
    a = data['row'].iat[0]
    if a != candidate:
        data.loc[data['row'] == candidate, 'row'] = -1
        data.loc[data['row'] == a, 'row'] = candidate
        data.loc[data['row'] == -1, 'row'] = a
    candidate_list = [candidate]
    for index in data.index:
        r = data.loc[index, :]
        b = r['row']
        if b in candidate_list:
            pass
        else:
            candidate += 1
            if b == candidate:
                pass
            else:
                data.loc[data['row'] == candidate, 'row'] = -1
                data.loc[data['row'] == b, 'row'] = candidate
                data.loc[data['row'] == -1, 'row'] = b
            candidate_list.append(candidate)
    return data


def divide_rows(data, y_max, n_cluster):
    label = data['row'].max()
    y_max_2 = y_max[np.where(data['row'] == label)]
    kmeans = AgglomerativeClustering(n_cluster)
    kmeans.fit(y_max_2.reshape(-1, 1))
    labels_2 = kmeans.labels_
    labels_2 = [i + label for i in labels_2]

    labels_2 = np.array(labels_2)
    data_2 = data.loc[data['row'] == label, :]
    data_2.loc[:, 'row'] = labels_2
    candidate = label
    data_2 = sort_row(data_2, candidate)
    max_label = max(labels_2)
    return data_2, max_label


def scatter_group(out_path, data, chro_data, add_annotation=0, x_min=0, x_max=0, showlegend=True):
    line_width = 5
    tickfont_size = 40
    axis_width = 3
    legend_size = 40
    title_size = 60
    data['max_y'] = data['y'].apply(lambda x: max(x))
    a = data['max_y'].max()
    b = data['max_y'].min()
    plot_rows = round(math.log10(a / b)) * 2
    plot_rows = min(plot_rows, 10)
    plot_rows = 8
    data = data.sort_values(by='height', ascending=False)

    y_max = []
    for index, row in data.iterrows():
        y_max.append(row.max_y)
    mean = np.mean(y_max, axis=0)
    std = np.std(y_max, axis=0)

    # 
    y_max = (y_max - mean) / std
    y_max = np.array(y_max)
    data['row'] = 2
    max_label = 2
    n_cluster = int((plot_rows - 1) / 3)
    if n_cluster <= 1:
        n_cluster = 2
    while max_label < plot_rows:
        if sum(data['row'] == data['row'].max()) > 10:
            data_2, max_label = divide_rows(data, y_max, n_cluster)
            data.loc[data_2.index, 'row'] = data_2['row']
        else:
            break
    data['row'][data['row'] > plot_rows] = plot_rows
    plot_rows = data['row'].max()
    
    #
    # df_row = pd.read_csv('GHB_group.csv')
    # df_row = df_row[['cname', 'new']]
    # data = pd.merge(data, df_row, on='cname')
    # data['row'] = data['new']
    # plot_rows = int(data['row'].max())

    fig = plt.figure(figsize=(16, 13))
    ax = fig.add_subplot(111, projection='3d')
    from matplotlib import cm
    data = data.sort_values(by='precursor')
    data['p_r'] = data['precursor'].apply(lambda x: str(round(x, 2)))
    p_r_unique = data['p_r'].unique()
    data['p_r_index'] = data['p_r'].apply(lambda x: list(p_r_unique).index(x))
    colors = cm.rainbow_r(pd.Series(p_r_unique).astype('category').cat.codes / len(p_r_unique))
    # 

    for index, row in data.iterrows():
        x = row['x']
        if row['row'] == 2:
            z = [y / 15 for y in row.y]
        elif row['row'] == 3:
            z = [y / 3 for y in row.y]
        else:
            z = row['y']
        y = [row['row']] * len(x)  # 
        ro = row['p_r_index']
        ax.plot(x, y, z, label=row['p_r'], color=colors[ro], linewidth=0.00001)

        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        verts = [list(zip(x, y, z))]
        poly = Poly3DCollection(verts, color=colors[ro], alpha=1)
        ax.add_collection3d(poly)
    # 
    ax.set_xlabel('time(min)')
    ax.set_ylabel('row')
    ax.set_zlabel('intensity')
    ax.set_ylim(8.5, 1.5)

    # 
    import matplotlib.colors as mcolors
    from matplotlib.cm import ScalarMappable
    norm = mcolors.Normalize(vmin=0, vmax=len(p_r_unique)-1)
    sm = ScalarMappable(cmap=cm.rainbow_r, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.1, orientation='vertical', shrink=0.8)
    ticks = [0, len(p_r_unique)-1]
    ticklabels = [p_r_unique[0], p_r_unique[-1]]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticklabels)
    cbar.set_label('m/z')
    # handles, labels = ax.get_legend_handles_labels()
    # unique_handles = [handles[labels.index(precursor)] for precursor in p_r_unique]
    # ax.legend(unique_handles, p_r_unique, loc='upper right', bbox_to_anchor=(1.5, 1), ncol=2, prop={'size': 10})
    ax.set_xlim(6, 23)
    # 
    ax.set_xticks(range(6, 24, 4))
    ax.set_zticks(range(0, int(1.2e8), int(0.3e8)))
    ax.set_yticks(range(1, 9, 2))
    # 
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    # 
    plt.savefig(f'{out_path}\example.svg', format='svg', dpi=300)
    # plt.show()
    plt.close()

    fig = make_subplots(rows=plot_rows, cols=1)
    exponent = math.floor(math.log10(abs(chro_data['y'].max())))
    if x_min == 0:
        x_min_chro = chro_data['x'].min()
        x_min = round(x_min_chro, 0)
    if x_max == 0:
        x_max_chro = chro_data['x'].max()
        x_max = round(x_max_chro, 0)

    fig.add_trace(
        go.Scatter(
            x=chro_data['x'],
            y=chro_data['y'] / (10 ** exponent),
            name='TIF',
            mode='lines',
            line=dict(
                width=line_width  # 
            )
        ), row=1, col=1
    )
    fig.update_yaxes(tickfont={'size': tickfont_size}, linewidth=axis_width, title={
        'text': f'Intensity / 1e{exponent}',
        'font': {'size': title_size}
    }, row=1, nticks=3)

    for row in range(2, plot_rows + 1):
        df_plot = data[data['row'] == row]
        max_y = df_plot['max_y'].max()
        exponent = math.floor(math.log10(abs(max_y)))
        for index, r in df_plot.iterrows():
            fig.add_trace(
                go.Scatter(
                    x=r.x,
                    y=[y / (10 ** exponent) for y in r.y],  # [y / max_y * 100 for y in r.y]
                    name=f'{r.precursor:.5f}',
                    mode='lines',
                    line=dict(
                        width=line_width  # 
                    )
                ), row=row, col=1
            )
            
            if add_annotation == 1:
                max_index = np.argmax(r.y)
                fig.add_annotation(x=r.x[max_index], y=r.y[max_index] / (10 ** exponent), text=str(r['cname']), showarrow=True,
                                arrowhead=line_width,
                                arrowsize=line_width, arrowwidth=line_width, xref='x', yref='y',
                                font={'size': tickfont_size}, row=row, col=1)
                
        fig.update_yaxes(tickfont={'size': tickfont_size}, linewidth=axis_width,     title={
                'text': f'Intensity / 1e{exponent}',
                'font': {'size': title_size}
            },title_standoff=30, row=row, nticks=3)
    data.loc[index, 'row'] = row

    # 
    for row in range(1, plot_rows + 1, 1):
        fig.add_trace(go.Scatter(
            x=[x_min, x_max],
            y=[0, 0],
            mode='lines',
            showlegend=False,
            line=dict(
                width=line_width  # 
            )
        ), row=row, col=1)

    for row in range(1, plot_rows + 1):
        if row == plot_rows:
            fig.update_xaxes(tickfont={'size': tickfont_size},     title={
                'text': f'time(min)',
                'font': {'size': title_size}
            },row=row, col=1)
        else:
            fig.update_xaxes(tickfont={'size': tickfont_size}, title_text="", visible=True, row=row, col=1)
        # if row >= 9:
        #     fig.update_yaxes(tickformat=".2f", row=row, col=1)
        # elif row >= 6:
        #     fig.update_yaxes(tickformat=".1f", row=row, col=1)
        # else:
        #     fig.update_yaxes(tickformat=".0f", row=row, col=1)
    fig.update_xaxes(range=(x_min, x_max), linewidth=axis_width)

    height = 141.4 * plot_rows * 3
    fig.update_layout(
        height=height, width=1000 * 3,
        template='simple_white',
        font={'color': '#000', 'family': 'Times New Roman'},
        title=dict(text='TIF and filtered chro', x=0.5, font_size=plot_rows * 8),
        hoverlabel=dict(font_size=tickfont_size, font_family='Times New Roman'),
        margin=dict(t=plot_rows * 15, b=plot_rows * 15, l=plot_rows * 5, r=plot_rows * 3),
        legend=dict(font=dict(size=legend_size)),
    )
    if showlegend is False:
        fig.update_layout(
            showlegend=False,
        )
    return fig


def gaussian(x, *param):
    return param[0] * np.exp(-np.power(x - param[1], 2.) / (2 * np.power(param[2], 2.)))  # 高斯公式


def gaussian_fit(fit_y):
    # 
    fit_x = np.arange(len(fit_y))
    zMatrix = np.matrix(np.log(fit_y))
    xMatrixT = np.matrix(np.reshape(np.concatenate((np.ones((len(fit_y))), fit_x, fit_x * fit_x)), (3, len(fit_y))))
    xMatrix = np.matrix(xMatrixT.T)
    bMatrix = ((xMatrixT * xMatrix).I * xMatrixT) * zMatrix.T  # 
    b2, b1, b0 = float(bMatrix[2][0]), float(bMatrix[1][0]), float(bMatrix[0][0])
    s = -1 / b2
    xMaxi = s * b1 / 2
    yMaxi = exp(b0 + xMaxi ** 2 / s)
    popt, pcov = curve_fit(gaussian, fit_x, fit_y, p0=[yMaxi, xMaxi, s])
    print(pcov)
    y = gaussian(fit_x, *popt)
    return y


def kernel_interpolation(y_sample, x1, sig):
    gaussian_kernel = lambda x, c, h: np.exp(-(x - x[c]) ** 2 / (2 * (h ** 2)))
    num = len(y_sample)
    w = np.zeros(num)
    int_matrix = np.asmatrix(np.zeros((num, num)))
    for i in range(num):
        int_matrix[i, :] = gaussian_kernel(x1, i, sig)
    w = int_matrix.I * np.asmatrix(y_sample).T
    return w


def kernel_interpolation_rec(w, x1, x2, sig):
    gkernel = lambda x, xc, h: np.exp(-(x - xc) ** 2 / (2 * (h ** 2)))
    num = len(x2)
    y_rec = np.zeros(num)
    for i in range(num):
        for k in range(len(w)):
            y_rec[i] = y_rec[i] + w[k] * gkernel(x2[i], x1[k], sig)
    return y_rec


def widthAtPcHeight(height, bslStartIdx, bslEndIdx, startIdx, endIdx, topIdx, rawData, heightPct, timeInterval):
    """

        :param rawData:
        :param heightPct:
        :param timeInterval:
        :return: pointL,pointR, width
        """
    #
    xheight = height * heightPct / 100
    bls = np.linspace(rawData[bslStartIdx], rawData[bslEndIdx], bslEndIdx - bslStartIdx + 1)
    if xheight + bls[topIdx - startIdx] <= min(rawData[startIdx], rawData[endIdx]):
        return None
    peakSliceL = rawData[startIdx:topIdx + 1]
    peakSliceR = rawData[topIdx:endIdx + 1]
    pointL = np.argmin(abs(peakSliceL - (xheight))) + startIdx
    pointR = np.argmin(abs(peakSliceR - (xheight))) + topIdx
    width = (pointR * timeInterval - pointL * timeInterval) / 60
    return width


def normfun(x, a, mu, sigma):
    pdf = np.exp(-((x - mu) ** 2) / (2 * sigma ** 2)) / a
    return pdf


def gaussian_func(x, mu, sigma, A):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fitting_the_curve(df_f, rawfile, fitting=1):
    count_less_than_8 = 0
    count_fit_fail = 0
    count_no_base = 0
    cutoff_num = 8
    extend_length = 50

    for index, row in df_f.iterrows():
        startTime = row['startTime']
        endTime = row['endTime']
        precursorsection = row['precursorsection']
        precursorchordata = rawfile.GetChroData(startTime, endTime, precursorsection, '', 'full ms', smoothingType=2,
                                                smoothingValue=15)
        x = precursorchordata[0][0]
        y = precursorchordata[0][1]
        x, y = zip(*[(xi, yi) for xi, yi in zip(x, y) if yi != 0])

        if fitting != 1:
            df_f.at[index, 'x'] = np.array(x)
            df_f.at[index, 'y'] = np.array(y)
        else:
            # 
            max_index = np.argmax(y)

            #
            left_x = x[:max_index + 1]
            left_y = y[:max_index + 1]

            right_x = x[max_index:]
            right_y = y[max_index:]
            if (len(left_x) < cutoff_num) and (len(right_x) >= cutoff_num):
                selected = 'right'
                selected_x = right_x
                selected_y = right_y
            elif (len(left_x) >= cutoff_num) and (len(right_x) < cutoff_num):
                selected = 'left'
                selected_x = left_x
                selected_y = left_y

            # 
            # 
            #
            #
            elif (len(left_x) >= cutoff_num) and (len(right_x) >= cutoff_num):
                _, left_p_value = normaltest(left_y)

                # 
                _, right_p_value = normaltest(right_y)

                # 
                if left_p_value < right_p_value:
                    selected = 'left'
                    selected_x = left_x
                    selected_y = left_y
                else:
                    selected = 'right'
                    selected_x = right_x
                    selected_y = right_y
            else:
                # 
                df_f.at[index, 'x'] = np.array(x)
                df_f.at[index, 'y'] = np.array(y)
                count_less_than_8 += 1
                continue
            # 
            mu_guess = selected_x[np.argmax(selected_y)]  # 
            sigma_guess = np.std(selected_x) * 1.5  #
            a_guess = np.max(selected_y)  #

            try:
                params, _ = curve_fit(gaussian_func, selected_x, selected_y, p0=[mu_guess, sigma_guess, a_guess],
                                      maxfev=2000)
                mu, sigma, a = params
                # 
            except RuntimeError:
                mu = None

            if mu is None:
                if selected == 'left':
                    selected = 'right'
                    selected_x = right_x
                    selected_y = right_y
                else:
                    selected = 'left'
                    selected_x = left_x
                    selected_y = left_y
                if len(selected_x) >= 8:
                    mu_guess = selected_x[np.argmax(selected_y)]  # 
                    sigma_guess = np.std(selected_x) * 1.2  # 
                    a_guess = np.max(selected_y)  #
                    try:
                        params, _ = curve_fit(gaussian_func, selected_x, selected_y,
                                              p0=[mu_guess, sigma_guess, a_guess],
                                              maxfev=2000)
                        mu, sigma, a = params
                    except RuntimeError:
                        count_fit_fail += 1
                        df_f.at[index, 'x'] = np.array(x)
                        df_f.at[index, 'y'] = np.array(y)
                        continue
                else:
                    count_fit_fail += 1
                    df_f.at[index, 'x'] = np.array(x)
                    df_f.at[index, 'y'] = np.array(y)
                    continue
            interpolated_x = []
            if selected == 'left':
                start_value = selected_x[0] - extend_length * 0.01
                for i in range(extend_length):
                    value = start_value + i * 0.01
                    interpolated_x.append(value)
                interpolated_x = np.concatenate((interpolated_x, selected_x))
            else:
                for i in range(extend_length):
                    value = selected_x[-1] + i * 0.01
                    interpolated_x.append(value)
                interpolated_x = np.concatenate((selected_x, interpolated_x))

            interpolated_y = gaussian_func(interpolated_x, mu, sigma, a)
            crossing_index = np.where(interpolated_y <= 100)
            if selected == 'left':
                if len(crossing_index[0]) > 0:
                    crossing_index = crossing_index[0][-1]
                else:
                    count_no_base += 1
                    df_f.at[index, 'x'] = np.array(x)
                    df_f.at[index, 'y'] = np.array(y)
                    continue
                if crossing_index < 50:
                    outerpolated_x = interpolated_x[crossing_index:50]
                    outerpolated_y = interpolated_y[crossing_index:50]
                else:
                    outerpolated_x = []
                    outerpolated_y = []
                interpolated_x = interpolated_x[crossing_index:]
                interpolated_y = interpolated_y[crossing_index:]
                new_x = np.concatenate((outerpolated_x, selected_x))
                new_y = np.concatenate((outerpolated_y, selected_y))
            else:
                if len(crossing_index[0]) > 0:
                    crossing_index = crossing_index[0][0]
                else:
                    count_no_base += 1
                    df_f.at[index, 'x'] = np.array(x)
                    df_f.at[index, 'y'] = np.array(y)
                    continue
                cutoff = len(interpolated_x) - 50
                if crossing_index < cutoff:
                    outerpolated_x = []
                    outerpolated_y = []
                else:
                    outerpolated_x = interpolated_x[cutoff:(crossing_index + 1)]
                    outerpolated_y = interpolated_y[cutoff:(crossing_index + 1)]
                interpolated_x = interpolated_x[:crossing_index]
                interpolated_y = interpolated_y[:crossing_index]
                new_x = np.concatenate((selected_x, outerpolated_x))
                new_y = np.concatenate((selected_y, outerpolated_y))

            df_f.at[index, 'fitting'] = 1
            symmetric_y = interpolated_y[::-1]  
            if selected == 'left':
                mid = interpolated_x[-1]
                symmetric_x = [mid + (mid - x) for x in interpolated_x]  
                symmetric_x = symmetric_x[::-1]
                result_x = np.concatenate((new_x, symmetric_x[1:]))  
                result_y = np.concatenate((new_y, symmetric_y[1:]))
            else:
                mid = interpolated_x[0]
                symmetric_x = [mid - (x - mid) for x in interpolated_x] 
                symmetric_x = symmetric_x[::-1]
                result_x = np.concatenate((symmetric_x[:-1], new_x))  
                result_y = np.concatenate((symmetric_y[:-1], new_y))

            # plt.plot(symmetric_x, symmetric_y)
            # plt.title('symmetric')
            # plt.show()

            # plt.plot(result_x, result_y, color='blue', label='result')
            # plt.plot(x, y, color='red', label='origin')
            # plt.legend()
            # plt.title('result')
            # plt.show()
            #
            # plt.plot(interpolated_x, interpolated_y, color='green', label='interpolated')
            # plt.plot(selected_x, selected_y, color='red', label='selected')
            # plt.legend()
            # plt.title('interpolated')
            # plt.show()

            df_f.at[index, 'x'] = result_x
            df_f.at[index, 'y'] = result_y

    print('count_fit_fail:', count_fit_fail, ', count_less_than_8:', count_less_than_8, 'count_no_base:',
          count_no_base)
    return df_f


def drop_near_precursor(df, ppm6):
    df = df.sort_values(by='precursor')  
    df['diff'] = df['precursor'].diff() / df['precursor'] * 1e6  
    df_result = df.loc[df['diff'].abs() > ppm6]  
    for index, row in df_result.iterrows():
        df['tolerance'] = abs(df['precursor'] - row['precursor']) / df['precursor'] * 1e6
        max_row = df.loc[df['tolerance'] < ppm6, 'TIC'].idxmax()
        df_result.loc[index, :] = df.loc[max_row, :]
    df_result = df_result.drop(columns='diff')  
    return df_result

def drop_odd_precursor(df_ms2_kif):
    df_ms2_kif.loc[:, 'precursor_int'] = df_ms2_kif['precursor'].apply(lambda x: math.floor(x))
    df_ms2_kif = df_ms2_kif[df_ms2_kif['precursor_int'] % 2 == 0]
    return df_ms2_kif  

def keyion_processor(file_name, keyion, ppm4, ppm5=5, ppm6=5, if_fitting=1, integrate=3, min_abundance=0, any_ion=False, neutral_loss_list=None):
    file = r'{}.raw'.format(file_name)
    rawfile = MSFileReader(file)

    chro_data = rawfile.GetChroData(rawfile.StartTime, rawfile.EndTime, '1-9999', '', 'full ms',
                                    smoothingType=2, smoothingValue=15)
    chro_df = pd.DataFrame({'x': chro_data[0][0], 'y': chro_data[0][1]})

    df_ms2, df_ms2_kif = initdata(rawfile, keyion, ppm4, min_abundance=min_abundance, any_ion=any_ion, neutral_loss_list=neutral_loss_list)
    df_ms2_kif = drop_odd_precursor(df_ms2_kif)
    print(f'{datetime.datetime.now()}: finished KIF with {len(df_ms2_kif)} spectra')
    
    df_f_all = drop_near_precursor(df_ms2, ppm6)
    df_f_all['x'] = ''
    df_f_all['y'] = ''
    for index, row in df_f_all.iterrows():
        minprecursor = round(row['precursor'] / (ppm5 / 1e6 + 1), 5)
        maxprecursor = round(row['precursor'] / (1 - ppm5 / 1e6), 5)
        precursorsection = str(minprecursor) + '-' + str(maxprecursor)
        chordata = rawfile.GetChroData(rawfile.StartTime, rawfile.EndTime, precursorsection, '', 'full ms',
                                       smoothingType=2, smoothingValue=15)
        df_f_all.loc[index, 'x'] = str(chordata[0][0])
        df_f_all.loc[index, 'y'] = str(chordata[0][1])

    df_fig = pd.DataFrame()
    for index, row in df_ms2_kif.iterrows():
    
        minprecursor = round(row['precursor'] / (ppm5 / 1e6 + 1), 5)
        maxprecursor = round(row['precursor'] / (1 - ppm5 / 1e6), 5)
        precursorsection = str(minprecursor) + '-' + str(maxprecursor)
        precursorRT = rawfile.RTFromScanNum(int(row['preScanNum']))
        chordata = rawfile.GetChroData(rawfile.StartTime, rawfile.EndTime, precursorsection, '', 'full ms',
                                       smoothingType=2, smoothingValue=15)
        timedata = list(chordata[0][0])
        rawdata = list(chordata[0][1])
        timeInterval = round((rawfile.EndTime - rawfile.StartTime) * 60 / rawfile.LastSpectrumNumber, 1)
        TimeData = np.asarray(timedata) * 60
        # if round(row['precursor'], 5) == 352.15439:
        #     plt.plot(timedata, rawdata)
        #     plt.show()
        err, curr = Ti.getTradPreResult(rawdata, timedata, timeInterval, 0, -1, 5, integrate)
        res = Et.calculateResult(curr, rawdata, TimeData, int(round(1 / timeInterval)))
        flag = 0
        for item in res:
            for peak in item['Peaks']:
                peakstart = timedata[peak['StartIdx']]
                peaktop = timedata[peak['ApexIdx']]
                peakend = timedata[peak['EndIdx']]
                # if round(row['precursor'], 5) == 352.11734:
                #    x = timedata[peak['StartIdx']:peak['EndIdx']]
                #    y = rawdata[peak['StartIdx']:peak['EndIdx']]
                #    plt.plot(x, y)
                #    plt.show()
                if precursorRT >= peakstart and precursorRT <= peakend:
                    startTime = peakstart
                    rt = peaktop
                    endTime = peakend
                    area = peak['Area']
    
                    height = peak['Height']
                    raw_w = rawdata[peak['StartIdx']:peak['EndIdx']]
                    raw_w = list(filter(lambda x: x > 0, raw_w))
                    width = len(raw_w)
                    flag = 1
                    y0 = rawdata[peak['ApexIdx']]
                    y1 = rawdata[peak['StartIdx']]
                    y2 = rawdata[peak['EndIdx']]
                    dheight_left = y0 - y1
                    dheight_right = y0 - y2
                    break
    
            if flag:
                break
    
        # plt.scatter(timedata, rawdata)
        # plt.show()
        # time_s = timedata[peak['StartIdx']:peak['EndIdx']]
        # raw_s = rawdata[peak['StartIdx']:peak['EndIdx']]
        # plt.scatter(time_s, raw_s)
        # plt.show()
        dheight_left_percent = 0 if height == 0 else round(dheight_left / height * 100, 1)
        dheight_right_percent = 0 if height == 0 else round(dheight_right / height * 100, 1)
        df_temp = pd.DataFrame([{
            'scanNum': row['scanNum'], 'TIC': row['TIC'],
            'startTime': startTime, 'RT': rt, 'endTime': endTime, 'timedata': timedata,
            'precursor': round(row['precursor'], 5), 'precursorsection': precursorsection,
            'area': area, 'height': height, 'width': width,
            'dheight_left': dheight_left, 'dheight_left(%)': dheight_left_percent,
            'dheight_right': dheight_right, 'dheight_right(%)': dheight_right_percent,
            'delta': row['delta'], 'keyion': row['keyion'], 'keyion_observe': row['keyion_observe'], 'rel_abundance': row['rel_abundance']}])
        df_fig = pd.concat([df_fig, df_temp])
    print(f'{datetime.datetime.now()}: finished ')
    
    df_fig.sort_values(by=['precursor'], ascending=[True], inplace=True)
    df_fig = df_fig.reset_index(drop=True)
    df_f = drop_near_rows(df_fig, ppm6)
    print(f'{datetime.datetime.now()}: finished ')
    
    df_f[['x', 'y', 'fitting']] = None
    df_f = fitting_the_curve(df_f, rawfile, if_fitting)
    df_f = df_f.reset_index(drop=True).reset_index(drop=False)
    df_f = df_f.rename(columns={'index': 'cname'}).drop(['timedata'], axis=1)
    df_f['formula'] = None

    return df_f_all, df_f, chro_df, df_ms2, df_ms2_kif


def drop_near_rows(df, ppm6):
    rows_to_remove = []

    for a in range(len(df)):
        for i in range(a + 1, len(df)):
            if ((math.isclose(df.loc[a, 'startTime'], df.loc[i, 'startTime'], abs_tol=0.1) and
                 math.isclose(df.loc[a, 'endTime'], df.loc[i, 'endTime'], abs_tol=0.1)) or
                math.isclose(df.loc[a, 'RT'], df.loc[i, 'RT'], abs_tol=0)) and \
                    math.isclose(df.loc[a, 'precursor'], df.loc[i, 'precursor'], rel_tol=ppm6 / 1e6):
                if df.loc[a, 'TIC'] < df.loc[i, 'TIC']:
                    rows_to_remove.append(a)
                else:
                    rows_to_remove.append(i)
    df = df.drop(rows_to_remove)

    return df
