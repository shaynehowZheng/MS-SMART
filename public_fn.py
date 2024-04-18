
# from ipython_genutils.py3compat import annotate
# import dash_table
import plotly.graph_objs as go
# import dash_core_components as dcc
# import dash_html_components as html
# import dash_bootstrap_components as dbc
import pandas as pd
from itertools import chain
import copy
import math
import numpy


def spScatter(data, type='line',range=[],width=0.6):
    fig = go.Figure()
    for i, item in enumerate(data['data']):

        xData = copy.deepcopy(list(item['x'])) 
        yData =copy.deepcopy(list(item['y'])) 
        name=item['name']
        text=item['text']
        # print(item['name'])

        if(type == 'bar' and (None not in xData)):
            
            resY = []
            resX = []
            newText=[]
            for index, v in enumerate(yData):
                # print(item)
                resY.append(yData[index:index+1])
                resX.append(xData[index:index+1])
                newText.append(text[index:index+1])

            for item in list(resY):
                item.append(0)
                item.append(None)
            for item in list(resX):
                item.append(item[0])
                item.append(None)
            for textItem in newText:
                textItem.append(textItem[0])
                textItem.append(None)

            xData = list(chain(*resX))
            yData = list(chain(*resY))
            text=list(chain(*newText))
            # print(yData)
        elif (type == 'line' and (None in xData)):
            newX=[]
            newY=[]
            newText=[]
            spX=numpy.array(xData).reshape(int(len(xData)/3),3)
            spY=numpy.array(yData).reshape(int(len(yData)/3),3)
            spText=numpy.array(text).reshape(int(len(text)/3),3)
            # print(spX)
            for item in spX:
                newX.append(item[0])
            for item in spY:
                newY.append(item[0])
            for item in spText:
                newText.append(item[0])
            xData=newX
            yData=newY
            text=newText

        # print(text)
        fig.add_trace(
            go.Scatter(
                x=xData,
                y=yData,
                name=name,
                connectgaps=False,
                mode='lines',
                line={'width': 1.4},
                text=text,
                # text1=data['topmass'],
                # hoverinfo='y',
                # visible=True if index == 0 else 'legendonly',
                # line=dict(color=color[i]),
                # hovertemplate=(data['ylabel'] if data['ylabel'] else 'y')+': %{y:}<br>'+(data['xlabel'] if data['xlabel'] else 'x')+': %{x:.6f}'
                hovertemplate='<br>'+(data['ylabel'] if data['ylabel'] else 'y')+': %{y:}<br>'
                # +(data['xlabel'] if data['xlabel'] else 'x')+': %{x:.6f}<br>'
                +'%{text}',
                hoverlabel={
                    'font':{'size':20}
                }
            )
        )
        
            
            
            
        
    # fig.add_annotation(x=6.192025, y=102855316.21533203, text='1', showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1,
    #                            xref='x', yref='y', font={'size': 12})
    if 'annotation' in data.keys():
        for item in data['annotation'].values:
            fig.add_annotation(x=item[0], y=item[1], text=str(item[1]), showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1,
                               xref='x', yref='y', font={'size': 12})

    fig.update_layout(
        # autosize=False,
        # width=1000,
        # height=550,
        template='simple_white',
        font={'color': '#000', 'family': 'Times New Roman'},
        xaxis=dict(title=data['xlabel']),
        yaxis=dict(title=data['ylabel'], fixedrange=True),
        hoverlabel=dict(
            font_size=14,
            font_family='Times New Roman'
        ),
        showlegend=True,
        # legend=dict(xanchor='center', x=0.9,yanchor='bottom',y=1.1),
        legend=dict(yanchor='top', y=0.95),
        margin=dict(t=15, b=0, l=10, r=10),
        hovermode="x unified"
    )
    # if (type == 'bar'):
    
    annX=data['data'][0]['x']
    annY=data['data'][0]['y']
    if (None in annX and None in annY):
        newX=[]
        newY=[]
        spX=numpy.array(annX).reshape(int(len(annX)/3),3)
        spY=numpy.array(annY).reshape(int(len(annY)/3),3)
        # print(spX)
        for item in spX:
            newX.append(item[0])
        for item in spY:
            newY.append(item[0])
        annX=newX
        annY=newY
    # print(annX)
    customdata = pd.DataFrame(
                {'x': annX, 'y': annY}).to_dict('records')
    ann=drawAnnotations(customdata,width,80000000,range)
            # print(ann)
    for item in ann:
        del item['overlap']
    fig.update_layout(
        annotations=ann
    )
    return fig


def massSpectrum(data,range=[],width=14.5):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=data['x'],
            y=data['y'],
            name=data['time'],
            # width=0.03
            connectgaps=False,
            mode='lines',
            line={'width': 1}
        )
    )
    customdata = pd.DataFrame(
        {'x': data['x'], 'y': data['y']}).to_dict('records')
    
    # print(drawAnnotations(customdata))

    fig.update_layout(
        template='simple_white',
        font={'color': '#000', 'family': 'Times New Roman'},
        hoverlabel=dict(
            font_size=14,
            font_family='Times New Roman'
        ),
        showlegend=True,
        xaxis=dict(title=data['xlabel']),
        yaxis=dict(title=data['ylabel'], fixedrange=True),
        # legend=dict(xanchor='center', x=0.9,yanchor='bottom',y=1.1),
        legend=dict(yanchor='top', y=0.94),
        margin=dict(t=5, b=0, l=10, r=10),
    )
    ann=drawAnnotations(customdata,width,10000000,range)
    for item in ann:
        del item['overlap']
    fig.update_layout(annotations=ann)
    return fig


def drawAnnotations(customdata: list, i, labelWidth,labelHeight, range=None):
    alllabels = []
    labelData = sorted(customdata, key=lambda label: label['y'], reverse=True)
    if(range != [] and range != None):
        labelData = list(filter(lambda item: (item['x'] > math.floor(
            range[0])) and (item['x'] < math.ceil(range[1])), labelData))
    if i == 1:
        for data in labelData:
            annotations = {
                'xref': 'x',
                'yref': 'y',
                'x': data['x'],
                'y': data['y'],
                'xanchor': 'center',
                'yanchor': 'bottom',
                'text': round(data['x'], 5),
                'overlap': False,
                'showarrow': False,
                'opacity': 1,
                'arrowhead': 0,
                'arrowwidth': 0.5,
                'arrowcolor': 'blue',
                'ax': 0,
                'ay': 0,
                'font': {
                    'family': 'Arial',
                    'size': 10,
                },

            }
            alllabels.append(annotations)
    else:
        for data in labelData:
            annotations = {
                'xref': 'x{}'.format(i),
                'yref': 'y{}'.format(i),
                'x': data['x'],
                'y': data['y'],
                'xanchor': 'center',
                'yanchor': 'bottom',
                'text': round(data['x'], 5),
                'overlap': False,
                'showarrow': False,
                'opacity': 1,
                'arrowhead': 0,
                'arrowwidth': 0.5,
                'arrowcolor': 'blue',
                'ax': 0,
                'ay': 0,
                'font': {
                    'family': 'Arial',
                    'size': 10,
                },

            }
            alllabels.append(annotations)
    showAnno = checkOverlap(alllabels,labelWidth,labelHeight)
    return showAnno

def checkOverlap(annotations:list,labelWidth,labelHeight):
    # labelWidth = 2
    # labelHeight = 100000000
    for self in annotations:
        a={'x':self['x'],'y':self['y']}
        self['showarrow']=False
        if (self['overlap']):
            continue;
        for that in annotations:
            if(self['x']!=that['x']):
                b={'x':that['x'],'y':that['y']}
                h=abs(a['x']-b['x'])
                v=abs(a['y']-b['y'])
                if(not (h>labelWidth or h==labelWidth or (h >math.ceil(labelWidth/2)and h<labelWidth  and v>labelHeight+1))):
                    if(that['overlap'] and self['y']<that['y'] and v>math.ceil(labelHeight/2)-2 and h<math.ceil(labelWidth/2)-1 ):
                        self['overlap']=True
                    that['overlap']=True
                # if(that['overlap']==False):
                #     if(h<labelWidth or v<labelHeight or (h<math.ceil(labelWidth/2) and v<labelHeight+1)):
                #         that['overlap']=True
    showedannotations=list(filter(lambda item:not item['overlap'],annotations))
    showedannotations=sorted(showedannotations, key=lambda label: label['x'])
    # for i,v in enumerate(showedannotations):
    #     if(i<len(showedannotations)-1):
    #         leftele = showedannotations[i]
    #         rightele = showedannotations[i + 1]
    #         if (not rightele):break;
    #         betweenAnn=list(filter(lambda item: (item['x'] >leftele['x'])  and (item['x'] < rightele['x']), annotations))
    #         betweenAnnMax = betweenAnn[0]
    #         if (not betweenAnnMax):continue
    #         checkOverlapMiddle(leftele, betweenAnnMax, rightele, labelWidth, labelHeight)
    overlapAnnotations=list(filter(lambda item:not item['overlap'],annotations))
    return overlapAnnotations

def checkOverlapMiddle(leftele, betweenAnnMax, rightele, width, height):
    verticelSpaceMin = 0;
    verticelSpaceMax = 0;
    verticelSpace1 = 0;
    verticelSpace2 = 0
    horizontalSpace1 = 0
    horizontalSpace2 = 0
    moveSpace = 0
    x0=leftele['x']
    x1=betweenAnnMax['x']
    x2=rightele['x']
    y0=leftele['y']
    y1=betweenAnnMax['y']
    y2=rightele['y']
    ymin = max(leftele['y'], betweenAnnMax['y'], rightele['y'])
    ymax = min(leftele['y'], betweenAnnMax['y'], rightele['y'])
    verticelSpaceMin = min(abs(y0 - y1), abs(y1 - y2))
    verticelSpaceMax = max(abs(y0 - y1), abs(y1 - y2))
    horizontalSpace1 = min(abs(x1 - x0), abs(x1 - x2))
    horizontalSpace2 = max(abs(x1 - x0), abs(x1 - x2))
    if(abs(x0-x1)<abs(x1-x2)):
        moveSpace = math.ceil(width / 2) - horizontalSpace1
        verticelSpace1 = abs(y0 - y1)
        verticelSpace2 = abs(y2 - y1)
    else:
        moveSpace = horizontalSpace1 - math.ceil(width / 2)
        verticelSpace2 = abs(y0 - y1)
        verticelSpace1 = abs(y2 - y1)
    if(betweenAnnMax['y'] != ymin and betweenAnnMax['y'] != ymax):
        if(verticelSpaceMin > height and verticelSpaceMax > height):
            betweenAnnMax['overlap'] = False
            betweenAnnMax['showarrow'] = True
            betweenAnnMax['ay'] = -2
            betweenAnnMax['ax'] = moveSpace
        elif (verticelSpace1 > height and horizontalSpace2 > width):
            betweenAnnMax['overlap'] = False
            betweenAnnMax['showarrow'] = True
            if(verticelSpace2 < height and horizontalSpace2 < math.ceil(width / 2) * 3):
                betweenAnnMax['ay'] = -10
            betweenAnnMax['ax'] = moveSpace

    if(betweenAnnMax['y'] == ymin):
        if(verticelSpaceMin > height and horizontalSpace2 > width):
            betweenAnnMax['overlap'] = False
            betweenAnnMax['showarrow'] = True
            betweenAnnMax['ay'] = -1
            betweenAnnMax['ax'] = moveSpace




def getInfoByScanNum(rawfile, scanNum):

    precursor = rawfile.GetFullMSOrderPrecursorDataFromScanNum(
        scanNum, rawfile.GetMSOrderForScanNum(scanNum)-2)[0]

    ionPeaks = rawfile.GetLabelData(scanNum)[0]
    ionPeaks = pd.DataFrame(list(ionPeaks), index=ionPeaks._fields).T
    # ionPeaks['relAbundance'] = ionPeaks.intensity/ionPeaks.intensity.max()*100
    # ionPeaks['_id'] = [bson.ObjectId() for i in range(len(ionPeaks))] #对每个离子新境ID
    ionPeaks = ionPeaks.drop(
        ['intensity', 'resolution', 'baseline', 'noise', 'charge'], axis=1)
    ionPeaks = ionPeaks.to_dict('records')
    return pd.Series([precursor, ionPeaks])





def getAllData(rawfile):
    """
    :Param  rawfile:MSFileReader
    :Returns  data:dataframe
    """

    data = pd.DataFrame(list(range(rawfile.FirstSpectrumNumber,
                                   rawfile.LastSpectrumNumber + 1)), columns=['scanNum'])
    data['msOrder'] = data.apply(
        lambda x: rawfile.GetMSOrderForScanNum(x), axis=1)
    data = data[data['msOrder'] > 1]
    data[['precursor', 'Ion']] = data.apply(
        lambda x: getInfoByScanNum(rawfile, x.scanNum), axis=1)

    return data


def findprecursor(mass, df_ms2):
    # res = {}
    precursor = []
    for i in range(len(mass)):
        temp = set()
        for index, row in df_ms2.iterrows():
            for item in row['Ion']:
                if round(item['mass'], 5) > mass[i]:
                    break
                elif abs(mass[i] - item['mass'])/item['mass']*1e6 <= 5:
                    temp.add(round(row.precursor, 5))
        precursor.append(temp)
        # res.update({mass[i]:precursor})

    t = precursor[0]
    for i in range(len(precursor)):
        t = t & precursor[i]

    precursor = list(t)

    precursor.sort()
    for i in range(len(precursor)-1, 0, -1):
        if abs(precursor[i]-precursor[i-1])/precursor[i-1]*1e6 <= 5:
            del precursor[i]

    return precursor


def getrelIntensity(y):
    df_data = pd.DataFrame(y, columns=['y'])
    df_data['relIntensity'] = df_data.y/df_data.y.max()*100

    return df_data['relIntensity'].to_list()

def getChorInfor(rawfile,x):
    scanNum = rawfile.ScanNumFromRT(x)
    basepeak = dict(rawfile.GetScanHeaderInfoForScanNum(scanNum))['BasePeakMass']
    masslist = rawfile.GetLabelData(scanNum)
    mass_x = list(masslist[0][0])
    mass_y = list(masslist[0][1])
    sort_mass_y = frozenset(mass_y)
    sort_mass_y = sorted(sort_mass_y,reverse=True)
    topmass = []
    for i in range(3):
        topmass.append(mass_x[mass_y.index(sort_mass_y[i])])

    return pd.Series([scanNum,basepeak,topmass])

def getAllChorData(rawfile,startTime,endTime,precursormass='50-1000', massTolerance=5, units=1,msfilter=''):

    df_result = pd.DataFrame()
    rawfile.SetMassTolerance(
    userDefined=True, massTolerance=massTolerance, units=units)
    allChorData = rawfile.GetChroData(startTime,endTime,precursormass,'',msfilter)
    df_result['x'] = allChorData[0][0]
    df_result['y'] = allChorData[0][1]
    df_result['relIntensity'] = df_result.apply(lambda res:res.y/df_result['y'].max()*100,axis=1)
    df_result[['scanNum','basepeak','topmass']] = df_result['x'].apply(lambda x :getChorInfor(rawfile,x))
    
    return df_result



# def getChorData(rawfile, precursor, massTolerance=550, units=0):

#     df_chorData = pd.DataFrame(columns=['X', 'Y', 'relIntensity', 'Name'])
#     # The type of tolerance value (amu = 2, mmu = 0, or ppm = 1)
#     rawfile.SetMassTolerance(
#         userDefined=True, massTolerance=massTolerance, units=units)

#     for i in range(len(precursor)):
#         chorData = rawfile.GetChroData(rawfile.GetStartTime(
#         ), rawfile.GetEndTime(), str(precursor[i]), '', 'Full ms')
#         x = list(chorData[0][0])
#         y = list(chorData[0][1])
#         df_chorData = df_chorData.append(
#             [{'X': x, 'Y': y, 'Name': str(precursor[i])}], ignore_index=True)
#         df_chorData['relIntensity'] = df_chorData['Y'].apply(
#             lambda x: getrelIntensity(x))

#     return df_chorData


def getSpectrumData(rawfile, time):
    scanNum = rawfile.ScanNumFromRT(time)
    msfilter = rawfile.GetFilterForScanNum(scanNum)
    # data = rawfile.GetLabelData(scanNum)
    # data = rawfile.GetMassFromScanNum(scanNum)
    data = rawfile.GetMassListFromScanNum(scanNum)
    x = data[0][0]
    y = data[0][1]
    # charge = data[0][5]
    max_y = max(y)
    relIntensity = []
    for i in range(len(y)):
        relIntensity.append(y[i]/max_y*100)
    # print(type(charge))
    return x, y, relIntensity,msfilter
    # return {
    #     x:list(x),
    #     y:list(y),
    #     charge:list(charge),
    #     relIntensity:relIntensity
    # }

