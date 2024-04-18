# -*- coding:utf-8 -*-
#@Time : 2020/10/29 11:13
#@Author: zgz
#@File : EventTool.py

"""
事件用到的一些方法
"""
import numpy as np
import copy
from PICore.common import peak
from PICore.DetectMethod import Trad_Integral as Ti

#最后峰的一些面积计算方法
def calculateResult(pre_result, rawData, timeData, inteval):
    resclusters = []
    tempClusterID = 0
    tempPeakID = 0
    for group in pre_result:
        # 判断起点和终点是否一致 如果一致则排除掉该峰

        bslType = group["bslType"]
        tempcluster = peak.PeakCluster(clusterID=tempClusterID)
        tempbslcoefs = group["bsl"]
        # python3 删除了has_key
        if ('bslY') in group:
            tempbsl = peak.Baseline(startpoint=int(tempbslcoefs[0]),
                                    endpoint=int(tempbslcoefs[1]),
                                    rawData=rawData,
                                    inteval=inteval,
                                    bslType = bslType,
                                    startY=group["bslY"][0], endY=group["bslY"][1]
                                    ) #传入的是频率
        else:
            tempbsl = peak.Baseline(startpoint=int(tempbslcoefs[0]),
                                    endpoint=int(tempbslcoefs[1]),
                                    rawData=rawData,
                                    inteval=inteval,
                                    bslType=bslType,
                                    startY=None, endY=None
                                    )
        tempcluster.setBaseline(tempbsl)
        bsl_data =tempbsl.generatebslData(rawData,inteval) #一个峰簇里面
        q = len(group["peak"])
        i = 0
        for currpeak in group["peak"]:
            tempPeak = peak.Peak(clusterID = int(tempClusterID), peakID = int(tempPeakID),
                                 startPoint = int(currpeak[0]), endPoint = int(currpeak[-1]), apexPoint = int(currpeak[1]),
                                 baseline = tempbsl)
            tempPeak.peakCalculation(rawData, timeData, inteval) #传入的是时间间隔
            if q == 1:
                tempPeak.PeakType = "BB" #这种情况时单峰 则基线和垂直线的交叉点不变 是None
            else:
                if i == 0 and q > 1:
                    tempPeak.PeakType = "BV"
                elif i == q - 1:
                    tempPeak.PeakType = "VB"
                else:
                    tempPeak.PeakType = "VV"
                tempPeak.findBLY(bsl_data, int(tempbslcoefs[0]))# 这种情况下是融合峰, 则变化基线和垂直线的交叉点
            tempPeakID = tempPeakID + 1
            i += 1
            tempcluster.Peaks.append(tempPeak)
        tempClusterID += 1
        resclusters.append(tempcluster)
    # 循环体结束
    res = peak.integralResult(resclusters)
    res.UpdatePercentArea()
    return res.ExportResDict(inteval)


def djuge(result,bslType,preResult):
    indexSub = []
    for i in range(len(result)-1):
        if result[i][-1] != result[i+1][0]:
            indexSub.append(i + 1)
            tempgroup = {"peak": [], "bsl": None, "bslType": bslType}
            tempgroup['peak'].append(result[i + 1])
            tempgroup['bsl'] = [result[i + 1][0], result[i + 1][-1]]
            preResult.append(tempgroup)
    return indexSub,preResult

#批量删除
def dealResult(Preresult,indexFlag):
    Preresult = [Preresult[i] for i in range(len(Preresult)) if(i not in indexFlag)] #批量删除列表元素
    return Preresult

#将处理事件中的融合峰根据实际情况断开 并添加到事件之外的结果中
def solveFusedPeakEvent(EventResult,preResult):
    for group in EventResult['peak']:
        bslType = EventResult['bslType']
        # 单峰
        if len(group) == 1:
            tempgroup = {"peak": [], "bsl": None, "bslType": bslType}
            tempgroup['peak'].append(group[0])  #这种方式形成了二维的
            tempgroup['bsl'] = [group[0][0],group[0][-1]]
            preResult.append(tempgroup)
        else:
                indexSub, preResult = djuge(group,bslType,preResult)
                if len(indexSub) >= 1:
                    temp = dealResult(group, indexSub)
                    indexSub, preResult = djuge(temp, bslType, preResult)
                    if len(indexSub)==0:
                        tempgroup = {"peak": [], "bsl": None, "bslType": bslType}
                        tempgroup['peak'] = temp  #这种也形成了二维的 因为temp 是就两个
                        tempgroup['bsl'] = [temp[0][0],temp[-1][-1]]
                        preResult.append(tempgroup)
                    #当融合峰没有断开的点时 则全部加进去
                else:
                    tempgroup = {"peak": [], "bsl": None, "bslType": bslType}
                    tempgroup['peak'] = group
                    tempgroup['bsl'] = [group[0][0], group[-1][-1]]
                    preResult.append(tempgroup)

    return preResult #将符合条件的和断开的事件结果进行整合

def djugefusedpeak(group):
    """
解决事件结果中存在基线穿透的问题
    :param group: 要断开的时间结果
    :return: 返回的是融合峰 和单峰的结果
    """
    difpeakList = []
    fusepeak =[]
    indexFlag = []
    sameList = []
    Flag = [0]*len(group) #生成和group一样大小的标签列表
    #这种情况下 如果第一个峰不和第二峰是融合峰的话则第一个结果不会被处理，如是融合峰，则会被处理，添加到融合峰中
    for i  in range(len(group)-1,0,-1):
        if group[i][0] != group[i-1][-1]:
            if Flag[i]==0:
                difpeakList.append(group[i])
            else:
                continue
        else:
            # 如果当前的i适合indexFlag中最后索引是连续的 则在sameList中接着添加， 如果不是则另起炉灶 创建一个新的融合峰
            if len(indexFlag) ==0: # 为空时表示刚开始
                if Flag[i] == 0:
                    sameList.append(group[i])
                    indexFlag.append(i)
                if Flag[i - 1] == 0:
                    sameList.append(group[i-1])
                    indexFlag.append(i-1)
            elif i == indexFlag[-1]: #表示连续的一段是融合峰时
                if Flag[i] == 0:
                    sameList.append(group[i])
                    indexFlag.append(i)
                if Flag[i - 1] == 0:
                    sameList.append(group[i - 1])
                    indexFlag.append(i - 1)
            else:
                fusepeak.append(sameList)
                sameList = []
                if Flag[i] == 0:
                    sameList.append(group[i])
                    indexFlag.append(i)
                if Flag[i - 1] == 0:
                    sameList.append(group[i - 1])
                    indexFlag.append(i - 1)
            Flag[i - 1] = 1
#最后结果的特殊处理
    if group[0][-1] != group[1][0]:
        difpeakList.append(group[0]) #二维的列表
    # 如果 fusepeak
    if len(fusepeak) ==0 and len(sameList) != 0:
        fusepeak.append(sameList) #三维的列表
    elif len(fusepeak) !=0 and len(sameList) !=0:
        tempfusePeakList = []
        for tempfusePeak in fusepeak:
            for item in tempfusePeak:
                tempfusePeakList.append(item)
        temp = [x for x in sameList if x in tempfusePeakList] #如果有元素 说明samelist的元素在fusepeak中，也说明此时sameList没有重新更新 如果没有则说明sameList有更新，但没有添加到fusepeak中
        if len(temp) ==0:
            fusepeak.append(sameList)
        else:
            pass
    return fusepeak,difpeakList



def DealHorizontalBSL(data,EventResult,preResult):
    # 事件结果存在的情况：1、只有一个单峰；2、只有一个融合峰；3、既有单峰，又有融合峰
    # 当事件结果中只有一个峰（单峰或者融合峰）时,基线的起终点就在峰上,这个峰可能是融合峰或者是单峰，做这步的操作是由于EventResult存储的峰是2维度的[[[]]]
    bslType = EventResult['bslType']  #
    bs = EventResult['bsl'][0] #事件基线的开始点
    be = EventResult['bsl'][-1] #时间基线的终止点
    bslData = HorizonBaseline(data, bs, be, bslType) #根据当前的事件结果生成临时的基线
    for item in range(len(EventResult['peak'])):
        group = EventResult['peak'][item]
        # 某个结果是单峰时
        if len(group) ==1:
            tempgroup = {"peak": [], "bsl": None, "bslType": bslType,"bslY":None}
            tempgroup['peak'].append(group[0])
            tempgroup['bsl'] = [group[0][0],group[0][-1]]
            tempgroup["bslY"] = [bslData[0],bslData[-1]] #用于存放基线不在原始谱图上 针对水平基线类的事件
            preResult.append(tempgroup)
        # 结果是融合峰事
        else:
            #当结果组中只有两个峰时
            if len(group)==2:
                # 当结果中两个峰不是融合峰,则断开
                if group[0][-1] != group[1][0]:
                    for i in range(len(group)):
                        tempgroup = {"peak": [], "bsl": None, "bslType": bslType, "bslY": None}
                        tempgroup['peak'].append(group[i])
                        tempgroup['bsl'] = [group[i][0], group[i][-1]]
                        tempgroup["bslY"] = [bslData[0], bslData[-1]]  # 用于存放基线不在原始谱图上 针对水平基线类的事件
                        preResult.append(tempgroup)
                else:
                    # 否则则当做是一个融合峰
                    tempgroup = {"peak": [], "bsl": None, "bslType": bslType, "bslY": None}
                    tempgroup['peak'] = group
                    tempgroup['bsl'] = [group[0][0], group[-1][-1]]
                    tempgroup["bslY"] = [bslData[0], bslData[-1]]  #
                    preResult.append(tempgroup)
            else:
                # 当结果组中峰的个数大于两个时
                fusepeak,difpeakList = djugefusedpeak(group) #fusepeak 是倒叙的
                # fusepeak = sorted(fusepeak, key=lambda x: x[0]) #将融合峰里面的排序是
                if len(fusepeak) !=0:
                    for i in range(len(fusepeak)):
                        temp = fusepeak[i]
                        tempgroup = {"peak": [], "bsl": None, "bslType": bslType, "bslY": None}
                        tempgroup['peak'] = temp #融合直接赋值即可
                        tempgroup['bsl'] = [temp[-1][0],temp[0][-1]]
                        tempgroup["bslY"] = [bslData[0], bslData[-1]]  #
                        preResult.append(tempgroup)
                if len(difpeakList) !=0:
                    for i in  range(len(difpeakList)):
                        temp = difpeakList[i]
                        tempgroup = {"peak": [], "bsl": None, "bslType": bslType, "bslY": None}
                        tempgroup['peak'].append(temp) #单峰则添加，作为融合峰一样的二维列表，单只有一个峰而已
                        tempgroup['bsl'] = [temp[0], temp[-1]]
                        tempgroup["bslY"] = [bslData[0], bslData[-1]]  #
                        preResult.append(tempgroup)
    return preResult

#处理事件中基线穿越的问题
def solveSingleEventbsl(EventResult, rawData):
    bslType = EventResult["bslType"]
    tempBslstart = EventResult['peak'][0][0][0]
    tempBslend = EventResult['peak'][-1][-1][-1]
    tempbsl = HorizonBaseline(rawData, tempBslstart, tempBslend, bslType)  #水平基线 有两种情况 一种是正向水平 一种是反向水平的设置
    for groupi in EventResult['peak']:
        tempReuslt = []
        for groupj in groupi:
            # 将处于基线下的峰删掉 不能直接动态的删除 这样子会根据原来的位置跳过某些数据不检测
            if rawData[groupj[1]] <tempbsl[0]:
                tempReuslt.append(groupj)
            #判断一个峰的前半部分的起点是否小于基线 如果小于则寻找大于等于基线的第一个点
            elif rawData[groupj[0]] < tempbsl[0] and rawData[groupj[-1]] < tempbsl[0]:
                start = groupj[0]
                mid = groupj[1]
                temprawData = rawData[start:mid+1]
                indexLine = np.linspace(start, mid, mid - start + 1)
                index = np.where(temprawData >=  tempbsl[0])  # 返回的是元组类型 元组存放一个数组 我们默认采用第一个
                dataindex = indexLine[index[0][0]]
                groupj[0] =int(dataindex-1)
                # 然后依次判断一个峰的后半部分是否小于基线 如果小于基线 则寻找小于基线的第一个点
                mid = groupj[1]
                end = groupj[2]
                temprawData = rawData[mid:end+1]
                indexLine = np.linspace(mid, end, end - mid + 1)
                index = np.where(temprawData <= tempbsl[0])  # 返回的是元组类型 元组存放一个数组 我们默认采用第一个
                dataindex = indexLine[index[0][0]]
                groupj[2] = int(dataindex-1)
            elif rawData[groupj[-1]] < tempbsl[0]:
                mid = groupj[1]
                end = groupj[2]
                temprawData = rawData[mid:end+1]
                indexLine = np.linspace(mid, end, end - mid + 1)
                index = np.where(temprawData <= tempbsl[0])  # 返回的是元组类型 元组存放一个数组 我们默认采用第一个
                dataindex = indexLine[index[0][0]]
                groupj[2] = int(dataindex - 1)
            elif rawData[groupj[0]] < tempbsl[0]:
                start = groupj[0]
                mid = groupj[1]
                temprawData = rawData[start:mid+1]
                indexLine = np.linspace(start, mid, mid - start + 1)
                index = np.where(temprawData >= tempbsl[0])  # 返回的是元组类型 元组存放一个数组 我们默认采用第一个
                dataindex = indexLine[index[0][0]]
                groupj[0] = int(dataindex - 1)
        if len(tempReuslt) != 0:
            for item in tempReuslt:
                groupi.remove(item)
     # 当前事件段内的峰都在基线下 则返回空值  否则就是返回正常的值
    temp =copy.copy(EventResult['peak'])
    for groupi in temp:
        if len(groupi)==0 :
            EventResult['peak'].remove(groupi)
    if len(EventResult['peak']) == 0:
        return None
    else:
        EventResult['bsl'][-1] = EventResult['peak'][-1][-1][-1]
        return EventResult

# 对事件之外的一些峰进行基线的处理。因为事件的开始和结束时间可能是在融合峰中，此时断开了融合，那么剩余的峰的基线可能会变化，所以也要进行处理一遍
def solveBaseline(EventResult, rawData):
    result = []
    # 从前往后扫描，出现基线穿越则断开，开始新的峰簇
    for group in EventResult:
        # 对事件之外来判断是否存在基线穿越时 不对之前的水平基线事件判断是否存在基线穿越的情况
        bslType = group["bslType"]
        # 单峰
        if bslType == "auto":
            if len(group['peak']) == 1:
                result.append(group)
            else:
                # 融合峰：逐个峰扫描，若两峰之间最低点低于基线，则从最低点处断开基线
                # 断开基线前为一个融合峰，后面为另一个融合峰
                tempBslstart = group['peak'][0][0]
                tempBslend = group['peak'][-1][-1]
                tempbsl = Ti.baseline(rawData, tempBslstart, tempBslend)
                tempgroup = {"peak":[], "bsl":None,"bslType":bslType}
                tempstartcheck = tempBslstart
                tempbslcoef = [tempBslstart, tempBslend]
                for i in range(len(group['peak'])-1):
                    if group['peak'][i][-1] == group['peak'][i+1][0]:
                        nextLow = group['peak'][i][-1]
                    else:
                        nextLow = group['peak'][i][-1] + (
                        rawData[group['peak'][i][-1]:group['peak'][i + 1][0]]).argmin()
                        group['peak'][i][-1] = nextLow
                        group['peak'][i+1][0] = nextLow
                    tempgroup["peak"].append(group['peak'][i])
                    if rawData[nextLow] < tempbsl[nextLow-tempstartcheck]:
                        tempgroup["bsl"] = [tempBslstart,nextLow]
                        result.append(tempgroup)
                        tempgroup = {"peak":[], "bsl":None,"bslType":bslType}
                        tempBslstart = group['peak'][i+1][0]
                        tempbsl = Ti.baseline(rawData, tempBslstart, tempBslend)
                        tempstartcheck = tempBslstart
                        tempbslcoef = [tempBslstart, tempBslend]
                tempgroup["peak"].append(group['peak'][-1])
                tempgroup["bsl"] = tempbslcoef
                result.append(tempgroup)
        else:
            result.append(group)
    return result


def selectResult(startTime,endTime,Preresult):
    #先找出在此时间段的峰,在此段峰中进行操作,在
    curres =[]
    indexFlag = []
    Preresult = sorted(Preresult, key=lambda x: x["bsl"][0])  # 按照基线的起始点排序
    for i in range(len(Preresult)):
    # 当事件的开始时间段在 一个融合峰中,则需要将满足条件的后部分添加进来，此时需要将原始的融合峰断开
        if len(Preresult[i]['peak'])<=1:
            if (Preresult[i]['peak'][0][0] >= startTime or Preresult[i]['peak'][0][0] <= startTime)and Preresult[i]['peak'][0][2]> startTime and (Preresult[i]['peak'][0][2] <= endTime or Preresult[i]['peak'][0][0] <= endTime) :
                curres.append(Preresult[i])
                indexFlag.append(i)
        else:
            if Preresult[i]['peak'][0][1] >= startTime and Preresult[i]['peak'][-1][2] <= endTime:
                curres.append(Preresult[i])
                indexFlag.append(i)
            else:
                 indexSub = []
                 for j in range(len(Preresult[i]['peak'])):
                     if (Preresult[i]['peak'][j][0] >= startTime and Preresult[i]['peak'][j][2] <= endTime) or (Preresult[i]['peak'][j][1] >= startTime and Preresult[i]['peak'][j][2] <= endTime):
                         indexSub.append(j)
                     elif (Preresult[i]['peak'][j][1] <= startTime and Preresult[i]['peak'][j][2] >= startTime):
                         indexSub.append(j)
                     elif (Preresult[i]['peak'][j][1] <= endTime and Preresult[i]['peak'][j][2] >= endTime):
                         indexSub.append(j)
                     elif (Preresult[i]['peak'][j][0] <= endTime and Preresult[i]['peak'][j][1] >= endTime):
                         indexSub.append(j)
                    #当整个融合峰都在这个事件时间段内
                 if len(Preresult[i]['peak'])== len(indexSub):
                     curres.append(Preresult[i])
                     indexFlag.append(i)
                     #当融合峰部分在事件时间段内  同时要将融合峰断开
                 elif len(indexSub) >=1:
                     if indexSub[0] != 0 and indexSub[-1] != len(Preresult[i]['peak'])-1 :
                         tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
                         tempgroup['peak'].extend(Preresult[i]['peak'][indexSub[-1]+1:]) #取值时包头不包尾
                         tempbslcoef = [tempgroup['peak'][0][0],tempgroup['peak'][-1][-1]]
                         tempgroup['bsl'] = tempbslcoef
                         del (Preresult[i]['peak'][indexSub[-1]+1:])
                         Preresult.append(tempgroup)
                         #处理在事件内的值
                     tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
                     for item in range(len(indexSub)):
                         tempgroup['peak'].append(Preresult[i]['peak'][indexSub[item]])

                     tempbslcoef = [tempgroup['peak'][0][0], tempgroup['peak'][-1][2]]
                     tempgroup['bsl'] = tempbslcoef
                     Preresult[i]['peak'] = dealResult(Preresult[i]['peak'], indexSub) # 批量删除
                     Preresult[i]['bsl'][0] = Preresult[i]['peak'][0][0]
                     Preresult[i]['bsl'][1] = Preresult[i]['peak'][-1][2]
                     curres.append(tempgroup)


    Preresult = dealResult(Preresult, indexFlag) #删除事件内的峰结果
    return curres,Preresult

#所有事件公用的方法 判断当前事件的时间段同时在一个峰的前面或者一个峰的后面
def getStartandEnd (startTime,endTime,selectRes,Preresult):
    if len(selectRes[0]['peak']) ==1 and len(selectRes) ==1 :
        if (startTime > selectRes[0]['peak'][0][0] and endTime < selectRes[0]['peak'][0][1]) or  (startTime > selectRes[0]['peak'][0][1] and endTime < selectRes[0]['peak'][0][2]):
            tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
            tempgroup['peak'].append(selectRes[0]['peak'][0])
            tempbslcoef = [selectRes[0]['peak'][0][0], selectRes[0]['peak'][0][2]]
            tempgroup['bsl'] = tempbslcoef
            Preresult.append(tempgroup)
            del (selectRes[0])
    return Preresult, selectRes

def HorizonBaseline(data, bs, be,bslType):
    """
    生成临时基线数据
    :param data: array_like(原数据)
    :param x1: int(基线起点)
    :param x2: int(基线终点)
    :return bsl: array_like()
    """
    if bslType =="Horizon":
        y1 = data[bs]
    else:
        y1 = data[be]
    num = be - bs + 1
    bsl = np.linspace(y1, y1, num)
    return bsl

def generatebasline(data, bs, be):
    """
    生成基线数据
    :param data: array_like(原数据)
    :param bs: int(基线起点)
    :param be: int(基线终点)
    :return bsl: array_like(基线上的点数组)
    """
    y1 = data[bs]
    y2 = data[be]
    num = be - bs + 1
    bsl = np.linspace(y1, y2, num)  # linspace相当于根据起止点生成这条直线上的数据
    return bsl
