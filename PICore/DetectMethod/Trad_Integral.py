# import sys
# import os
# sys.path.append('../')

import numpy as np
from scipy import signal
from PICore.common import peak
from PICore.common.base1 import peakwidth

def __autoThreshold(rawData, detectWidth, dataFreq): 
    # 默认阈值，当输入阈值为Null时使用
    # 二阶导数噪声 * 峰宽参数
    # secondDerivativeFilted = signal.medfilt(-np.diff(rawData, 2))[:]
    # return np.std(secondDerivativeFilted)*detectWidth*3
    # 根据采样频率来进行设置自动阈值的点数设置
    pointNum = int(20*(dataFreq))
    secondDerivativeFilted = signal.medfilt(-np.diff(rawData, 2))[-pointNum:]
    return  np.std(secondDerivativeFilted)*detectWidth*3

def __preFilter(rawData):
    # 前置滤波器
    b, skiprow=signal.butter(8, 0.9, 'lowpass')  
    rawData=signal.filtfilt(b, skiprow, rawData)
    return rawData

def __autoPeakWidth(rawData,dataFreq):
    pw_index,wn = peakwidth(rawData)
    return pw_index,int(wn / dataFreq)

def __bundleData1(rawData,dataFreq,detectWidth):
    # 进行原始数据点分组，根据输入峰宽参数、采样频率计算分组间斜率
    bundleLength = int(detectWidth * dataFreq // 15)
    if bundleLength < 1:
        bundleLength = 1
    rawLength = len(rawData)
    # 若最后一组补不足，则去掉不足的部分    
    residual = rawLength % bundleLength
    if isinstance(rawData, np.ndarray):
        rawData = rawData.tolist()
    if residual != 0:
        for i in range(residual):
            rawData.pop()
    bundledRawData = np.array(rawData).reshape(len(rawData)//bundleLength,bundleLength)
    avgBundle = bundledRawData.mean(axis=1)
    # 相邻两组间斜率
    bundledSlope = []
    for i in range(len(avgBundle)-1):
        bundledSlope.append((avgBundle[i+1]-avgBundle[i])/(bundleLength / dataFreq))
    # 返回分组、组间斜率、每组点数
    return bundledRawData, bundledSlope, bundleLength

def __bundleData(rawData,timeInterval,detectWidth):
    # 进行原始数据点分组，根据输入峰宽参数、采样频率计算分组间斜率
    bundleLength = int(detectWidth /timeInterval // 15)
    if bundleLength < 1:
        bundleLength = 1
    rawLength = len(rawData)
    # 若最后一组补不足，则去掉不足的部分
    residual = rawLength % bundleLength
    if isinstance(rawData, np.ndarray):
        rawData = rawData.tolist()
    if residual != 0:
        for i in range(residual):
            rawData.pop()
    bundledRawData = np.array(rawData).reshape(len(rawData)//bundleLength,bundleLength)
    avgBundle = bundledRawData.mean(axis=1)
    # 相邻两组间斜率
    bundledSlope = []
    for i in range(len(avgBundle)-1):
        bundledSlope.append((avgBundle[i+1]-avgBundle[i])/(bundleLength *timeInterval))
    # 返回分组、组间斜率、每组点数
    return bundledRawData, bundledSlope, bundleLength

# 进行初步峰检测
# 此处prepeak 为[峰起点索引、峰顶点索引、峰终点索引]形式
def __findPrePeaks(bundledRawData, bundledSlope, bundleLength, threshold,startPoint=0):
    prepeakStartIdx = []
    prepeakTopIdx = []
    prePeaks = []
    liftoffFlag = False
    touchdownFlag = False
    for i in range(len(bundledSlope)-1):
        # 情况1：找峰起点
        if not liftoffFlag and not touchdownFlag:
            # 1.1 宽度较小的峰（分组中，出现连续的一个>起点阈值，一个<终点阈值斜率）
            if bundledSlope[i] >= threshold and bundledSlope[i+1] < 0:
                tr = bundledRawData[i].argmin()
                startidx = i*bundleLength + tr + startPoint
                prepeakStartIdx.append(startidx)
                # 从当前起点最低的出发，找局部最高值，避免出现由于分组原因导致的本组起点前的点高于本峰的最高值
                topidx = startidx + np.concatenate([bundledRawData[i][tr:],bundledRawData[i+1]]).argmax()
                prepeakTopIdx.append(topidx)
                liftoffFlag = False
                touchdownFlag = True
            else:
                avgslope = np.abs(bundledSlope[i] + bundledSlope[i+1])/2
                # 常规情况，连续两个一阶导数均值 > 阈值
                if avgslope >= threshold:
                    startidx = i*bundleLength + bundledRawData[i].argmin() + startPoint
                    prepeakStartIdx.append(startidx)
                    # prepeakStartIdx.append(i)
                    liftoffFlag = True
                    touchdownFlag = False
                else:
                    continue
        # 情况2：找峰顶点
        elif liftoffFlag and not touchdownFlag:
            if bundledSlope[i] >=0 and bundledSlope[i+1] < 0:
                topidx = i*bundleLength + np.concatenate(bundledRawData[i:i+2]).argmax() + startPoint
                prepeakTopIdx.append(topidx)
                # prepeakTopIdx.append(i)
                liftoffFlag = False
                touchdownFlag = True
            else:
                continue
        # 情况3：找峰终点
        elif touchdownFlag and not liftoffFlag:
            if bundledSlope[i] < -threshold:
                # 3.1 突然出现的斜率变正，从当前顶点到下一个变正的区间内找最低点作为峰
                if bundledSlope[i+1] >= 0:
                    tr = prepeakTopIdx.pop() - startPoint
                    z = tr // bundleLength
                    c = tr % bundleLength
                    endidx = tr + np.concatenate([bundledRawData[z][c:], np.concatenate(bundledRawData[z+1:i+2])]).argmin() + startPoint
                    prePeaks.append([prepeakStartIdx.pop(), tr + startPoint, endidx])
                    touchdownFlag = False
                    liftoffFlag = False
                elif bundledSlope[i+1] < threshold:
                    continue
            else:
                # 连续两个斜率 > -阈值，从当前顶点到下一个变正的区间内找最低点作为峰
                if bundledSlope[i+1] >= -threshold:
                    tr = prepeakTopIdx.pop() - startPoint
                    z = tr // bundleLength
                    c = tr % bundleLength
                    endidx = tr + np.concatenate([bundledRawData[z][c:], np.concatenate(bundledRawData[z+1:i+2])]).argmin() + startPoint
                    prePeaks.append([prepeakStartIdx.pop(), tr + startPoint, endidx])                    
                    liftoffFlag = False
                    touchdownFlag = False
                else:
                    continue
    return prePeaks

def baseline(data, bs, be):
    """
    生成临时基线数据
    :param data: array_like(原数据)
    :param x1: int(基线起点)
    :param x2: int(基线终点)
    :return bsl: array_like()
    """
    y1 = data[bs]
    y2 = data[be]
    num = be - bs + 1
    bsl = np.linspace(y1, y2, num)
    return bsl

def __solveFusedPeaks(prePeaks):
    # 融合峰处理
    result = []
    tempcluster = [prePeaks[0]]
    for i in range(len(prePeaks)-1):
        # 判断前后两峰间距离与两峰宽比值
        w3 = prePeaks[i+1][0] - prePeaks[i][2]  
        w1 = prePeaks[i][2] - prePeaks[i][0]
        w2 = prePeaks[i+1][2] - prePeaks[i+1][0]
        if max(w1,w2) >= 3*w3:
            tempcluster.append(prePeaks[i+1])
        else:
            result.append(tempcluster)
            tempcluster = [prePeaks[i+1]]
    result.append(tempcluster)
    return result


def __solveBaseline(fusedPeaks, rawData):
    result = []
    # 从前往后扫描，出现基线穿越则断开，开始新的峰簇
    for group in fusedPeaks:
        tempBslstart = group[0][0]
        tempBslend = group[-1][-1]
        tempbsl = baseline(rawData, tempBslstart, tempBslend)
        if len(group) == 1:
            # 单峰
            # 确定起终点前后不会与基线交叉
            temppeakpos = rawData[group[0][0] : group[0][2] + 1] - tempbsl[0 : group[0][2] - group[0][0] + 1]
            newst = group[0][0] + temppeakpos[0 : group[0][1] - group[0][0] + 1].argmin()
            newed = group[0][1] + temppeakpos[group[0][1] - group[0][0] : group[0][2] - group[0][0] + 1].argmin()
            group[0][0], group[0][2] = newst, newed
            res = {"peak":[group[0]], "bsl":[newst, newed],"bslType":"auto"}
            result.append(res)
        else:
            # 融合峰：逐个峰扫描，若两峰顶点之间出现低于基线的点，则从基线下最低点处断开基线
            # 断开基线前为一个融合峰，后面为另一个融合峰
            # 首先判断起点处是否出现交叉
            temppeakpos = rawData[group[0][0] : group[1][2] + 1] - tempbsl[0 : group[1][2] - group[0][0] + 1]
            newst = group[0][0] + temppeakpos[0 : group[0][1] - group[0][0] + 1].argmin()
            group[0][0] = newst
            tempBslstart = group[0][0]
            tempbsl = baseline(rawData, tempBslstart, tempBslend)
            tempgroup = {"peak":[], "bsl":None,"bslType":"auto"}
            tempstartcheck = tempBslstart
            temppeakpos = rawData[tempstartcheck : group[-1][-1] + 1] - tempbsl[0 : group[-1][-1] - tempstartcheck + 1]
            tempbslcoef = [tempBslstart, tempBslend]
            for i in range(len(group)-1):
                nextLow = group[i][1] + (temppeakpos[group[i][1] - tempstartcheck : group[i+1][1] - tempstartcheck + 1]).argmin()
                group[i][-1] = nextLow
                group[i+1][0] = nextLow
                tempgroup["peak"].append(group[i])
                if temppeakpos[nextLow - tempstartcheck] < 0:
                    tempgroup["bsl"] = [tempBslstart,nextLow]
                    result.append(tempgroup)
                    tempgroup = {"peak":[], "bsl":None,"bslType":"auto"}
                    tempBslstart = group[i+1][0]
                    tempbsl = baseline(rawData, tempBslstart, tempBslend)
                    tempstartcheck = tempBslstart
                    temppeakpos = rawData[tempstartcheck : group[-1][-1] + 1] - tempbsl[0 : group[-1][-1] - tempstartcheck + 1]
                    tempbslcoef = [tempBslstart, tempBslend]
            tempgroup["peak"].append(group[-1])
            tempgroup["bsl"] = tempbslcoef
            result.append(tempgroup)
    return result


def __calculateResult(pre_result, rawData, timeData, timeInterval):
    resclusters = []
    tempClusterID = 0
    tempPeakID = 0
    for group in pre_result:
        tempcluster = peak.PeakCluster(clusterID=tempClusterID)
        tempbslcoefs = group["bsl"]
        tempbsl = peak.Baseline(startpoint=int(tempbslcoefs[0]),
                                endpoint=int(tempbslcoefs[1]),
                                curvetype="linear",
                                curvepoints=[rawData[tempbslcoefs[0]], rawData[tempbslcoefs[1]]])
        tempcluster.setBaseline(tempbsl)
        q = len(group["peak"])
        i = 0
        for currpeak in group["peak"]:
            tempPeak = peak.Peak(clusterID = int(tempClusterID), peakID = int(tempPeakID), 
                                 startPoint = int(currpeak[0]), endPoint = int(currpeak[-1]), apexPoint = int(currpeak[1]),
                                 baseline = tempbsl)
            tempPeak.peakCalculation(rawData, timeData, timeInterval)
            if q == 1:
                tempPeak.PeakType = "BB"
            else:
                if i == 0 and q > 1:
                    tempPeak.PeakType = "BV"
                elif i == q-1:
                    tempPeak.PeakType = "VB"
                else:
                    tempPeak.PeakType = "VV"
            tempPeakID = tempPeakID + 1
            i += 1
            # res.append(json.loads(json.dumps(tempPeak, default=lambda o: o.__dict__, sort_keys=True, indent=4)))
            tempcluster.Peaks.append(tempPeak)
        tempClusterID  += 1
        resclusters.append(tempcluster)
    res = peak.integralResult(resclusters)
    res.UpdatePercentArea()
    return res.ExportJsonStr()

# def __adjustPrepeaks(rawData, prePeaks):
#     # 对峰进行调整，确保基线在本峰内不穿透
#     for peak in prePeaks:
#         peak[0] = peak[0] + np.argmin(rawData[peak[0]:peak[1]])
#         peak[2] = peak[1] + np.argmin(rawData[peak[1]:peak[2]])
#     return prePeaks

def getTradPreResult(rawData, timeData, timeInterval, startTime, endTime, threshold = None, detectWidth = None):
    # dataFreq = int(round(1 / timeInterval))
    dataFreq = ((1 / timeInterval))
    if detectWidth is None:
        # 根据采样频率来获取后面一段数据来进行计算阈值
        pw_index,detectWidth = __autoPeakWidth(rawData, dataFreq)
    if threshold is None:
        threshold = __autoThreshold(rawData, detectWidth, dataFreq)

    # start = int((startTime - timeData[0])//timeInterval)
    # end = int((endTime - timeData[0])//timeInterval)
    start = int(startTime // timeInterval)
    end = int(endTime // timeInterval)
    rawDataSliced = rawData[start:end+1]
    # bundledRawData, bundledSlope, bundleLength = __bundleData(rawDataSliced, dataFreq, detectWidth)
    bundledRawData, bundledSlope, bundleLength = __bundleData(rawDataSliced, timeInterval, detectWidth)
    prePeaks = __findPrePeaks(bundledRawData, bundledSlope, bundleLength, threshold, start)
    if len(prePeaks) ==0:
        err= -60
        return err, None
    fusedPeaks = __solveFusedPeaks(prePeaks)
    pre_result = __solveBaseline(fusedPeaks, rawData)
    return 0, pre_result

def __traditonalIntegral(rawData, timeData, timeInterval, threshold = 0, detectWidth = 15):
    dataFreq = 1 / timeInterval
    if threshold is None:
        threshold = __autoThreshold(rawData, detectWidth, dataFreq)
    bundledRawData, bundledSlope, bundleLength = __bundleData(rawData, dataFreq, detectWidth)
    prePeaks = __findPrePeaks(bundledRawData, bundledSlope, bundleLength, threshold)
    # prePeaks = __adjustPrepeaks(rawData, prePeaks)
    fusedPeaks = __solveFusedPeaks(prePeaks)
    pre_result = __solveBaseline(fusedPeaks, rawData)
    result = __calculateResult(pre_result, rawData, timeData, timeInterval)
    return result

