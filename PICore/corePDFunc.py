# import sys
# import os
# sys.path.append('../')


from PICore.common.peak import *
import PICore.EventFramework as Ef
from PICore.Events import EventTool as ET
from PICore import ParamCheck as PC


def AdjustBSL(Presult,rawData,inteval):
    """
    讲过滤后的结果重新调整基线
    :param Presult: 过滤后的结果
    :return: 调整基线后的结果
    """
    for index in range(len(Presult) - 1, -1, -1):
        # 当只有个峰的时候 检查一下基线的起终点是否和单峰的起终点是否一致
        if len(Presult[index]['Peaks'])==1:
            PeakStartIdx = Presult[index]['Peaks'][0]['StartIdx']
            PeakEndIdx = Presult[index]['Peaks'][0]['EndIdx']
            if PeakStartIdx==Presult[index]['Baseline']['startIdx'] and PeakEndIdx==Presult[index]['Baseline']['endIdx']:
                pass
            else:
                # 如果不是水平情况下 则让基线的起终点等于峰的起终点 并重新计算一下斜率和漂移
                if Presult[index]['bslType']=="auto":
                    Presult[index]['Baseline']['startIdx'] = PeakStartIdx
                    Presult[index]['Baseline']['endIdx'] = PeakEndIdx
                    Presult[index]['Baseline']["slope"] = rawData[PeakEndIdx]-rawData[PeakStartIdx]/PeakEndIdx-PeakStartIdx #更新斜率 self.drift/((self.endPoint-self.startPoint)*inteval/3600)
                    Presult[index]['Baseline']["drift"] = (rawData[PeakEndIdx]-rawData[PeakStartIdx])/((PeakEndIdx-PeakStartIdx)*inteval/3600)
                else:
                # 否则则只要更新基线的起终点
                    Presult[index]['Baseline']['startIdx'] = PeakStartIdx
                    Presult[index]['Baseline']['endIdx'] = PeakEndIdx
        else:
            # 当存在融合峰
            PeakStartIdx = Presult[index]['Peaks'][0]['StartIdx']
            PeakEndIdx = Presult[index]['Peaks'][-1]['EndIdx']
            # 当第一峰的起点和基线的起点一致 并且 最后一个峰的终点和基线的终点一致时
            if PeakStartIdx == Presult[index]['Baseline']['startIdx'] and PeakEndIdx==Presult[index]['Baseline']['endIdx']:
                pass
            else:
                pass

def CalPercentArea(res):
    """
过滤后的也是需要重新计算百分面积的
    :param res: 过滤后的事件结果
    :return:返回重新计算后的百分面积
    """
    sumarea = 0
    sumheight = 0
    for cluster in res:
        for peak in cluster["Peaks"]:
            sumarea += peak["Area"]
            sumheight += peak["Height"]

    for cluster in res:
        for peak in cluster["Peaks"]:
            if sumarea ==0:
                peak["PercentArea"] ==0
            else:
                peak["PercentArea"] = peak["Area"]/sumarea*100
            if sumheight ==0:
                peak["PercentHeight"] ==0
            else:
                peak["PercentHeight"] = peak["Height"] / sumheight * 100
    return res

def Validata(res):
    PercentsumArea = 0
    Percentsumheight = 0
    for cluster in res:
        for peak in cluster["Peaks"]:
            PercentsumArea += peak["PercentArea"]
            Percentsumheight += peak["PercentHeight"]
    print(PercentsumArea)
    print(Percentsumheight)

# 手动积分
def ManualPeakDetection(rawData,timeInterval,timeParam):
    """
    # 手动积分算法
    :param rawData: 原始数据
    :param timeInterval: 采样频率
    :param detectionMethod: 手动积分的参数设法
    :return: status: int(运行状态代码，<0时为错误)
    :return: results: string(Json字符串形式的峰结果)
    """
    if rawData is None or len(rawData) <= 1:
        return -11  # 输入原始采样数据有误
    if timeInterval is None or timeInterval < 1:
        return -12, None  # 输入时间间隔错误
    inteval = timeInterval / 1000.0
    timeData = np.linspace(0, len(rawData) * inteval, num=len(rawData), endpoint=False)
    timeEvent = timeParam['params']
    err = PC.checkManualParam(timeEvent, timeData)
    if err < 0:
        return err,None
    # err, time = PC.loadMaunalParams(timeParam)
    Preresult = Ef.ManualIntegral(rawData, inteval, timeEvent)
    if Preresult != None:
        res = ET.calculateResult(Preresult, rawData, timeData, inteval)
        if len(res) == 0:
            return -60, None
    else:
        return -80,None #当前时间段不能进进行手动积分
    return 0, res




def runPeakDetection(rawData, timeInterval, detectionMethod):
    """
    根据指定的积分方法进行峰检测与积分操作
    :param: rawData: float[](所有采样数据点)
    :param: timeData: float[](所有采样时间)
    :param: detectionMethod: string(Json字符串形式的输入积分方法参数)
    :return: status: int(运行状态代码，<0时为错误)
    :return: results: string(Json字符串形式的峰结果)
    """

    if rawData is None or len(rawData) <=1: 
        return -11 # 输入原始采样数据有误
    if timeInterval is None or timeInterval < 1:
        return -12, None # 输入时间间隔错误
    inteval = timeInterval / 1000.0
    timeData = np.linspace(0, len(rawData)*inteval, num=len(rawData), endpoint=False)
    err, algo, pieces, events = PC.loadJsonParams(detectionMethod)
    if err < 0:
        return err, None # 输入方法json无法解析
    err = PC.checkParam(algo, pieces, events, timeData)
    if err < 0:
        return err, None # 输入方法json数据错误
    err, res1= Ef.runPiecewiseIntegral(rawData, timeData, inteval, algo, pieces)
    if err < 0:
        return err ,None
    err,Preresult,runResFilterEvents = Ef.SortEvents(events, algo, res1, inteval, rawData,timeData) #处理事件


    if err < 0 or Preresult ==None:
        return err, None # 积分算法内部出现错误
    try:
        res = ET.calculateResult(Preresult, rawData, timeData, inteval)
    except:
        return -30, None # 积分算法内部出现错误
    #当用户在没有输入积分参数的情况下，res结果是为None的 所以即使输入了适应性参数，也不会给计算，因为逻辑就是错误的
    if res ==None:
        return -30,None # 积分算法内部出现错误
    res = Ef.runResFilter(res, runResFilterEvents,inteval)  # 峰过滤事件

    # 如果没有过滤事件 则结果就是原本的结果 如果有则重新计算峰的百分高度和百分面积
    if len(runResFilterEvents) ==0:
        res = res
    else:
        res = CalPercentArea(res) #
    # Validata(res)  #验证百分面积和百分高度是否为100

    if len(res) ==0:
        return -60,None
    return 0, res

'''
def runApex(rawData, timeData, params):
    """
    输入原始数据、原始采样频率、积分参数，进行Apex积分(目前正在施工中，还不能进行峰事件的设定)
    总体设计按照：先进行分段事件峰检测，向后合并结果，再进行分段峰修改事件。
    :param: rawData: float[](所有采样数据点)
    :param: timeData: float[](采样时间轴)
    :param: allparams: dict(已经解析的积分方法参数)
    :return: status: int(运行状态代码，<0时为错误)
    :return: results: string(Json字符串形式的峰结果)
    """
    timeInteval = Series(timeData).diff()[1:].mean()
    try:
        StartTime = params["startTime"]
        # EndTime = params["endTime"]
        peakWidth = params["peakWidth"]
        apexThreshold = params["apexThreshold"]
        startThreshold = params["startThreshold"]
        endThreshold = params["endThreshold"]
        # minArea = params["minArea"]
        # minHeight = params["minHeight"]
        # events = params["events"]
    except:
        return -14, None # 缺少一个或多个必要积分参数输入
    try:
        res=Ai.main(raw_data = rawData, 
                    raw_step = timeInteval, 
                    Integral_time= StartTime,
                    apexthre= apexThreshold,
                    liftoff=startThreshold,
                    touchdown=endThreshold,
                    width=peakWidth)
        return 0, resultParsing(res,rawData)
    except:
        return -30, None # 积分算法内部出现错误

def runTrad(rawData, timeData, params):
    """
    输入原始数据、原始采样频率、积分参数，进行Apex积分(目前正在施工中，还不能进行峰事件的设定)
    总体设计按照：先进行分段事件峰检测，向后合并结果，再进行分段峰修改事件。
    
    :param: rawData: float[](所有采样数据点)
    :param: timeData: float[](采样时间轴)
    :param: allparams: dict(已经解析的积分方法参数)
    :return: status: int(运行状态代码，<0时为错误)
    :return: results: string(Json字符串形式的峰结果)
    """
    try:
        peakWidth = params["peakWidth"]
        Threshold = params["Threshold"]
        # minArea = params["minArea"]
        # minHeight = params["minHeight"]
        # events = params["events"]
    except:
        return -14, None # 缺少一个或多个必要积分参数输入
    try:
        res=Ti.traditonalIntegral(rawData=rawData, timeData=timeData, threshold=Threshold,
                                    detectWidth=peakWidth)
        # out = []
        # for item in res:
        #     out.append(json.loads(json.dumps(item, default=lambda o: o.__dict__, sort_keys=True, indent=4)))
        # out = json.dumps(out)
        return 0, res
    except:
        return -30, None # 积分算法内部出现错误
'''