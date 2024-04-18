# -*- coding:utf-8 -*-
#@Time : 2020/11/9 11:51
#@Author: zgz
#@File : ParamCheck.py

import simplejson

# 事件的一些参数
__TRAD_PARAMS = ['startTime', 'endTime', 'peakWidth', 'tgtThreshold']
__APEX_PARAMS = ['startTime', 'endTime', 'peakWidth', 'apexThreshold','startThreshold', 'endThreshold']
__EVENT_PARAMS = ['startTime', 'endTime', 'eventType','paramValue']
__EVENT_COMNAMES =['AllowNegative','NoIntg','SetMinArea','SetMinHeight','SetMaxWidth','SetMaxHeight','TgtSkim']  #共有事件的名称
__EVENT_TRAD_NAMES_BSL =['ForceBslByPeak','ForceBslByTime','FwdHorizonByTime','RevrsHorizonByTime','FwdHorizonByPeak','RevrsHorizonByPeak','VtoV'] #传统积分特有事件的名称
__EVENT_APEX_NAMES =['DetectShoulder','GausSkim',"VtoV"]  #Apex积分特有事件的名称
__EVENT_TRAD_NAMES_PEAK = ['ForcePeak']
__EVENT_TRAD_NAMES_Dropline = ['ForceDropline']

def loadJsonParams(jsonparams):
    try:
        readparams = simplejson.loads(jsonparams)
        algo = readparams["algorithm"]
        pieces = readparams["params"]
        events = readparams["events"]
        return 0, algo, pieces, events
    except:
        return -13, None, None, None # 输入积分方法json无法解析或无相应字段


def loadMaunalParams(jsonparams):
    try:
        readparams = simplejson.loads(jsonparams)
        timeParam = readparams["params"]
        return 0, timeParam
    except:
        return -13,None

def checkManualParam(timeEvent,timeData):
    for item in timeEvent:
        if item["startTime"] is not None and item["endTime"] is not None:
            continue
        else:
            return -31  # 有时间段参数缺失

    preend = float('-inf')
    for item in sorted(timeEvent,key=lambda x:x["startTime"]):
        # 检查是否出现时间段交叉
        if item["startTime"] > preend and item["startTime"] < item["endTime"] and item["startTime"] >= timeData[0] and item["endTime"]<timeData[-1] :
            preend = item["endTime"]
        else:
            return -33  # 输入时间段有交叉或错误
    return 0



def checkParam(algo, pieces, events, timeData):
############# 1、检验积分参数！！！！！！！
    if algo == "Traditional":
        parnames = __TRAD_PARAMS

    elif algo == "Apex":
        parnames = __APEX_PARAMS
        eventNames = __EVENT_APEX_NAMES

    else:
        return -34 #输入的积分方法无效
    for item in pieces:
        # 检查参数是否与算法一致、或时间参数缺损
        if set(list(item.keys())) == set(parnames):
            if item["startTime"] is not None and item["endTime"] is not None:
                continue
            else:
                return -31  # 有时间段参数缺失
        else:
            return -32  # 参数与算法不一致

    preend = float('-inf')
    for item in sorted(pieces, key=lambda x: x["endTime"]):
        # print(item)
        # 检查是否出现时间段交叉
        if item["endTime"] == -1:
            continue
        elif item["startTime"] > preend and item["startTime"] < item["endTime"] and item["startTime"] >= timeData[0]:
            preend = item["endTime"]
        else:
            return -33  # 输入时间段有交叉或错误



############ 2、检验事件参数!!!!!!!!!
   #检验事件参数类型是否在事件范围

    if algo == "Traditional":
        eventNames = __EVENT_TRAD_NAMES_BSL+__EVENT_COMNAMES+__EVENT_TRAD_NAMES_PEAK+__EVENT_TRAD_NAMES_Dropline
    elif algo == "Apex":
        eventNames = __EVENT_APEX_NAMES+__EVENT_COMNAMES

    for item_name in events:
        # print(item_name)
        if item_name['eventType'] in eventNames:
            continue
        else:
            return -35  # 输入的事件参数有问题

    # 检验事件的时间参数是否存在
    eventsparnames = __EVENT_PARAMS
    for item in events:
        if set(list(item.keys())) == set(eventsparnames):
            if item["startTime"] is not None and item["endTime"] is not None:
                continue
            else:
                return -31  # 有时间段参数缺失
        else:
            return -32  # 参数与算法不一致

    #判断事件的时间参数是否存在异常 如终止时间小于开始时间
    for item_event in events:
        if item_event["startTime"] < item_event["endTime"]:
            continue
        else:
            return -36 #时间参数存在 异常

    return 0

def SameEventParamCheck(events,eventName,timeData):
    # 将事件分类 判断修改同一事件的时间是否交叉
    # 同一事件类型的时间 也判断是否存在交叉
    sameEventList = []
    preend = float('-inf')
    for item in events:
        if item['eventType'] in eventName:
            sameEventList.append(item)
    for item in sorted(sameEventList, key=lambda x: x["endTime"]):
        # 检查是否出现时间段交叉
        if item["startTime"] > preend and item["startTime"] < item["endTime"] and item["startTime"] >= timeData[0]:
            preend = item["endTime"]
        else:
            return -33  # 输入时间段有交叉或错误
    return 0


def EventClassify(events,timeData):
    BslEvents = []
    DroplineEvent = []
    PeakEvent = []
    runResFilterEvent = []
    for item in events:
        # 将基线类的事件分类并判断是否存在基线交叉 且也能判断同一事件的事件时间是否存在交叉
        #将基线归类
        if item['eventType'] in __EVENT_TRAD_NAMES_BSL:
            BslEvents.append(item)
        # 对事件进行事件参数的检验

        # 强迫峰事件
        elif item['eventType'] in __EVENT_TRAD_NAMES_PEAK:
            PeakEvent.append(item)

        # 强迫基线
        elif item['eventType'] in __EVENT_TRAD_NAMES_Dropline:
            DroplineEvent.append(item)

        #峰过滤事件  不做时间交叉检验
        elif item['eventType'] in __EVENT_COMNAMES:
            runResFilterEvent.append(item)
        # 对同一事件进行参数检验
    err = SameEventParamCheck(DroplineEvent, __EVENT_TRAD_NAMES_Dropline, timeData)
    if err < 0:
        return err, None,None,None,None
    err = SameEventParamCheck(PeakEvent, __EVENT_TRAD_NAMES_PEAK, timeData)
    if err < 0:
        return err, None,None,None,None
    err = SameEventParamCheck(BslEvents, __EVENT_TRAD_NAMES_BSL, timeData)
    if err < 0:
        return err, None,None,None,None


    return 0,BslEvents,PeakEvent,DroplineEvent,runResFilterEvent