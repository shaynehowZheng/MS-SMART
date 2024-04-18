from PICore.DetectMethod import Trad_Integral as Ti
from PICore.DetectMethod import Apex_Integral as Ai
from PICore.Events import TradEvents as TE
from PICore.Events import EventTool as ET
from PICore.Events import CommonEvents as Ce
from PICore import ParamCheck as PC


def runPiecewiseIntegral(rawData, timeData,inteval, algo, integralPieces):
    preRes = []
    pieceIndx = 0
    for piece in integralPieces:
        pieceIndx   =pieceIndx+1
        if piece["endTime"] ==-1:
            piece["endTime"] = timeData[-1]
        else:
            piece["endTime"] =piece["endTime"]
        if algo == "Traditional":
            err,curr = Ti.getTradPreResult(rawData, timeData, inteval, piece["startTime"], piece["endTime"], piece["tgtThreshold"], piece["peakWidth"])
            # 当存在分段积分时，如果第一个参数导致没有结果，则不 影响后面的分段积分的检测和结果
            if len(integralPieces) == 1:
                if err < 0:
                    return err, None
                else:
                    preRes += curr
            elif pieceIndx <= len(integralPieces):
                if err <0:
                    continue
                else:
                    preRes += curr
        elif algo == "Apex":
            err,curr = Ai.getApexPreResult(rawData, timeData,inteval, piece["startTime"], piece["endTime"],piece["peakWidth"], piece["apexThreshold"],piece["startThreshold"], piece["endThreshold"])
            if len(integralPieces) == 1:
                if err < 0:
                    return err, None
                else:
                    preRes += curr
            elif pieceIndx <= len(integralPieces):
                if err < 0:
                    continue
                else:
                    preRes += curr
        else:
            return -30, None
    if len(preRes)==0:
        return -60, None
    return 0, preRes

def ManualIntegral(rawData,timeInterval,TimeEvent):
    """
    手动积分 直接选择
    :param rawData:
    :param timeInterval:
    :param timeParam:
    :return:
    """
    ForcePeakRes = []
    for signal in TimeEvent:
        startTime = int(signal['startTime'] / timeInterval)
        endTime = int(signal['endTime'] / timeInterval)
        temp= TE.ForcePeak(startTime, endTime, rawData)
        ForcePeakRes.append(temp[0])
    # 当所选的强迫峰形成了一个完整的负峰 则返回一个NONE 即表示不允许用户这么操作
    return ForcePeakRes




def runResFilter(pre_result, filterPieces,inteval):
    '''
      峰过滤事件，可分为相同时间段的峰过滤和不同时间段的峰过滤，相同时间段的过滤可看作是and的关系
      不同时间段可看作是or的关系
      '''
    filterPieces = sorted(filterPieces, key=lambda x: x["endTime"])
    res = pre_result

    for piece in filterPieces:

        if piece["eventType"] == "SetMinArea":
            res = Ce.Areafilter(res, piece["startTime"], piece["endTime"], piece["paramValue"], inteval)

        if piece["eventType"] == "SetMinHeight" or piece["eventType"] == "SetMaxHeight":
            res = Ce.Heightfilter(res, piece["startTime"], piece["endTime"], piece["eventType"], piece["paramValue"],
                                  inteval)

        if piece["eventType"] == "SetMaxWidth":
            res = Ce.Widthfilter(res, piece["startTime"], piece["endTime"], piece["paramValue"], inteval)

    return res

def SortEvents(events,algo,Preresult,timeInterval,rawData,timeData):
    # 首先将事件进行分类
    err,BslEvents,PeakEvents,DroplineEvents,runResFilterEvents = PC.EventClassify(events,timeData)
    if err <0:
        return err,None,None #这个返回和真实的值返回配套
    Preresult = BslEvent(BslEvents, algo, Preresult, timeInterval, rawData)  # 修改基线类事件
    Preresult = sorted(Preresult, key=lambda x: x["bsl"][0])  # 按照基线的起始点排序
    Preresult = DroplineEvent(DroplineEvents, algo, Preresult, timeInterval, rawData)  # 强迫垂线事件
    Preresult = PeakEvent(PeakEvents, algo, Preresult, timeInterval, rawData) #强迫峰事件


    # if len(BslEvents) >=1:
    #     Preresult = BslEvent(BslEvents, algo, Preresult, timeInterval, rawData)  # 修改基线类事件
    # elif len(DroplineEvents) >=1:
    #     Preresult = DroplineEvent(DroplineEvents, algo, Preresult, timeInterval, rawData)  # 强迫垂线事件
    # elif len(PeakEvents)>=1:
    #     Preresult = PeakEvent(PeakEvents, algo, Preresult, timeInterval, rawData)  # 强迫峰事件
    # elif len(runResFilterEvents) >=1:
    #     Preresult = runResFilter(Preresult, runResFilterEvents, timeInterval)  # 峰过滤事件

    return 0,Preresult,runResFilterEvents


def BslEvent(events,algo,Preresult,timeInterval,rawData):
    """
    对输入积分事件进行分段处理，目前还需要设计一下，主要要分清楚需要重新跑积分的事件和直接在原有峰结果上修改的事件，以及他们之间的优先级关系，谁覆盖谁。
    修改积分参数类：后出现的事件覆盖前出现的，后修改的参数覆盖前修改的
    修改峰事件：修改峰的事件在峰检测完成后执行，后出现的事件覆盖前出现的
    返回的是分段结果
    """
    # 修改基线 修改基线中的子事件理论上是不能有交叉的 这段可能改掉
    for single in events:
        startTime = int(single['startTime'] / timeInterval)
        endTime = int(single['endTime'] / timeInterval)

        if algo == "Traditional":
            eventType = single["eventType"]
            # 按时间正向水平修改基线
            if eventType == "FwdHorizonByTime":
                ByTimeselectRes,Preresult = ET.selectResult(startTime, endTime, Preresult) #筛选出在事件时间中检测出的峰
                #当事件段内没有峰的时候
                if len(ByTimeselectRes) == 0:
                    continue
                else:
                    # 首先判断当前事件的时间的起终点是不是在一个峰的同一侧
                    Preresult, selectRes = ET.getStartandEnd(startTime, endTime, ByTimeselectRes, Preresult)
                    if len(selectRes) >=1:
                        FHBTPreresult,PreRes = TE.FwdHorizonByTime(startTime,endTime,ByTimeselectRes,Preresult) # 事件中峰的起终点的处理
                        if PreRes == None: #当所选的时间段内没有合适的结果
                            continue
                        Eventresult = ET.solveSingleEventbsl(PreRes, rawData) #对事件的结果进行检查处理 判断是否有基线穿透的问题
                        FHBTResult = ET.solveBaseline(FHBTPreresult, rawData)  # 同时也对事件之外的结果进行处理 判断是否有基线穿透的问题
                        # 当事件的结果为None 值时表示根据用户的需要 得不到相应合理的值 则返回None
                        if Eventresult == None:
                            Preresult = FHBTResult
                        else:
                            # Preresult = ET.solveFusedPeakEvent(Eventresult, FHBTResult)  # 对事件结果进行基线判断之后 可能存在融合峰变单峰的问题 此时应将融合峰断开
                            Preresult = ET.DealHorizontalBSL(rawData,Eventresult, FHBTResult) #对水平基线类事件做基线的特色处理，可能基线不在原始谱图对应的结果下面
                    else:
                        Preresult = ET.solveBaseline(Preresult, rawData) #将所选的区域弄成一个单独的峰 同时对事件外的峰进行判断 看是否存在基线穿透的问题


            # 按时间反向水平修改基线
            elif eventType == "RevrsHorizonByTime":
                RevByTimeselectRes, Preresult = ET.selectResult(startTime, endTime, Preresult)
                if len(RevByTimeselectRes) ==0:
                    continue
                else:
                    # 首先判断当前事件的时间的起终点是不是在一个峰的同一侧
                    Preresult, selectRes = ET.getStartandEnd(startTime, endTime, RevByTimeselectRes, Preresult)
                    if len(selectRes) >= 1:
                        RHBTPreresult,PreRes = TE.RevrsHorizonByTime(startTime,endTime,RevByTimeselectRes,Preresult)
                        if PreRes == None:  # 当所选的时间段内没有合适的结果
                            continue
                        Eventresult = ET.solveSingleEventbsl(PreRes, rawData) #对事件的结果进行检查处理 判断是否有基线穿透的问题
                        RHBTResult = ET.solveBaseline(RHBTPreresult, rawData)  # 同时也对事件之外的结果进行处理 判断是否有基线穿透的问题
                        # 当事件的结果为None 值时表示根据用户的需要 得不到相应合理的值 则返回None
                        if Eventresult == None:
                            Preresult = RHBTResult
                        else:
                            # Preresult = ET.solveFusedPeakEvent(Eventresult, RHBTResult)   # 对事件结果进行基线判断之后 可能存在融合峰变单峰的问题 此时应将融合峰断开
                            Preresult = ET.DealHorizontalBSL(rawData, Eventresult, RHBTResult)  # 对水平基线类事件做基线的特色处理，可能基线不在原始谱图对应的结果下面
                    else:
                        Preresult = ET.solveBaseline(Preresult, rawData) #将所选的区域弄成一个单独的峰 同时对事件外的峰进行判断 看是否存在基线穿透的问题


            # 按峰正向水平修改基线
            elif eventType == "FwdHorizonByPeak":
                ByPeakselectRes,Preresult = ET.selectResult(startTime, endTime, Preresult)
                if len(ByPeakselectRes)==0:
                    continue
                else:
                    # 首先判断当前事件的时间的起终点是不是在一个峰的同一侧
                    Preresult, selectRes = ET.getStartandEnd(startTime, endTime, ByPeakselectRes, Preresult)
                    if len(selectRes) >= 1:
                        FHBPPreresult,PreRes = TE.FwdHorizonByPeak(startTime,endTime,ByPeakselectRes,Preresult)
                        if PreRes == None:  # 当所选的时间段内没有合适的结果
                            continue
                        Eventresult = ET.solveSingleEventbsl(PreRes, rawData)  # 对事件的结果进行检查处理 判断是否有基线穿透的问题
                        FHBPResult = ET.solveBaseline(FHBPPreresult, rawData)  # 同时也对事件之外的结果进行处理 判断是否有基线穿透的问题
                        # 当事件的结果为None 值时表示根据用户的需要 得不到相应合理的值 则返回None
                        if Eventresult == None:
                            Preresult = FHBTResult
                        else:
                            # Preresult = ET.solveFusedPeakEvent(Eventresult,FHBPResult)  # 对事件结果进行基线判断之后 可能存在融合峰变单峰的问题 此时应将融合峰断开
                            Preresult = ET.DealHorizontalBSL(rawData, Eventresult, FHBPResult)
                    else:
                        Preresult = ET.solveBaseline(Preresult, rawData) #将所选的区域弄成一个单独的峰 同时对事件外的峰进行判断 看是否存在基线穿透的问题


            # 按峰反向水平修改基线
            elif eventType == "RevrsHorizonByPeak":
                RevByPeakselectRes, Preresult = ET.selectResult(startTime, endTime, Preresult)
                if len(RevByPeakselectRes) ==0:
                    continue
                else:
                    # 首先判断当前事件的时间的起终点是不是在一个峰的同一侧
                    Preresult, selectRes = ET.getStartandEnd(startTime, endTime, RevByPeakselectRes, Preresult)
                    if len(selectRes) >= 1:
                        RHBPPreresult,PreRes = TE.RevrsHorizonByPeak(startTime,endTime,RevByPeakselectRes,Preresult)
                        if PreRes == None:  # 当所选的时间段内没有合适的结果
                            continue
                        Eventresult = ET.solveSingleEventbsl(PreRes, rawData)  # 对事件的结果进行检查处理 判断是否有基线穿透的问题
                        RHBPResult = ET.solveBaseline(RHBPPreresult, rawData)  # 同时也对事件之外的结果进行处理 判断是否有基线穿透的问题
                        # 用户输入的时间内没有合适的时间段
                        if Eventresult == None:
                            Preresult = RHBPResult
                        else:
                            # Preresult = ET.solveFusedPeakEvent(Eventresult, RHBPResult)  # 对事件结果进行基线判断之后 可能存在融合峰变单峰的问题 此时应将融合峰断开
                            Preresult = ET.DealHorizontalBSL(rawData, Eventresult, RHBPResult)
                    else:
                        Preresult = ET.solveBaseline(Preresult, rawData) #将所选的区域弄成一个单独的峰 同时对事件外的峰进行判断 看是否存在基线穿透的问题

            # 按时间强迫基线
            elif eventType == "ForceBslByTime":
                ForceBslByTimeselectRes, Preresult = ET.selectResult(startTime, endTime, Preresult)
                if len(ForceBslByTimeselectRes) ==0:
                    continue
                else:
                    # 首先判断当前事件的时间的起终点是不是在一个峰的同一侧
                    Preresult, selectRes = ET.getStartandEnd(startTime, endTime, ForceBslByTimeselectRes, Preresult)
                    if len(selectRes) >= 1:
                        FBBTPreresult,PreSelect = TE.ForceBslByTime(startTime,endTime,ForceBslByTimeselectRes,Preresult)
                        #用户输入的时间点在结果中没有合适的时间段结果
                        if PreSelect == None:
                            Preresult = FBBTPreresult
                        else:
                            Preresult = ET.solveBaseline((FBBTPreresult+PreSelect), rawData)  # 对结果检测是否存在基线穿越的问问题
                    else:
                        Preresult = ET.solveBaseline(Preresult, rawData) #将所选的区域弄成一个单独的峰 同时对事件外的峰进行判断 看是否存在基线穿透的问题

            #   按峰强迫基线
            elif eventType == "ForceBslByPeak":
                ForceBslByPeakselectRes, Preresult = ET.selectResult(startTime, endTime, Preresult)
                if len(ForceBslByPeakselectRes) ==0:
                    continue
                else:
                    Preresult, selectRes = ET.getStartandEnd(startTime, endTime, ForceBslByPeakselectRes, Preresult)
                    if len(selectRes) >= 1:
                        FBBPPreresult,PreSelect= TE.ForceBslByPeak(startTime,endTime,ForceBslByPeakselectRes,Preresult)
                        if PreSelect == None:
                            Preresult = FBBPPreresult
                        else:
                           Preresult = ET.solveBaseline((FBBPPreresult +PreSelect), rawData) #对结果检测是否存在基线穿越的问问题
                    else:
                        Preresult = ET.solveBaseline(Preresult, rawData) #将所选的区域弄成一个单独的峰 同时对事件外的峰进行判断 看是否存在基线穿透的问题

            # 谷到谷事件
            elif eventType == "VtoV":
                Preresult = Ce.VtoV(Preresult, startTime, endTime, rawData) #

    return Preresult

#强迫垂线事件
def DroplineEvent(events,algo,Preresult,timeInterval,rawData):
    for single in events:
        startTime = int(single['startTime'] / timeInterval)
        endTime = int(single['endTime'] / timeInterval)
        if algo == "Traditional":
            eventType = single["eventType"]
            if eventType == "ForceDropline":
               Preresult =TE.ForceDropline(Preresult, startTime, rawData)
    return Preresult


# 强迫峰事件
def PeakEvent(events,algo,Preresult,timeInterval,rawData):
    for single in events:
        startTime = int(single['startTime'] /timeInterval)
        endTime = int(single['endTime'] / timeInterval)
        if algo == "Traditional":
            eventType = single["eventType"]
            if eventType == "ForcePeak":
                ForcePeakselectRes, Preresult = ET.selectResult(startTime, endTime, Preresult)
                ForcePeakRes = TE.ForcePeak(startTime, endTime, rawData)
                # 当所选的强迫峰形成了一个完整的负峰 则返回一个NONE 即表示不允许用户这么操作
                if ForcePeakRes ==None:
                    Preresult = Preresult
                else:
                    Preresult = ET.solveBaseline(Preresult, rawData)
                    Preresult = Preresult+ForcePeakRes

    return Preresult



