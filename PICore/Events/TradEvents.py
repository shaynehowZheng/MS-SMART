
from PICore.Events import EventTool as ET


### 修改基线类
# 按时间正向水平修改基线
def FwdHorizonByTime(startTime, endTime, selectRes,Preresult):
    """
    :param startTime: 事件的开始时间
    :param endTime: 事件的结束时间
    :param finalRes: 事件时间段中的积分结果
    :return:
    """
# 基线的开始时间分为三种情况：1、当前时间段的起始时间在峰起点之前，则基线起点和峰起点一致,
    if startTime <= selectRes[0]['peak'][0][0]:
        bslStart = selectRes[0]['peak'][0][0]
    #2、当前时间段的起始时间在一个峰的保留时间前半部分
    elif startTime > selectRes[0]['peak'][0][0] and startTime < selectRes[0]['peak'][0][1]:
        bslStart = startTime
    else:
        Preresult, selectRes, bslStart = getbslStart(startTime,selectRes,Preresult)

# 基线的终止时间大体有两种情况：1、当前时间段内只包含一个峰的时候 2、当前时间段包含多个峰的时候且多个峰为融合峰的时候,
    # 判断最后一个融合峰与基线终点的关系。
    # 基线的终止时间分为三种情况：1、当前时间段的终止时间在最后一个峰的终点之后，则基线终点和最后一个峰的终止时间
    if endTime >= selectRes[-1]['peak'][-1][-1]:
        bslEnd = selectRes[-1]['peak'][-1][-1]
    # 2、当前时间段的终止时间在一个峰的保留时间后半部分
    elif  endTime > selectRes[-1]['peak'][-1][1] and endTime < selectRes[-1]['peak'][-1][2]:
        bslEnd = selectRes[-1]['peak'][-1][1]
    else:
        Preresult, selectRes, bslEnd = getbslEnd(endTime,selectRes,Preresult)
    #修改基线的起终点
    #把得到的结果处理为融合峰 并且修改
    if len(selectRes)!=0:
        selectRes[0]['peak'][0][0] = bslStart
        selectRes[-1]['peak'][-1][2] = bslEnd
        tempRes = []
        for i in range(len(selectRes)):
            tempRes.append(selectRes[i]["peak"])
        finalRes = {'peak':tempRes,'bsl':[bslStart,bslEnd],'bslType':'Horizon'}
    else:
        finalRes = None
    return Preresult,finalRes





# 按时间反向水平修改基线 与正向的区别在于以时间的终点为准来画基线 所以判断时与正向相反
def RevrsHorizonByTime(startTime, endTime,selectRes,Preresult):

    # 基线的终止时间分为三种情况：1、当前时间段的终止时间在最后一个峰的终点之后，则基线终点和最后一个峰的终止时间
    if endTime >= selectRes[-1]['peak'][-1][-1]:
        bslEnd = selectRes[-1]['peak'][-1][-1]
    # 2、当前事件的时间段的终止时间在最后一个峰的后半部分
    elif endTime > selectRes[-1]['peak'][-1][1] and endTime < selectRes[-1]['peak'][-1][2]:
        bslEnd = endTime
    else:
        Preresult, selectRes, bslEnd = getbslEnd(endTime,selectRes,Preresult)


    if len(selectRes):
        # 基线的开始时间分为三种情况：1、当前时间段的起始时间在峰起点之前，则基线起点和峰起点一致,
        if startTime <= selectRes[0]['peak'][0][0]:
            bslStart = selectRes[0]['peak'][0][0]
        # 2、当前时间段的起始时间在一个峰的保留时间前半部分
        elif startTime > selectRes[0]['peak'][0][0] and startTime < selectRes[0]['peak'][0][1]:
            bslStart = selectRes[0]['peak'][0][0]
        else:
            # 3、在一个峰的保留时间后半部分 第一个是单峰的时候
            Preresult, selectRes, bslStart = getbslStart(startTime, selectRes, Preresult)

    if len(selectRes)!=0:
        selectRes[0]['peak'][0][0] = bslStart
        selectRes[-1]['peak'][-1][2] = bslEnd
        tempRes = []
        for i in range(len(selectRes)):
            tempRes.append(selectRes[i]["peak"])
        finalRes = {'peak':tempRes,'bsl':[bslStart,bslEnd],'bslType':'RevHorizon'}
    else:
        finalRes = None
    return Preresult,finalRes

# 按峰正向水平修改基线
def FwdHorizonByPeak(startTime, endTime,selectRes,Preresult):
    # 基线的开始时间分为三种情况：1、当前时间段的起始时间不在峰起点之前，则基线起点和峰起点一致,
    if startTime <= selectRes[0]['peak'][0][0]:
        bslStart = selectRes[0]['peak'][0][0]
    # 2、当前时间段的起始时间在一个峰的保留时间前半部分
    elif startTime > selectRes[0]['peak'][0][0] and startTime < selectRes[0]['peak'][0][1]:
        bslStart = selectRes[0]['peak'][0][0]
    else:
        # 3、在一个峰的保留时间后半部分 第一个是单峰的时候
        Preresult, selectRes, bslStart = getbslStart(startTime,selectRes,Preresult)

    # 基线的终止时间大体有两种情况：1、当前时间段内只包含一个峰的时候 2、当前时间段包含多个峰的时候且多个峰为融合峰的时候,
    # 判断最后一个融合峰与基线终点的关系。有两种情况 在峰的前半部分 和在峰的后半部分 且在单峰和融合峰时处理的情况不一样样
    if endTime >= selectRes[-1]['peak'][-1][-1]:
        bslEnd = selectRes[-1]['peak'][-1][-1]
    # 2、当前时间段的起始时间在一个峰的保留时间后半部分
    elif  endTime > selectRes[-1]['peak'][-1][1] and endTime < selectRes[-1]['peak'][-1][2]:
        bslEnd = selectRes[-1]['peak'][-1][2]
    else:
        #  当最后的结果是单峰时
        Preresult, selectRes, bslEnd = getbslEnd(endTime,selectRes,Preresult)

 # 将结果进行整合 处理成一个融合峰方便后续的 判断基线是否穿越的问题
    if len(selectRes)!=0:
        selectRes[0]['peak'][0][0] = bslStart
        selectRes[-1]['peak'][-1][2] = bslEnd
        tempRes = []
        for i in range(len(selectRes)):
            tempRes.append(selectRes[i]["peak"])
        finalRes = {'peak': tempRes, 'bsl': [bslStart, bslEnd], 'bslType': 'Horizon'}
    else:
        finalRes = None
    return Preresult, finalRes

def RevrsHorizonByPeak(startTime, endTime,selectRes,Preresult):
    # 基线的开始时间分为三种情况：1、当前时间段的起始时间不在峰起点之前，则基线起点和峰起点一致,
    if startTime <= selectRes[0]['peak'][0][0]:
        bslStart = selectRes[0]['peak'][0][0]
    # 2、当前时间段的起始时间在一个峰的保留时间前半部分
    elif startTime > selectRes[0]['peak'][0][0] and startTime < selectRes[0]['peak'][0][1]:
        bslStart = selectRes[0]['peak'][0][0]
    else:
        # 3、在一个峰的保留时间后半部分 第一个是单峰的时候
        Preresult, selectRes, bslStart = getbslStart(startTime,selectRes,Preresult)

    # 基线的终止时间大体有两种情况：1、当前时间段内只包含一个峰的时候 2、当前时间段包含多个峰的时候且多个峰为融合峰的时候,
    # 判断最后一个融合峰与基线终点的关系。有两种情况 在峰的前半部分 和在峰的后半部分 且在单峰和融合峰时处理的情况不一样样
    if len(selectRes)!=0:
        if endTime >= selectRes[-1]['peak'][-1][-1]:
            bslEnd = selectRes[-1]['peak'][-1][-1]
        # 2、当前时间段的起始时间在一个峰的保留时间前半部分
        elif  endTime > selectRes[-1]['peak'][-1][1] and endTime < selectRes[-1]['peak'][-1][2]:
            bslEnd = selectRes[-1]['peak'][-1][2]
        else:
            #  当最后的结果是单峰时
            Preresult, selectRes, bslEnd = getbslEnd(endTime, selectRes, Preresult)

    # 将结果进行整合 处理成一个融合峰方便后续的 判断基线是否穿越的问题
    if len(selectRes)!=0:
        selectRes[0]['peak'][0][0] = bslStart
        selectRes[-1]['peak'][-1][2] = bslEnd
        tempRes = []
        for i in range(len(selectRes)):
            tempRes.append(selectRes[i]["peak"])
        finalRes = {'peak': tempRes, 'bsl': [bslStart, bslEnd], 'bslType': 'RevHorizon'}
    else:
        finalRes = None
    return Preresult, finalRes


def ForceBslByTime(startTime, endTime,selectRes,Preresult):
# 基线的开始时间分为三种情况：1、当前时间段的起始时间在峰起点之前，则基线起点和峰起点一致,
    if startTime <= selectRes[0]['peak'][0][0]:
        bslStart = selectRes[0]['peak'][0][0]
    #2、当前时间段的起始时间在一个峰的保留时间前半部分
    elif startTime > selectRes[0]['peak'][0][0] and startTime < selectRes[0]['peak'][0][1]:
        bslStart = startTime
    else:
       # 3、在一个峰的保留时间后半部分
       Preresult, selectRes, bslStart = getbslStart(startTime,selectRes,Preresult)
    if len(selectRes) !=0:
        # 基线的终止时间大体有两种情况：1、当前时间段内只包含一个峰的时候 2、当前时间段包含多个峰的时候且多个峰为融合峰的时候,
        # 判断最后一个融合峰与基线终点的关系。
        if endTime >= selectRes[-1]['peak'][-1][-1]:
            bslEnd = selectRes[-1]['peak'][-1][-1]
        # 2、当前时间段的起始时间在一个峰的保留时间前半部分
        elif  endTime > selectRes[-1]['peak'][-1][1] and endTime < selectRes[-1]['peak'][-1][2]:
            bslEnd = selectRes[-1]['peak'][-1][2]
        else:
            #  当最后的结果是单峰时
            Preresult, selectRes, bslEnd = getbslEnd(endTime, selectRes, Preresult)

    #修改基线的起终点
    if len(selectRes)!=0:
        selectRes[0]['bsl'][0] = bslStart
        selectRes[-1]['bsl'][-1] = bslEnd
        selectRes[0]['peak'][0][0] = bslStart
        selectRes[-1]['peak'][-1][2] = bslEnd
    else:
        selectRes =None
    return Preresult, selectRes


def ForceBslByPeak(startTime, endTime,selectRes,Preresult):
    #不是水平的 所以不融合成水平基线
    # 基线的开始时间分为三种情况：1、当前时间段的起始时间在峰起点之前，则基线起点和峰起点一致,
    if startTime <= selectRes[0]['peak'][0][0]:
        bslStart = selectRes[0]['peak'][0][0]
    # 2、当前时间段的起始时间在一个峰的保留时间前半部分
    elif startTime > selectRes[0]['peak'][0][0] and startTime < selectRes[0]['peak'][0][1]:
        bslStart = selectRes[0]['peak'][0][0]
    else:
        # 3、在一个峰的保留时间后半部分 第一个是单峰的时候
        Preresult, selectRes, bslStart = getbslStart(startTime,selectRes,Preresult)

    # 基线的终止时间大体有两种情况：1、当前时间段内只包含一个峰的时候 2、当前时间段包含多个峰的时候且多个峰为融合峰的时候,
    # 判断最后一个融合峰与基线终点的关系。
    if len(selectRes)!=0:
        if endTime >= selectRes[-1]['peak'][-1][-1]:
            bslEnd = selectRes[-1]['peak'][-1][-1]
        # 2、当前时间段的起始时间在一个峰的保留时间前半部分
        elif  endTime > selectRes[-1]['peak'][-1][1] and endTime < selectRes[-1]['peak'][-1][2]:
            bslEnd = selectRes[-1]['peak'][-1][2]
        else:
            #  当最后的结果是单峰时
            Preresult, selectRes, bslEnd = getbslEnd(endTime, selectRes, Preresult)
            # 当前时间段在积分时间内 则峰的终止时间在事件的终止时间之前
    #修改基线 当所在的时间段内还有事件结果时
    if len(selectRes)!=0:
        selectRes[0]['bsl'][0] = bslStart
        selectRes[-1]['bsl'][-1] = bslEnd
        selectRes[0]['peak'][0][0] = bslStart
        selectRes[-1]['peak'][-1][2] = bslEnd
    else:
        selectRes = None
    return Preresult, selectRes

# 强迫垂线
def ForceDropline(Preresult,DroplinePoint,rawdata):
    #强迫垂线：包括两种情况：1、如果在峰的前半部分 则寻找当前峰起点到垂线点的最高点 如果最高点事垂线点 则垂线点也是峰顶点 起点保持不变
    #2、在峰的后半部分亦是如此进行

    for i in range(len(Preresult)):
        #当垂线点在一个单峰中
        for j in range(len(Preresult[i]['peak'])):
            # 当垂线点在峰的前半部分
            if Preresult[i]['peak'][j][0] <= DroplinePoint and DroplinePoint <= Preresult[i]['peak'][j][1]:
                startPoint = Preresult[i]['peak'][j][0]
                tempData = rawdata[startPoint:DroplinePoint+1]
                peakIndex = tempData.argmax()
                #当最大值在边界的时候 如果末尾加1则会超出边界值了
                if peakIndex+1 ==len(tempData):
                    peakPiont = Preresult[i]['peak'][j][0] + peakIndex
                else:
                    peakPiont = Preresult[i]['peak'][j][0]+peakIndex+1
                if peakPiont ==DroplinePoint:
                    temppeak = [startPoint, peakPiont, peakPiont]
                else:
                    temppeak = [startPoint, peakPiont, DroplinePoint]
                Preresult[i]['peak'][j][0] = DroplinePoint
                Preresult[i]['peak'].insert(j,temppeak) #如果在前半部分 则在插入到峰中前一个中
                break
                # 当垂线带你在单峰的后半部分
            elif Preresult[i]['peak'][j][1] <= DroplinePoint and DroplinePoint <= Preresult[i]['peak'][j][2]:
                endPoint = Preresult[i]['peak'][j][2]
                tempData = rawdata[DroplinePoint:endPoint+1]
                peakIndex = tempData.argmax()
                if peakIndex+1 ==len(tempData):
                    peakPiont = DroplinePoint+ peakIndex
                else:
                    peakPiont = DroplinePoint+ peakIndex + 1

                if peakPiont == DroplinePoint:
                    temppeak = [peakPiont, peakPiont,endPoint]
                else:
                    temppeak = [DroplinePoint, peakPiont, endPoint]
                Preresult[i]['peak'][j][2] = DroplinePoint
                Preresult[i]['peak'].insert(j+1, temppeak)  # 如果在前后半部分 则在插入到单峰的第二个
                break
    return Preresult



def ForcePeak(startTime, endTime,rawData):
    #强迫峰的情况可分为两种：1、当所选的时间内有峰 则不管时间段内的峰怎么样 直接将事件时间的起终点作为基线的起终点 然后寻找这段时间内最大的值作为峰顶点 但用户在这种情况应该操作的实际意义不大
    # 在没有检测到的峰的地方拉一条基线，其中包括没有积分区域，则我只在这段区域找一个最大值作为峰顶点(若没有最大值)
    #当强迫峰事件段内有峰的时候
    #强迫峰的起点

    #当强迫峰事件段内没有峰的情况
    finalRes= []
    peakStart = startTime
    peakEnd = endTime
    tempData = rawData[peakStart:peakEnd+1]
    peakBsl = ET.generatebasline(rawData,peakStart,peakEnd)

    ############################## 不做判断 直接形成强迫峰##################
    index = tempData.argmax()
    peakPoint = peakStart+index+1
    tempRes = {'peak':[],'bsl':None,'bslType':'auto'}
    tempRes['peak'].append([peakStart,peakPoint,peakEnd])
    tempRes['bsl'] = [peakStart,peakEnd]
    finalRes.append(tempRes)

    ############################## 用于判断是否形成了负峰##################
    # all() 函数用于判断给定的可迭代参数 iterable 中的所有元素是否都为 TRUE，如果是返回 True，否则返回 False。
    # 元素除了是 0、空、None、False 外都算 True。
    # if tempData.all() < peakBsl.all():
    #     finalRes = None
    # else:
    #     index = tempData.argmax()
    #     peakPoint = peakStart+index+1
    #     tempRes = {'peak':[],'bsl':None,'bslType':'auto'}
    #     tempRes['peak'].append([peakStart,peakPoint,peakEnd])
    #     tempRes['bsl'] = [peakStart,peakEnd]
    #     finalRes.append(tempRes)
    return finalRes




#共用的方法
# 1、获取事件的起始点
def getbslStart(startTime,selectRes,Preresult):
    if len(selectRes)==1 and  len(selectRes[0]['peak']) == 1:
        #当结果只有一峰时，且时间在峰的后半部分
        if startTime >= selectRes[0]['peak'][0][1] and startTime < selectRes[0]['peak'][0][2]:
            bslStart = None
            tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
            tempgroup['peak'].append(selectRes[0]['peak'][0])
            tempbslcoef = [selectRes[0]['peak'][0][0], selectRes[0]['peak'][0][2]]
            tempgroup['bsl'] = tempbslcoef
            Preresult.append(tempgroup)
            # 此时应该是断开 而不是删除
            del (selectRes[0])  # 则把不符合的过滤掉
    # 3、 第一个是单峰的时候,且当前时间段的事件中只有一个峰的时候
    elif len(selectRes[0]['peak']) == 1 and startTime >= selectRes[0]['peak'][0][1] and startTime < selectRes[0]['peak'][0][2]:
        bslStart = selectRes[1]['peak'][0][0]
        tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
        tempgroup['peak'].append(selectRes[0]['peak'][0])
        tempbslcoef = [selectRes[0]['peak'][0][0], selectRes[0]['peak'][0][2]]
        tempgroup['bsl'] = tempbslcoef
        Preresult.append(tempgroup)
        # 此时应该是断开 而不是删除
        del (selectRes[0])  # 则把不符合的过滤掉
    else:
        # 第一个是融合峰的时候
        if startTime >= selectRes[0]['peak'][0][1] and startTime <= selectRes[0]['peak'][0][2]:
            bslStart = selectRes[0]['peak'][1][0]
            tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
            tempgroup['peak'].append(selectRes[0]['peak'][0])
            tempbslcoef = [selectRes[0]['peak'][0][0], selectRes[0]['peak'][0][2]]
            tempgroup['bsl'] = tempbslcoef
            Preresult.append(tempgroup)
            # 此时应该是断开 而不是删除
            del (selectRes[0]['peak'][0])  # 则把不符合的过滤掉
            selectRes[0]['bsl'][0] = selectRes[0]['peak'][0][0]  # 修改基线
    return Preresult, selectRes ,bslStart


#2、获取事件的终点
def getbslEnd(endTime,selectRes,Preresult):
    #  当最后的结果是单峰时
    # 最最后是单峰的时候，且结果中只有一个峰时，小于某个峰的保留时间 并大于这个峰的开始时间，则结果中没有满足的条件的结果峰
    if len(selectRes)==1 and len(selectRes[0]['peak']) == 1:
        if endTime >= selectRes[-1]['peak'][0][0] and endTime <= selectRes[-1]['peak'][0][1]:
            bslEnd = None
            tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
            tempgroup['peak'].append(selectRes[-1]['peak'][-1])
            tempbslcoef = [selectRes[-1]['peak'][-1][0], selectRes[-1]['peak'][-1][2]]
            tempgroup['bsl'] = tempbslcoef
            Preresult.append(tempgroup)
            del (selectRes[-1])
    elif len(selectRes[-1]['peak']) == 1:
        # 当最后是单峰的时候 小于某个峰的保留时间 并大于这个峰的开始时间 基线的终点等于前一个峰的终止时间(结果中有多个峰)
        if endTime >= selectRes[-1]['peak'][0][0] and endTime <= selectRes[-1]['peak'][0][1]:
            bslEnd = selectRes[-2]['peak'][-1][2]
            tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
            tempgroup['peak'].append(selectRes[-1]['peak'][-1])
            tempbslcoef = [selectRes[-1]['peak'][-1][0], selectRes[-1]['peak'][-1][2]]
            tempgroup['bsl'] = tempbslcoef
            Preresult.append(tempgroup)
            del (selectRes[-1])  # 删除最后一个不作为融合峰的一个
        # 当最后是单峰的时候 大于某个峰的保留时间 并小于这个峰的终止时间 基线的终点等于这个峰的终止时间
        elif endTime >= selectRes[-1]['peak'][0][1] and endTime <= selectRes[-1]['peak'][0][2]:
            bslEnd = selectRes[-1]['peak'][-1][2]
    # 当最后的结果是融合峰时
    else:
        # 当最后是融合峰的时候  大于某个峰的保留时间 并小于这个峰的终止时间 基线的终点这个峰的终止时间
        if endTime >= selectRes[-1]['peak'][-1][1] and endTime <= selectRes[-1]['peak'][-1][2]:
            bslEnd = selectRes[-1]['peak'][-1][2]
        # 当最后是融合峰的时候  小于某个峰的保留时间 并大于这个峰的开始时间 基线的终点前一个峰的终止时间
        elif endTime >= selectRes[-1]['peak'][-1][0] and endTime <= selectRes[-1]['peak'][-1][1]:
            bslEnd = selectRes[-1]['peak'][-2][2]
            tempgroup = {'peak': [], 'bsl': None, 'bslType': "auto"}
            tempgroup['peak'].append(selectRes[-1]['peak'][-1])
            tempbslcoef = [selectRes[-1]['peak'][-1][0], selectRes[-1]['peak'][-1][2]]
            tempgroup['bsl'] = tempbslcoef
            Preresult.append(tempgroup)
            del (selectRes[-1]['peak'][-1])  # 删除最后一个不作为融合峰的一个

    return Preresult, selectRes, bslEnd

