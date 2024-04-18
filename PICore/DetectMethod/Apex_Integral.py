import numpy as np
from scipy import signal
from PICore.common import peak

def runpreFilter(rawData):
    # 前置滤波器
    b, skiprow=signal.butter(8, 0.6, 'lowpass')  
    rawData=signal.filtfilt(b, skiprow, rawData)
    #rawData = signal.medfilt(rawData,5)
    return rawData

def autopeakWidth(rawData,topindex=None):
    """
    函数:自动峰宽计算
    :param rawData: array_like(需要计算峰宽的数组)
    :param topindex: int(需要计算峰宽的峰顶点的index,为'None'时计算自动峰宽)
    :return peakwidth: int(峰宽值取整)
    """
    if topindex is None:
        derivative_2 = -np.diff(rawData, 2)#二阶导
        derivative_2 = np.hstack((derivative_2[0:2],derivative_2))
        topindex = np.argmax(derivative_2)#寻找二阶导最大的索引值，即为最高峰顶点
    leftinf_index, rightinf_index, _ = inflectionDetection(rawData, topindex)#找拐点
    t_peakwidth = (rightinf_index - leftinf_index) * 4.89549 / 2. #峰宽值为拐点之间的时间差*4.89549 / 2
    peakwidth = int(np.round(t_peakwidth)) #对峰宽值取整
    return peakwidth

def filter_parameter(rawData):
    """
    函数：滤波系数计算(必须为奇数(为偶数则+1)且不小于7)
    :param rawData: array_like(需要计算滤波系数的数组)
    :return flt_par: int(计算得到的滤波系数)
    """
    auto_pw = autopeakWidth(rawData)
    flt_par = int(np.round(auto_pw / 1.5))
    if flt_par % 2 == 0:
        flt_par += 1
    if flt_par <= 7:
        flt_par = 7

    # print('filter parameter =', flt_par)
    return flt_par

def autoapexThreshold(rawData,dataFreq):
    """
    函数：自动阈值计算
    * 阈值为噪声区域标准差六倍
    :param rawData: array_like(需要计算阈值的数组)
    :return thre: float(通过输入噪声区域计算得到的阈值)
    """
    # thredata=np.array(rawData[20:136])
    pointNum = int(20 * dataFreq)
    thredata = np.array(rawData[:pointNum])
    # thredata = np.array(rawData[-pointNum:])
    thredata=signal.resample(thredata,106)
    if dataFreq < 5:
        thre = np.std(thredata) * 4.
    else:
        thre = np.std(thredata) * 6.
    return thre

def inflectionDetection(rawData, topindex, derivative_2=None):
    """
    函数：求拐点及拐点间隔
    * 仅考虑正峰情况
    :param rawData: array_like(输入需要求拐点的数组)
    :param topindex: int(Apex峰顶点的index)
    :param derivative_2: array_like, optional(二阶导数，缺省时自动计算)
    :return leftinf_index : int(左侧拐点index)
    :return rightinf_index : int(右侧拐点index)
    :return rightinf_index - leftinf_index: int(拐点距离)
    """
    n = len(rawData)
    if derivative_2 is None:
        #求二阶导
        derivative_2 = -np.diff(rawData, 2)
        derivative_2 = np.hstack((derivative_2[0:2],derivative_2))
    rightinf_index=topindex +1
    leftinf_index=topindex -1
    for i in range(topindex - 1, 1, -1):
        if i != 1:#在数据结尾处强制终止
            if derivative_2[i] == 0 or derivative_2[i] * derivative_2[min(i+1,n-2)] < 0:#寻找左侧拐点索引
                #加入判断
                if derivative_2[max(i - 1,1)] * derivative_2[min(i + 2,topindex - 1)] > 0:
                    continue
                leftinf_index = i
                break
            #if d_2[max(i - 2,1)] * d_2[min(i + 3,x - 1)] > 0:
                #continue
            #判断点二阶导是否为0
    for j in range(topindex, n - 2, 1):
        if j == n - 3:#在数据结尾处强制终止
            return leftinf_index , j , j - leftinf_index
        if derivative_2[j] == 0 or derivative_2[max(j - 1,1)] * derivative_2[j] < 0:#寻找右侧拐点索引
            if derivative_2[max(j - 1,topindex)] * derivative_2[min(j + 2,n - 2)] > 0:
                continue
            rightinf_index=j
            break
    #print(lp + 1, rp + 1, rp-lp)
    return leftinf_index , rightinf_index , rightinf_index - leftinf_index

def generatebasline(data,bs,be):
    """
    生成基线数据
    :param data: array_like(原数据)
    :param bs: int(基线起点)
    :param be: int(基线终点)
    :return bsl: array_like(基线上的点数组)
    """
    y1=data[bs]
    y2=data[be]
    num=be-bs+1
    bsl=np.linspace(y1,y2,num) #linspace相当于根据起止点生成这条直线上的数据
    return bsl

def slope(data,strpoint,endpoint,step=1):
    """
    函数：计算某曲线上给定两点连线的斜率以及漂移 
    :param data: array_like(原数据数组)
    :param strpoint: int(计算斜率的起点)
    :param endpoint: int(计算斜率的终点)
    :param step: float, optional(数据点间的时间间隔，默认值为1)
    :return slp: float(计算得到的线段斜率)
    """
    if strpoint == endpoint:
        return 0, 0

    elif strpoint < endpoint:
        x = [strpoint * step, endpoint * step]
        y = [data[strpoint], data[endpoint]]

    else:
        x = [endpoint * step, strpoint * step]
        y = [data[endpoint], data[strpoint]]
    slp = (y[1] - y[0]) / (x[1] - x[0])
    return slp

def findleftBaseline(rawData,lp,rp,pmin,str_thre):
    """
    函数：寻找基线左边的端点
    :param rawData: array_like(原数据数组)
    :param lp: int(开始时基线左边端点)
    :param rp: int(开始时基线右边端点)
    :param pmin:int(左边最小的索引值)
    :param str_thre: float(左边端点的阈值)
    :return leftpoint: int(满足阈值的基线左边端点)
    """
    derivative_1=np.diff(rawData,1)
    slp=slope(rawData,lp,rp) 
    leftpoint=lp
    for i in range(lp,pmin,-1):
        slp=slope(rawData,i,rp) #求基线斜率
        tem_lp=derivative_1[i]-slp #切线斜率与基线斜率的差值
        if tem_lp<=str_thre: #与阈值进行比较
            if derivative_1[i-1]-slope(rawData,i-1,rp)<=str_thre:
                if derivative_1[i-2]-slope(rawData,i-2,rp)<=str_thre:
                    if derivative_1[i-3]-slope(rawData,i-3,rp)<=str_thre:
                        if derivative_1[i-4]-slope(rawData,i-4,rp)<=str_thre:
                            if derivative_1[i-5]-slope(rawData,i-5,rp)<=str_thre:
                                #print('l',derivative_1[i-5],slope(rawData,i-5,rp))
                                leftpoint=i
                                break
        if i<=pmin+1:
            leftpoint=i
            break
    return leftpoint,tem_lp

def findrightBaseline(rawData,lp,rp,pmax,end_thre):
    """
    函数：寻找基线右边的端点
    :param rawData: array_like(原数据数组)
    :param lp: int(开始时基线左边端点)
    :param rp: int(开始时基线右边端点)
    :param pmax:int(右边最大的索引值)
    :param end_thre: float(右边端点的阈值)
    :return rightpoint: int(满足阈值的基线右边端点)
    """
    derivative_1=np.diff(rawData,1)
    slp=slope(rawData,lp,rp)
    rightpoint= rp
    for i in range(rp,pmax,1):
        slp=slope(rawData,lp,i) #求基线斜率
        tem_rp=-derivative_1[i]+slp #切线斜率与基线斜率的差值
        if tem_rp<=end_thre: #与阈值经行比较
            if -derivative_1[i+1]+slope(rawData,lp,i+1)<=end_thre:
                if -derivative_1[i+2]+slope(rawData,lp,i+2)<=end_thre:
                    if -derivative_1[i+3]+slope(rawData,lp,i+3)<=end_thre:
                        if -derivative_1[i+4]+slope(rawData,lp,i+4)<=end_thre:
                            if -derivative_1[i+5]+slope(rawData,lp,i+5)<=end_thre:
                                #print('r',derivative_1[i+5],slope(rawData,lp,i+5))
                                rightpoint=i
                                break
        if i>=pmax-1:
            rightpoint=i
            break
    return rightpoint,tem_rp

def runApexDetection(rawData,Integral_zone,filter_par2,apexThreshold):
    '''
    函数：峰顶点检测
    :param rawData: array_like(用于计算apex峰顶点的原数据)
    :param Integral_zone: list(积分区域)
    :param filter_par2: int(全局二阶过滤系数)
    :param apexThreshold: float(全局峰顶点检测阈值)
    :return: apex_peak : array_like(所有apex峰顶点的index)
    '''
    left_zone=Integral_zone[0]

    derivative_1 = np.diff(rawData,1)
    tem_derivative_2=-np.diff(derivative_1.T,1) 
    derivative_2=np.hstack((tem_derivative_2[0:2],tem_derivative_2)) #补全最前面缺少的两个二阶导数点
    filter_der2=signal.savgol_filter(derivative_2, filter_par2, 2) #滤波后的二阶导数

    arg_order = int(filter_par2 / 7)
    if arg_order < 7:
        arg_order = 7
    tem_peakTopIdx=signal.argrelmax(filter_der2,order=arg_order)

    #tem_peakTopIdx=signal.find_peaks(rawData,height=threshold(filter_der2)/4) 
    #find_peaks是对原始数据进行寻找极值，使用该函数的原因是对于平缓的极值，argrelmax可能会漏掉局部极值

    apex_index = tem_peakTopIdx[0]
    
    #apex_index =apex_index[(apex_index<np.asarray(right_zone))&(apex_index>np.asarray(left_zone))]
    
    apex_peak = []#峰顶点下标数组
    i=0

    while True:
        if filter_der2[apex_index[i]] > apexThreshold:
            if i==0:
                minl=max(apex_index[i]-10,0) + np.argmin(rawData[max(apex_index[i]-10,0):apex_index[i]])
                minr=apex_index[i] + np.argmin(rawData[apex_index[i]:apex_index[i+1]])
            elif i == len(apex_index) - 1:
                minl=apex_index[i-1] + np.argmin(rawData[apex_index[i-1]:apex_index[i]])
                minr=apex_index[i] + np.argmin(rawData[apex_index[i]:min(apex_index[i]+10,len(rawData)-1)])
            else:
                minl=apex_index[i-1] + np.argmin(rawData[apex_index[i-1]:apex_index[i]])
                minr=apex_index[i] + np.argmin(rawData[apex_index[i]:apex_index[i+1]])
            useb=minl + np.argmax(rawData[minl:minr])
            
            if useb - minl >= 2 and minr - useb >= 2:
                apex_peak.append(useb + left_zone)
            else:
                apex_index=np.delete(apex_index,i)
                i-=1
        i+=1
        if i>=len(apex_index)-1:
            break

    apex_peak=(list(set(apex_peak)))
    apex_peak.sort()

    return apex_peak

def solveBaseline(rawData,peakTopIdx,startThreshold=0.,endThreshold=0.005):
    """
    函数：用于第一次确定基线，由Apex点左右的拐点连线作为预备基线
    :param rawData: array_like(需要计算基线的原数据)
    :param peakTopIdx: array_like(顶点index数组)
    :param startThreshold: float, optional(峰起点阈值，默认值为0)
    :param endThreshold: float, ptional(峰终点阈值，默认值为0.05)
    :return peak : array_like(峰特征数组，包括峰顶点、基线起止点，eg:[[基线起点，顶点，基线终点]])
    """
    peak=[]
    derivative_1 = np.diff(rawData,1)
    tem_derivative_2=-np.diff(derivative_1.T,1)
    derivative_2=np.hstack((tem_derivative_2[0:2],tem_derivative_2)) #补全最前面缺少的两个二阶导数点

    #一次仅处理一个峰
    for i in range(len(peakTopIdx)):
    
        infl_l,infl_r,_=inflectionDetection(rawData,peakTopIdx[i],derivative_2) #求左右拐点
        infl_slope=slope(rawData,infl_l,infl_r)
        #将斜率转为角度后计算后再化为斜率
        str_thre=np.tan((np.arctan(derivative_1[infl_l])-np.arctan(infl_slope))*startThreshold)
        end_thre=np.tan((np.arctan(-derivative_1[infl_r])+np.arctan(infl_slope))*endThreshold)

        leftpoint=infl_l
        rightpoint=infl_r
        # 判断起点和终点是否一致 如果一致则排除掉该峰
        if rightpoint == peakTopIdx[i] or leftpoint==peakTopIdx[i]:
            continue
        else:
            while True:
                rightpoint,t_rp = findrightBaseline(rawData,leftpoint,rightpoint,len(rawData)-10,end_thre)
                leftpoint,t_lp = findleftBaseline(rawData,leftpoint,rightpoint,10,str_thre)
                if t_lp <= str_thre and t_rp <= end_thre:
                    if -derivative_1[rightpoint+1] + slope(rawData,leftpoint,rightpoint+1)<=end_thre:
                        if -derivative_1[rightpoint+2] + slope(rawData,leftpoint,rightpoint+2)<=end_thre:
                            if -derivative_1[rightpoint+3] + slope(rawData,leftpoint,rightpoint+3)<=end_thre:
                                if -derivative_1[rightpoint+4] + slope(rawData,leftpoint,rightpoint+4)<=end_thre:
                                    if -derivative_1[rightpoint+5] + slope(rawData,leftpoint,rightpoint+5)<=end_thre:
                                        break
                if (rightpoint >= len(rawData) -11) or (leftpoint <= 11):
                    break
                #

            tempeak=[]
            tempeak.append(leftpoint)
            tempeak.append(peakTopIdx[i])
            tempeak.append(rightpoint)
            peak.append(tempeak)

    return peak

def mergeBasline(peak):
    '''
    函数：合并基线
    :param peak: array_like(存储单个峰特征的list，eg：[[基线起点，峰顶点，基线终点]，[基线起点,峰顶点，基线终点]，...])
    :return fusedpeaks: array_like(基线融合之后的list，eg:[[基线起点，[峰顶点],基线终点]，[基线起点，[峰顶点]，基线终点]，....])
    '''

    #按照起点大小来排序,eg:[[1,3],[5,8],[3,6]],res:[[1,3],[3,6],[5,8]]
    peak=sorted(peak,key = lambda x: x[0])
    start,end = peak[0][0] , peak[0][-1]

    temp_top=[peak[0][1]]
    fusedpeaks=[]

    for i in range(1,len(peak)):
        s,e = peak[i][0],peak[i][-1]
        if s<= end:
            end=max(end,e)
            temp_top.append(peak[i][1])
        else:
            temp_top.sort()
            fusedpeaks.append([start,temp_top,end])
            temp_top=[peak[i][1]]
            start,end = s,e

    fusedpeaks.append([start,temp_top,end])
    return fusedpeaks

def sloveFusedBasline(rawData,fusedpeaks):
    """
    函数：处理融合峰的基线，主要是检测是否有基线穿越情况
    :param rawData: array_like(需要计算基线的原数据)
    :param fusedpeaks: array_like(基线融合之后的list,eg:[[基线起点，[峰顶点],基线终点],[基线起点，[峰顶点]，基线终点])
    :return result : array_like(处理融合峰后的结果)
    """
    result = []
    
    # 从前往后扫描，出现基线穿越则断开，开始新的峰簇
    for group in fusedpeaks:
        group[1].sort()
        tempBslstart = group[0]
        tempBslend = group[-1]
        tempbsl = generatebasline(rawData, tempBslstart, tempBslend)
        temppeakpos = rawData[group[0] : group[-1] + 1] - tempbsl[0 : group[-1] - group[0] + 1]
        # 单峰
        if len(group[1]) == 1:

            newst = group[0] + temppeakpos[0 : group[1][0] - group[0] + 1].argmin()
            newed = group[1][0] + temppeakpos[group[1][0] - group[0] : group[2] - group[0] + 1].argmin()
            group[0], group[2] = newst, newed
            res = {"peak":[[group[0],group[1][0],group[2]]], "bsl":[newst, newed],'bslType':'auto'}
            result.append(res)
        else:
            # 融合峰：逐个峰扫描，若两峰之间最低点低于基线，则从最低点处断开基线
            # 断开基线前为一个融合峰，后面为另一个融合峰
            newst = group[0] + temppeakpos[0 : group[1][0] - group[0] + 1].argmin()
            group[0] = newst
            tempBslstart = group[0]

            newed = group[1][-1] + temppeakpos[group[1][-1] - group[0] : group[-1] - group[0] + 1].argmin()
            group[-1] = newed
            tempBslend = group[-1]
            
            tempbsl = generatebasline(rawData, tempBslstart, tempBslend)
            tempgroup = {"peak":[], "bsl":None,'bslType':'auto'}

            temptop=[]

            #首先找到融合基线段中的最低点
            globallow = group[1][0] + np.argmin(rawData[group[1][0]:group[1][-1]])
            #temppeakpos = rawData[tempstartcheck : group[-1][-1] + 1] - tempbsl[0 : group[-1][-1] - tempstartcheck + 1]
            #globallow = group[1][0] + (temppeakpos[group[1][0] - tempBslstart : group[1][-1] - tempBslstart + 1]).argmin()
            if rawData[globallow] < tempbsl[globallow-tempBslstart]:#最低点在基线下方，从中间断开，形成两段
                
                #将基线分成两段
                tempbslcoef = [tempBslstart, globallow, tempBslend]
                
                #找到基线中的峰顶点
                arraygroup=np.asarray(group[1])
                temptop.append((arraygroup[(arraygroup >= tempBslstart) & (arraygroup < globallow)]).tolist())
                temptop.append((arraygroup[(arraygroup > globallow) & (arraygroup <= tempBslend)]).tolist())

                nextlow=[]
                
                for i in range(2):#划分两段基线的起止点
                    tempgroup = {"peak":[], "bsl":None,'bslType':'auto'}
                    if i == 0:#前一段融合峰情况，划分基线的起止点
                        tempgroup["bsl"] = [tempbslcoef[0],tempbslcoef[1]]
                        start = tempbslcoef[0]
                        end = tempbslcoef[1]
                        bslstartcheck = start
                        firststart = start
                        tempbsl = generatebasline(rawData, start, end)
                    else:#后一段融合峰情况，划分基线的起止点
                        tempgroup["bsl"] = [tempbslcoef[1],tempbslcoef[2]]
                        start = tempbslcoef[1]
                        end = tempbslcoef[2]
                        bslstartcheck = start
                        firststart = start
                        tempbsl = generatebasline(rawData, start, end)

                    k=len(nextlow)
                    #寻找峰顶点之间的最低点
                    for j in range(len(temptop[i])-1):
                        #group[i][1] +(temppeakpos[group[i][1] - tempstartcheck : group[i+1][1] - tempstartcheck + 1]).argmin()
                        #peaklow = temptop[i][j] + (temppeakpos[temptop[i][j] - start : temptop[i][j+1] - start + 1]).argmin()
                        nextlow.append(temptop[i][j] + np.argmin(rawData[temptop[i][j]:temptop[i][j+1]]))
                        if rawData[nextlow[k]] <= tempbsl[nextlow[k]-start]:
                        #if rawData[peaklow] <= tempbsl[peaklow-start]:
                            #nextlow.append(peaklow)
                            tempgroup["peak"].append([bslstartcheck,temptop[i][j],nextlow[k]])
                            tempgroup["bsl"] = [firststart,nextlow[k]]
                            result.append(tempgroup)
                            tempgroup = {"peak":[], "bsl":None,'bslType':'auto'}
                            bslstartcheck = nextlow[k]
                            firststart= bslstartcheck
                        else :
                            #nextlow.append(temptop[i][j] + np.argmin(rawData[temptop[i][j]:temptop[i][j+1]]))
                            tempgroup["peak"].append([bslstartcheck,temptop[i][j],nextlow[k]])
                            #tempgroup["bsl"] = [firststart[0],end]
                            bslstartcheck = nextlow[k]
                        k+=1
                    if len(temptop[i]) == 1:#如果只有一个峰顶点
                        res = {"peak":[[start,temptop[i][-1],end]], "bsl":[start,end],'bslType':'auto'}
                        result.append(res)
                        firststart=end
                    else:#将最后一个峰加进去
                        tempgroup["peak"].append([nextlow[-1],temptop[i][-1],end])
                        tempgroup["bsl"] = [firststart,end]
                        result.append(tempgroup)
                        #firststart=[]
                        tempgroup = {"peak":[], "bsl":None,'bslType':'auto'}
            else:#最低点在基线上方
                
                tempbslcoef = [tempBslstart, tempBslend]
                bslstartcheck = tempBslstart
                firststart = tempBslstart

                nextlow=[]#峰顶点之间的最低点
                for i in range(len(group[1])-1):
                        #nextlow.append(group[1][i] + (temppeakpos[group[1][i] - tempBslstart : group[1][i+1] - tempBslstart + 1]).argmin())#峰顶点中的最低点
                        nextlow.append(group[1][i] + np.argmin(rawData[group[1][i]:group[1][i+1]]))#峰顶点中的最低点
                        #peaklow = group[1][i] + (temppeakpos[group[1][i] - tempBslstart : group[1][i+1] - tempBslstart + 1]).argmin()
                        if rawData[nextlow[i]] <= tempbsl[nextlow[i]-tempBslstart]:
                        #if rawData[peaklow] <= tempbsl[peaklow-tempBslstart]:
                            #nextlow.append(peaklow)
                            tempgroup["peak"].append([bslstartcheck,group[1][i],nextlow[i]])
                            tempgroup["bsl"] = [firststart,nextlow[i]]
                            result.append(tempgroup)
                            tempgroup = {"peak":[], "bsl":None,'bslType':'auto'}
                            bslstartcheck = nextlow[i]
                            firststart=bslstartcheck
                            
                        else :
                            #nextlow.append(group[1][i] + np.argmin(rawData[group[1][i]:group[1][i+1]]))
                            
                            tempgroup["peak"].append([bslstartcheck,group[1][i],nextlow[i]])
                            #tempgroup["bsl"] = [firststart[0],end]
                            bslstartcheck = nextlow[i]
                            #firststart.append(bslstartcheck)
                            
                tempgroup["peak"].append([nextlow[-1],group[1][-1],tempBslend])
                tempgroup["bsl"] = [firststart,tempBslend]
                result.append(tempgroup)
                #firststart=[]
                tempgroup = {"peak":[], "bsl":None,'bslType':'auto'}

    return result


def getApexPreResult(rawData, timeData,timeInterval, startTime, endTime, peakWidth = None,apexThreshold=None,startThreshold = 0.,endThreshold=0.005 ):

    start = int((startTime - timeData[0])//timeInterval)
    end = int((endTime - timeData[0])//timeInterval)
    dataFreq = (1 / timeInterval)
    if start < 50:
        start = 50
    if end > len(rawData)-50:
        end = len(rawData)
        # end = len(rawData)-50
    Integral_zone = [start,end]

    if peakWidth != None:
        peakWidth = (int(peakWidth / timeInterval) // 2) * 2 + 1
        filter_par2 = peakWidth
    else:
        filter_par2 = filter_parameter(rawData)  # 计算二阶过滤系数
# 计算峰顶点的阈值 如果存在
    if apexThreshold == None:
        derivative_1 = np.diff(rawData,1)
        tem_derivative_2=-np.diff(derivative_1.T,1) 
        derivative_2=np.hstack((tem_derivative_2[0:2],tem_derivative_2)) #补全最前面缺少的两个二阶导数点
        filter_der2=signal.savgol_filter(derivative_2, filter_par2, 2) #滤波后的二阶导数
        apexThreshold = autoapexThreshold(filter_der2,dataFreq)
    else:
        apexThreshold = apexThreshold
    rawDataSliced = runpreFilter(rawData[start:end])
    peakTopIdx = runApexDetection(rawDataSliced,Integral_zone,filter_par2,apexThreshold)
    if len(peakTopIdx) == 0:
        err =-60
        return err,None
    peak = solveBaseline(rawData,peakTopIdx,startThreshold,endThreshold)
    fusedpeaks = mergeBasline(peak)
    pre_result = sloveFusedBasline(rawData,fusedpeaks)

    #print(pre_result)
    return 0, pre_result

# def RunApexIntegral(rawData, timeData,timeInterval,peakWidth = None,apexThreshold=None,startThreshold = 0.,endThreshold=0.005):
#     #根据输入参数运行Apex积分
#
#     dataFreq = int(round(1/timeInterval))
#     Integral_zone = [100,len(rawData)-100]
#
#     if apexThreshold == None:
#         derivative_1 = np.diff(rawData,1)
#         tem_derivative_2=-np.diff(derivative_1.T,1)
#         derivative_2=np.hstack((tem_derivative_2[0:2],tem_derivative_2)) #补全最前面缺少的两个二阶导数点
#         if peakWidth!=None:
#             peakWidth = (int(peakWidth/timeInterval)//2)*2+1
#             filter_par2=peakWidth
#         else:
#             filter_par2=filter_parameter(rawData) #计算二阶过滤系数
#         filter_der2=signal.savgol_filter(derivative_2, filter_par2, 2) #滤波后的二阶导数
#         apexThreshold = autoapexThreshold(filter_der2)/3
#
#     rawDataSliced = runpreFilter(rawData[100:len(rawData)-100])
#     peakTopIdx = runApexDetection(rawDataSliced,Integral_zone,filter_par2,apexThreshold)
#     peak = solveBaseline(rawData,peakTopIdx,startThreshold,endThreshold)
#     fusedpeaks = mergeBasline(peak)
#     pre_result = sloveFusedBasline(rawData,fusedpeaks)
#
#
#     result = calculateResult(pre_result, rawData, timeData, dataFreq)
#     #print(result)
#     return result