import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import integrate
import pandas as pd
from ..basic import base1 as BasicPart
from ..event_apex import event as pe

    
def find_apex(raw_data, Integral_zone,width, threshold = None):
    """
    函数：Apex峰顶点定位(输入为数组，一次处理所有数据)
    :param raw_data: array_like(用于计算apex峰顶点的原数据)
    :return: apex_index : array_like(所有apex峰顶点的index)
    :return: filtered_2 : array_like(经过滤波后的二阶导数)
    """
    # 计算需要使用的数据
    left_zone=Integral_zone[::2]
    right_zone=Integral_zone[1::2]
    get_i=raw_data
    numal=len(get_i)
    flt_par2 = BasicPart.filter_parameter(get_i)
    if width!=None:
        flt_par2=width
    #print(flt_par2)
    deri_1 = np.diff(get_i,1 )
    filtered_1 = deri_1
    deri_2 = -np.diff(filtered_1.T, 1)
    #print(len(deri_2))
    deri_2=np.hstack((deri_2[0:1],deri_2))
    #print(len(deri_2))
    filtered_2 = signal.savgol_filter(deri_2, flt_par2, 2)
    if threshold is not None:
        if threshold <= 0:
            thre_2 = BasicPart.threshold(filtered_2, 20, 136)
            thre= BasicPart.threshold(get_i, 20, 136)
        else:
            thre_2 = threshold
    else:
        thre_2 = BasicPart.threshold(filtered_2, 20, 136)
        thre= BasicPart.threshold(get_i, 20, 136)
            
    arg_order = int(flt_par2 / 5)
    if arg_order < 5:
        arg_order = 5
    temp = signal.argrelmax(filtered_2, order=arg_order)
    apex_index = temp[0]
    media_index=apex_index[apex_index<100]
    
    apex_index =apex_index[(apex_index<len(raw_data)-100)&(apex_index>100)]
    #diffapex=np.diff(apex_index)
    #threindex=np.argmax(diffapex)#[(apex_index<len(raw_data)-100)&(apex_index>100)]
    #thre_2=BasicPart.threshold(filtered_2,apex_index[threindex]+3*diffapex[threindex]/8,apex_index[threindex+1]-3*diffapex[threindex]/8)
    apex_peak = []
    i=0
    minlist=[]
    #minrlist=[]
    '''
    while True:
        #for i in range(1,len(apex_index)-1):
        if i==0:
            minl=apex_index[i]-10+np.argmin(get_i[apex_index[i]-10:apex_index[i]])
            minr=apex_index[i]+np.argmin(get_i[apex_index[i]:apex_index[i+1]])
        elif i == len(apex_index)-1:
            minl=apex_index[i-1]+np.argmin(get_i[apex_index[i-1]:apex_index[i]])
            minr=apex_index[i]+np.argmin(get_i[apex_index[i]:apex_index[i]+10])
        else:
            minl=apex_index[i-1]+np.argmin(get_i[apex_index[i-1]:apex_index[i]])
            minr=apex_index[i]+np.argmin(get_i[apex_index[i]:apex_index[i+1]])
        useb=minl+np.argmax(get_i[minl:minr])
        if ((useb>left_zone)&(useb<right_zone)).any():
            if useb-minl>3 and minr-useb>3:
                #if apex_index[i]-min(get_i[minl],get_i[minr])> thre:
                    apex_peak.append(useb)
                    minlist.append(minl)
                    minlist.append(minr)
                    #minrlist.append(minr)
            else :
                apex_index=np.delete(apex_index,i)
                i-=1
        i+=1
        if i>=len(apex_index)-1:
            break
    threlist= np.diff(np.array(minlist),1)[1::2]
    threindex=np.argmax(threlist)*2+1
    thre_2=BasicPart.threshold(filtered_2, minlist[threindex]+threlist.max()/4, minlist[threindex+1]-threlist.max()/4)
    apex_peak=np.array(apex_peak)
    apex_peak=apex_peak[filtered_2[apex_peak]>thre_2]
    '''
    while True:
    #for i in range(1,len(apex_index)-1):
        if filtered_2[apex_index[i]] >thre_2/4:
            if i==0:
                minl=apex_index[i]-10+np.argmin(get_i[apex_index[i]-10:apex_index[i]])
                minr=apex_index[i]+np.argmin(get_i[apex_index[i]:apex_index[i+1]])
            elif i == len(apex_index)-1:
                minl=apex_index[i-1]+np.argmin(get_i[apex_index[i-1]:apex_index[i]])
                minr=apex_index[i]+np.argmin(get_i[apex_index[i]:apex_index[i]+10])
            else:
                minl=apex_index[i-1]+np.argmin(get_i[apex_index[i-1]:apex_index[i]])
                minr=apex_index[i]+np.argmin(get_i[apex_index[i]:apex_index[i+1]])
            useb=minl+np.argmax(get_i[minl:minr])
            if ((useb>left_zone)&(useb<right_zone)).any():
                if useb-minl>3 and minr-useb>3:
                    #if apex_index[i]-min(get_i[minl],get_i[minr])> thre:
                        apex_peak.append(useb)
                else :
                    apex_index=np.delete(apex_index,i)
                    i-=1
        i+=1
        if i>=len(apex_index)-1:
            break
    #'''
    '''
    for i in range(len(re_index)-1):#apex_index:
        if filtered_2[re_index[i]]<0:
            if -filtered_2[re_index[i]] > thre_2:
        #print(re_index[i])
                useb=re_index[i]+np.argmax(get_i[re_index[i]:re_index[1+i]])+1
        #print(useb)
                if filtered_2[useb] > thre_2:
                    apex_peak.append(useb)
        
    for i in apex_index:
        #d60=[t for t in range(max(i-1,0),min(i+2,numal))]
        if filtered_2[i] > thre_2:
        
            
            maxnumber=min(i+(arg_order),numal-1)
            minnumber=max(i-(arg_order),0)
            zone=np.arange(minnumber,maxnumber)
            useb=i+np.argmax(raw_data[zone])-arg_order
            #print(useb,maxnumber,minnumber)
            if useb==maxnumber or useb==minnumber:
                #print("t")
                maxnumber=min(i+2*arg_order,numal-1)
                minnumber=max(i-2*arg_order,0)
                zone=np.arange(minnumber,maxnumber)
                useb=i+np.argmax(raw_data[zone])-arg_order
                if useb==maxnumber or useb==minnumber:
                    continue
            #if useb!=i:
                #print(useb,i)
            useb=max(useb,0)
            useb=min(useb,numal-1)
            apex_peak.append(useb)
            '''
            #print(apex_peak)
            
                #print(apex_peak)
                #re_index[abs(re_index-i)>arg_order]
            #useb=np.argmax(raw_data[d60])
            #print(useb)
            #if raw_data[useb+i-1]>=raw_data[max(useb+i-1-1,0)]:
                #if raw_data[useb+i-1]>=raw_data[min(useb+i+1-1,numal)]:
            #if raw_data[i]-raw_data[min(i+2,numal)]>=0 or raw_data[i+1]-raw_data[min(i+3,numal)]>=0:
                #if raw_data[i]-raw_data[max(i-1,0)]>=0 or raw_data[i-1]-raw_data[max(i-3,0)]>=0:
                    #if raw_data[min(i+1,numal)]>raw_data[min(i+2,numal)]:
                        #if raw_data[min(i-1,numal)]>raw_data[min(i-2,numal)]:
            #apex_peak.append(i)
            # _, _, inf_width = BasicPart.infl(raw_data, i, filtered_2)
            # if inf_width >= 3 and raw_data[i] > thre:
            #apex_peak.append(i)

    apex_peak=(list(set(apex_peak)))
    apex_peak.sort()
    return apex_peak, filtered_2

def funtemp_r(raw_data,lp,rp,rmax,end_thre):
    deri_11= np.diff(raw_data, 1)
    #print(rmax)
    k,_ = BasicPart.slope(raw_data, lp, rp)
    r=rp
    for j in range(rp,rmax,1):  # 预备基线右端点向右延伸，直到满足阈值要求
        k, _ = BasicPart.slope(raw_data, lp, j)
        temp_r = -deri_11[j] + k
        if temp_r <= end_thre:
            if -deri_11[j+1] + BasicPart.slope(raw_data, lp, j+1)[0] <= end_thre:
                if -deri_11[j+2] + BasicPart.slope(raw_data, lp, j+2)[0] <= end_thre:
                    r=j
                    #print('ddd')
                    break
        if j>= rmax-1:
            r=j
            #print('sm')
            break
    return k,r

def funtemp_l(raw_data,lp,rp,lmin,str_thre):
    deri_11= np.diff(raw_data, 1)
    k,_ = BasicPart.slope(raw_data, lp, rp)
    l=lp
    for i in range(lp,lmin, -1):  # 右端点满足后延伸左端点，直到满足阈值要求
        #print(i)
        k, _ = BasicPart.slope(raw_data, i, rp)
        temp_l = deri_11[i] - k
        if temp_l <= str_thre:
            if deri_11[i-1] - BasicPart.slope(raw_data, i-1,rp)[0]  <= str_thre:
                if deri_11[i-2] - BasicPart.slope(raw_data,i-2, rp)[0]  <= str_thre:
                    l=i
                    break
        if i <=lmin+1:
            l=i
            #print('em')
            break
    return k,l

def apex_bsl(zone_data, x, name, d_2, step=1., liftoff=0.0, touchdown=0.):
    """
    函数：用于第一次确定基线，由Apex点左右的拐点连线作为预备基线(输入为数组，一次仅处理一个峰)
    :param a: array_like(需要计算基线的原数据)
    :param x: int(需要计算基线的Apex点index)
    :param name: int(峰编号，用作标识，函数中不做处理)
    :param step: float, optional(计算终点阈值的参数，默认值为1)
    :param d_2: array_like(二阶导数)
    :param liftoff: float, optional(峰起点阈值，默认值为0)
    :param touchdown: float, ptional(峰终点阈值，默认值为0)
    :return Apex峰数据 : DataFrame('name', 'start', 'end', 'baseline start', 'baseline end', 'inflection left',
                                   'inflection right', 'Apex peak', 'slope', 'drift','peak type1', 'peak type2')
    """
    # 函数内部使用的导数
    deri_11 = np.diff(zone_data, 1)
    data_num = len(zone_data)-1
    # 一些初始需要计算的值
    infl_l, infl_r, _ = BasicPart.infl(zone_data, x, d_2)  # inf值需保留返回以供后续使用
    lp = infl_l
    rp = infl_r
    ini_slp, _ = BasicPart.slope(zone_data, lp, rp)
    str_thre = np.tan((np.arctan(deri_11[lp]) - np.arctan(ini_slp)) * liftoff)  # 将斜率化为角度计算后再化为斜率
    end_thre = np.tan((np.arctan(-deri_11[rp]) + np.arctan(ini_slp)) * touchdown)
    # 循环主体
            #elif a[max(i-1,0)]-a[i]>0:
                #if a[max(i-2,0)]-a[max(i-1,0)]<0:
                    #continue
                #if a[max(i-5,0)]-a[max(i-4,0)]<0:
                    #continue
                #print(4)
                #return k,i
    con=0
    while True: #con<=5:
        con+=1
        if con>10:
            print('f',rp,lp)
        k,rp = funtemp_r(zone_data,lp,rp,data_num-10,end_thre)
        k,lp = funtemp_l(zone_data,lp,rp,10,str_thre)
        temp_r = -deri_11[rp] + k
        if temp_r <= end_thre:
            if -deri_11[rp+1] + BasicPart.slope(zone_data, lp, rp+1)[0] <= end_thre:
                if -deri_11[rp+2] + BasicPart.slope(zone_data, lp, rp+2)[0] <= end_thre:
                    break
        if (rp>= data_num -11) or (lp <=11):
            break
        #t+=1
        #elif a[min(j+1,data_num - 1)]-a[j]>0:
            #if a[min(j+2,data_num - 1)]-a[min(j+1,data_num - 1)]<0:
                    #continue
            #if a[min(j+5,data_num - 1)]-a[min(j+4,data_num - 1)]<0:
                    #continue
            #print(6)
            #break
    slp, dft = BasicPart.slope(zone_data, lp, rp, step)
    return name, lp, rp, lp, rp, infl_l, infl_r, x, slp, dft, 'B', 'B',zone_data[lp],zone_data[rp]
    '''
    for j in range(j, data_num - 1, 1):  # 预备基线右端点向右延伸，直到满足阈值要求
        k, _ = BasicPart.slope(a, i, j)
        temp_r = -deri_11[j] + k
        if temp_r <= end_thre:  # 右端点判定
            if -deri_11[min(j+1,data_num - 2)] + BasicPart.slope(a, i, min(j+1,data_num - 2))[0]> end_thre:
                continue
            if -deri_11[min(j+2,data_num - 2)] + BasicPart.slope(a, i, min(j+2,data_num - 2))[0]> end_thre:
                continue
            if -deri_11[min(j+5,data_num - 2)] + BasicPart.slope(a, i, min(j+5,data_num - 2))[0]> end_thre:
                continue
            for i in range(i, 0, -1):  # 右端点满足后延伸左端点，直到满足阈值要求
                k, _ = BasicPart.slope(a, i, j)
                temp_l = deri_11[i] - k

                if temp_l <= str_thre:  # 左端点判定
                    if deri_11[max(i-1,0)] - BasicPart.slope(a, max(i-1,0), j)[0] > str_thre:
                        continue
                    if deri_11[max(i-2,0)] - BasicPart.slope(a, max(i-2,0), j)[0] > str_thre:
                        continue
                    if deri_11[max(i-5,0)] - BasicPart.slope(a, max(i-5,0), j)[0] > str_thre:
                        continue
                    k, _ = BasicPart.slope(a, i, j)
                    temp_r = -deri_11[j] + k

                    if temp_r <= end_thre:  # 左端点满足后再次检查右端点
                        if -deri_11[min(j+1,data_num - 2)] + BasicPart.slope(a, i, min(j+1,data_num - 2))[0] > end_thre:
                            continue
                        if -deri_11[min(j+2,data_num - 2)] + BasicPart.slope(a, i, min(j+2,data_num - 2))[0] > end_thre:
                            continue
                        if -deri_11[min(j+5,data_num - 2)] + BasicPart.slope(a, i, min(j+5,data_num - 2))[0]> end_thre:
                            continue
                        slp, dft = BasicPart.slope(a, i, j, step)
                        #print(name, i, j, i, j, infl_l, infl_r, x, slp, dft, 'B', 'B')
                        return name, i, j, i, j, infl_l, infl_r, x, slp, dft, 'B', 'B'
                    else:
                        break  # 右端点不满足则继续向左侧延伸寻找基线左端点
    '''
def apex_data_generator(raw_data, pre_peak, d_2,step=1., liftoff=0.0, touchdown=0.):
    """
    生成包含['name', 'start', 'end', 'baseline start', 'baseline end', 'inflection left','inflection right',
    'Apex peak', 'slope', 'drift','peak type1', 'peak type2']的ApexData
    :param raw_data: array_like(原数据)
    :param pre_peak: array_like(所有apex峰顶点)
    :param d_2: array_like(二阶导数)
    :param step: float, optional(原数据时间间隔，默认值为1)
    :param liftoff: float, optional(峰起点阈值，默认值为0)
    :param touchdown: float, optional(峰终点阈值，默认值为5)
    :return ApexData: DataFrame(包含所有信息的DataFramne)
    """
    name = 0
    apex_data = []
    for x in pre_peak:
        temp = apex_bsl(raw_data, x, name, d_2, step, liftoff, touchdown)
        if temp is not None:
            apex_data.append(temp)
            name += 1

    # print('apex_data', apex_data)
    ApexData = pd.DataFrame(apex_data, columns=['name', 'str', 'end', 'bslstr', 'bslend', 'infl','infr', 'Apex', 'slp', 'dft', 'type1', 'type2','bslstrh','bslendh'])
    return ApexData

def apex_data_modify(raw_data, ApexData):
    media=[]
    length = len(ApexData)
    #NewData=ApexData.copy()
    for i in range(1,length-1):
        #media.append(ApexData.iloc[i,:])
        if (ApexData.iloc[i+1, 1]-ApexData.iloc[i, 2]<3) and (ApexData.iloc[i, 2]-ApexData.iloc[i, 7]<=5):
            ApexData.iloc[i+1, 1]=ApexData.iloc[i, 1]
            ApexData.iloc[i+1, 3]=ApexData.iloc[i, 3]
            ApexData.iloc[i+1, 11] =ApexData.iloc[i, 11]
            ApexData.iloc[i+1,7]=np.argmax(raw_data[ApexData.iloc[i+1, 1]:ApexData.iloc[i+1, 2]])+ApexData.iloc[i+1, 1]
            #NewData=NewData.drop(i)
            media.append(i)
        elif (ApexData.iloc[i, 1]-ApexData.iloc[i-1, 2]<3) and (ApexData.iloc[i, 7]-ApexData.iloc[i, 1]<=5):
            #print(ApexData.iloc[i-1, 2],ApexData.iloc[i, 2])
            ApexData.iloc[i-1, 2]=ApexData.iloc[i, 2]
            ApexData.iloc[i-1, 4]=ApexData.iloc[i, 4]
            ApexData.iloc[i-1, 10] =ApexData.iloc[i, 10]
            ApexData.iloc[i-1,7]=np.argmax(raw_data[ApexData.iloc[i-1, 1]:ApexData.iloc[i-1, 2]])+ApexData.iloc[i-1, 1]
            #media[-1]=ApexData.iloc[i-1,:]
            media.append(i)
            #print(len(media),i)
            #print(len(media),i)
        #'''
        if ApexData.iloc[i-1, 7]==ApexData.iloc[i, 7]:
            ApexData.iloc[i-1, 2]=ApexData.iloc[i, 2]
            ApexData.iloc[i-1, 4]=ApexData.iloc[i, 4]
            ApexData.iloc[i-1, 10] =ApexData.iloc[i, 10]
            ApexData.iloc[i-1,7]=np.argmax(raw_data[ApexData.iloc[i-1, 1]:ApexData.iloc[i-1, 2]])+ApexData.iloc[i-1, 1]
            media.append(i)
    #ApexData=pd.concat(media,axis=1).T
    media=list(set(media))
    for i in media:
        ApexData=ApexData.drop(i)
    ApexData.index=range(len(ApexData))
    ApexData['name']=range(len(ApexData))
    if ApexData.iloc[-1, 7]==ApexData.iloc[-2, 7]:
        ApexData.iloc[-1, 2]=ApexData.iloc[-2, 2]
        ApexData.iloc[-1, 4]=ApexData.iloc[-2, 4]
        ApexData.iloc[-1, 10] =ApexData.iloc[-2, 10]
        ApexData.iloc[-1,7]=np.argmax(raw_data[ApexData.iloc[-1, 1]:ApexData.iloc[-1, 2]])+ApexData.iloc[-1, 1]
    ApexData=ApexData.drop(len(ApexData)-2)
    ApexData.index=range(len(ApexData))
    ApexData['name']=range(len(ApexData))
    for i in range(len(ApexData)):
        st = ApexData.iloc[i, 1]
        en = ApexData.iloc[i, 2]
        bs = ApexData.iloc[i, 3]
        be = ApexData.iloc[i, 4]
        apex = ApexData.iloc[i, 7]
        mn=np.argmin(raw_data[apex:be+1])+apex+1
        ms=np.argmin(raw_data[bs:apex])+bs-1
        ApexData.iloc[i, 4] = mn
        ApexData.iloc[i, 2]=mn
        ApexData.iloc[i, 3] = ms
        ApexData.iloc[i, 1]=ms
    #print(Apex_peak)
    return ApexData

def fuse_bsl(a, l, r, il, ir, step=1., liftoff=0., touchdown=0.0):
    """
    函数：用于再次确定基线，由融合峰的峰族最左侧子峰的左拐点和峰族最右侧子峰的右拐点连线作为预备基线
    :param a: array_like(需要计算基线的原数据)
    :param l: int(峰族最左侧子峰的起点index)
    :param r: int(峰族最右侧子峰的终点index)
    :param il: int(峰族最左侧子峰的左拐点index)
    :param ir: int(峰族最右侧子峰的右拐点index)
    :param step: float, optional(数据点时间间隔，默认值为1)
    :param liftoff: float, optional(计算起点阈值的参数，默认值为0)
    :param touchdown: float, optional(计算终点阈值的参数，默认值为0)
    :return i : int(峰族基线起点index)
    :return j : int(峰族基线终点index)
    :return slp : float(峰族基线斜率)
    :return dft : float(峰族基线漂移)
    :return flag : Bool(是否需重新定位基线)
    """
    data_num = len(a)
    deri_11 = np.diff(a, 1)

    # 一些初始需要计算的值
    lp = l
    rp = r
    ini_slp, _ = BasicPart.slope(a, il, ir)  # 最外侧峰的左右拐点连线作为预备基线
    str_thre = np.tan((np.arctan(deri_11[il]) - np.arctan(ini_slp)) * liftoff)  # 将斜率化为角度计算后再化为斜率
    end_thre = np.tan((np.arctan(-deri_11[ir]) + np.arctan(ini_slp)) * touchdown)
    slp, dft = BasicPart.slope(a, lp, rp)
    ini_s = np.tan(np.arctan(deri_11[lp]) - np.arctan(slp))
    ini_e = np.tan(np.arctan(-deri_11[rp]) + np.arctan(slp))

    # 判断是否满足条件，满足时直接返回原值
    if ini_s <= str_thre and ini_e <= end_thre:
        # print('pass')
        return l, r, slp, dft, False

    # 不满足条件时进入循环重新定位基线
    '''
    def funtemp_r(j):
        for j in range(j, data_num - 1, 1):  # 预备基线右端点向右延伸，直到满足阈值要求
            k, _ = BasicPart.slope(a, i, j)
            temp_r = -deri_11[j] + k
            if temp_r <= end_thre:
                return k,j
    def funtemp_l(i):
        for i in range(i, 0, -1):  # 右端点满足后延伸左端点，直到满足阈值要求
            k, _ = BasicPart.slope(a, i, j)
            temp_l = deri_11[i] - k
            if temp_l <= str_thre:
                return k,i
    while t<=5:
        k,j = funtemp_r(j)
        k,i = funtemp_l(i)
        temp_r = -deri_11[j] + k
        t+=1
        if temp_r <= end_thre:
            break
    slp, dft = BasicPart.slope(a, i, j, step)
    return i, j, slp, dft, True
    '''
    con=0
    while True:#con<=5:
        con+=1
        if con>10:
            print('s',rp,lp)
        k,rp = funtemp_r(a,lp,rp,len(deri_11) -10,end_thre)
        k,lp = funtemp_l(a,lp,rp,10,str_thre)
        temp_r = -deri_11[rp] + k
        if temp_r <= end_thre:
            if -deri_11[rp+1] + BasicPart.slope(a, lp, rp+1)[0] <= end_thre:
                if -deri_11[rp+2] + BasicPart.slope(a, lp, rp+2)[0] <= end_thre:
                    break
        if (rp>= len(deri_11)-11) or (lp <=11):
            break
    slp, dft = BasicPart.slope(a, lp, rp, step)
    return lp, rp, slp, dft, True
    '''
    for j in range(j, data_num - 1, 1):
        if j == data_num - 2:
            #print('nooooo')
            return i, j, slp, dft, True
        # print(j)

        k, _ = BasicPart.slope(a, i, j)
        temp_r = -deri_11[j] + k

        if temp_r <= end_thre:  # 右侧条件满足
            #print('1')
            for i in range(i, 0, -1):

                k, _ = BasicPart.slope(a, i, j)
                temp_l = deri_11[i] - k

                if temp_l <= str_thre:  # 左侧条件满足
                    #print('2')
                    k, _ = BasicPart.slope(a, i, j)
                    temp_r = -deri_11[j] + k

                    if temp_r <= end_thre:  # 检查右侧条件，满足时输出结果，不满足时进入循环
                        #print('3')
                        #print(i, j)
                        slp, dft = BasicPart.slope(a, i, j, step)
                        #print('finish')
                        return i, j, slp, dft, True
                    else:
                        break
    
    for j in range(j, data_num - 1, 1):  # 预备基线右端点向右延伸，直到满足阈值要求
        k, _ = BasicPart.slope(a, i, j)
        temp_r = -deri_11[j] + k
        if temp_r <= end_thre:  # 右端点判定
            if -deri_11[min(j+1,data_num - 1)] + BasicPart.slope(a, i, min(j+1,data_num - 1))[0]<= end_thre:
                continue
            if -deri_11[min(j+2,data_num - 1)] + BasicPart.slope(a, i, min(j+2,data_num - 1))[0]<= end_thre:
                continue
            if -deri_11[min(j+5,data_num - 1)] + BasicPart.slope(a, i, min(j+5,data_num - 1))[0]<= end_thre:
                continue
            for i in range(i, 0, -1):  # 右端点满足后延伸左端点，直到满足阈值要求
                k, _ = BasicPart.slope(a, i, j)
                temp_l = deri_11[i] - k

                if temp_l <= str_thre:  # 左端点判定
                    if deri_11[max(i-1,0)] - BasicPart.slope(a, max(i-1,0), j)[0] <= str_thre:
                        continue
                    if deri_11[max(i-2,0)] - BasicPart.slope(a, max(i-2,0), j)[0] <= str_thre:
                        continue
                    if deri_11[max(i-5,0)] - BasicPart.slope(a, max(i-5,0), j)[0] <= str_thre:
                        continue
                    k, _ = BasicPart.slope(a, i, j)
                    temp_r = -deri_11[j] + k

                    if temp_r <= end_thre:  # 左端点满足后再次检查右端点
                        if -deri_11[min(j+1,data_num - 1)] + BasicPart.slope(a, i, min(j+1,data_num - 1))[0]<= end_thre:
                            continue
                        if -deri_11[min(j+2,data_num - 1)] + BasicPart.slope(a, i, min(j+2,data_num - 1))[0]<= end_thre:
                            continue
                        if -deri_11[min(j+5,data_num - 1)] + BasicPart.slope(a, i, min(j+5,data_num - 1))[0]<= end_thre:
                            continue
                        slp, dft = BasicPart.slope(a, i, j, step)
                        
                        return i, j, slp, dft, True
                    else:
                        break  # 右端点不满足则继续向左侧延伸寻找基线左端点
    '''

def fuse_operate(ApexData, raw_data, d2, step=1., liftoff=0., touchdown=0.):
    """
    处理未进行峰融合的DataFrame，输出包含处理完全数据的DataFrame
    :param ApexData: DataFrame(初步寻找Apex峰顶点后的DataFrame)
    :param raw_data: array_like(原数据)
    :param d2: array_like(原数据二阶导数，用于峰分割)
    :param step: float, optional(数据点时间间隔，默认值为1)
    :param touchdown: float, optional(起始阈值，默认值为0)
    :param liftoff: float, optional(末尾阈值，默认值为0)
    :return ApexData: DataFrame(处理后的DataFrame)
    """
    t=0
    length = len(ApexData)  # DataFrame的长度不会改变
    while True:
        t+=1
        #print(t)
        # 判断
        fuse_num = 0
        for i in range(1, length):
            b1 = ApexData.iloc[i - 1, 4]
            a2 = ApexData.iloc[i, 3]
            if t>20:
                print(b1,a2)
            if b1 > a2:  # 前一个峰的基线终点在后一个峰的基线起点后面视为基线交叉
                if (ApexData.iloc[i - 1, 11] != 'V' )|(ApexData.iloc[i, 10] != 'V'):
                    fuse_num += 1
                    # 更新峰类型
                    ApexData.iloc[i - 1, 11] = 'V'
                    ApexData.iloc[i, 10] = 'V'
                    # 更新前峰终点后峰起点
                    p1 = ApexData.iloc[i - 1, 7]
                    p2 = ApexData.iloc[i, 7]
                    #print(p1,p2)
                    temp = pe.skim_vertical(raw_data, d2, p1, p2)
                    ApexData.iloc[i - 1, 2] = temp
                    ApexData.iloc[i, 1] = temp

        if fuse_num == 0:
            # 重新计算斜率及漂移
            bs = ApexData.iloc[:, 3]
            be = ApexData.iloc[:, 4]
            for i in range(len(ApexData)):
                slp_temp, dft_temp = BasicPart.slope(raw_data, bs[i], be[i], step)
                ApexData.iloc[i, 8] = slp_temp
                ApexData.iloc[i, 9] = dft_temp
            return ApexData

        # 根据基线交叉情况更新完毕峰类型后找到所有峰族并更新峰基线起终点
        cluster = pe.search_cluster(ApexData, 10, 11)

        # 在每个峰族中根据重新定位的基线更新峰族起始峰和结尾峰的起点终点
        for i in cluster:
            if len(i) > 1:
                str_p = i[0]
                end_p = i[1]
                # 基线扩展
                cluster_s = ApexData.iloc[str_p, 1]
                cluster_e = ApexData.iloc[end_p, 2]
                cluster_is = ApexData.iloc[str_p, 5]
                cluster_ie = ApexData.iloc[end_p, 6]
                #print('start with:', cluster_s, cluster_e, cluster_is, cluster_ie)
                new_bs, new_be, new_slp, new_dft, flag = fuse_bsl(raw_data, cluster_s, cluster_e, cluster_is, cluster_ie,step, liftoff, touchdown)
                # 更新峰信息
                if flag is True:  # 基线扩展后更新峰族起始峰和结尾峰的起点终点
                    # print('flag =', flag)
                    ApexData.iloc[str_p: end_p + 1, 3] = new_bs
                    ApexData.iloc[str_p: end_p + 1, 4] = new_be
                    ApexData.iloc[str_p: end_p + 1, 8] = new_slp
                    ApexData.iloc[str_p: end_p + 1, 9] = new_dft
                    ApexData.iloc[str_p, 1] = new_bs
                    ApexData.iloc[end_p, 2] = new_be
                    ApexData.iloc[str_p, 12]=raw_data[new_bs]
                    ApexData.iloc[str_p, 13]=raw_data[new_be]

                else:
                    ApexData.iloc[str_p: end_p + 1, 3] = new_bs
                    ApexData.iloc[str_p: end_p + 1, 4] = new_be
                    ApexData.iloc[str_p: end_p + 1, 8] = new_slp
                    ApexData.iloc[str_p: end_p + 1, 9] = new_dft

def V_to_V(ApexData,raw_data,timelist,step,liftoff=0., touchdown=0.0):
    if timelist[0]>len(raw_data):
        return ApexData
    index=ApexData.index
    timelist=[int(b) for b in np.array(timelist)/step]
    st = ApexData.iloc[:, 1]
    en = ApexData.iloc[:, 2]
    bs = ApexData.iloc[:, 3]
    be = ApexData.iloc[:, 4]
    apex=ApexData.iloc[:,7]
    data_num = len(raw_data)
    deri_11 = np.diff(raw_data, 1)
    media=[]
    ApexData.index=apex
    if len(timelist)%2 != 0:
        timelist.append(int(ApexData.iloc[-1,7]/step)+1)
    if timelist[0]>0:
        mediadata=ApexData.loc[0:timelist[0],:].copy()
        mediadata.index=mediadata.iloc[:,0]
        media.append(mediadata)
    for i in range(len(timelist)//2):
        if i!=0:
            if timelist[2*i-1]+1>=apex.iloc[-1]:
                last=timelist[2*i-1]+1
                break
            mediadata=ApexData.loc[timelist[2*i-1]+1:timelist[2*i],:].copy()
            mediadata.index=mediadata.iloc[:,0]
            media.append(mediadata)
        mediadata=ApexData.loc[timelist[2*i]+1:timelist[2*i+1],:].copy()
        if len(mediadata)==0:
            continue
        num=len(mediadata)
        mediadata.index=mediadata.iloc[:,0]
        cluster = pe.search_cluster(mediadata, 10, 11)
        for each in cluster:
            if len(each) > 1:
                for t in range(each[0],each[1]+1):
                    #print(mediadata)
                    #print(halfleft,halfright)
                    lpoint=mediadata.loc[t,'infl']
                    rpoint=mediadata.loc[t,'infr']
                    strpoint=mediadata.loc[t,'str']
                    if mediadata.loc[t,'str']!=mediadata.loc[t,'Apex']:
                        strpoint=np.argmin(raw_data[mediadata.loc[t,'str']:mediadata.loc[t,'Apex']])+mediadata.loc[t,'str']
                    endpoint=np.argmin(raw_data[mediadata.loc[t,'Apex']:mediadata.loc[t,'end']+1])+mediadata.loc[t,'Apex']
                    ini_slp, _ = BasicPart.slope(raw_data, lpoint, rpoint)
                    str_thre = np.tan((np.arctan(deri_11[lpoint]) - np.arctan(ini_slp)) * liftoff)  # 将斜率化为角度计算后再化为斜率
                    end_thre = np.tan((np.arctan(-deri_11[rpoint]) + np.arctan(ini_slp)) * touchdown)
                    con=0
                    while True:
                        con+=1
                        if con>10:
                            print('v',lpoint,rpoint,endpoint,strpoint)
                        k,rpoint= funtemp_r(raw_data,lpoint,rpoint,endpoint+1,end_thre)
                        k,lpoint= funtemp_l(raw_data,lpoint,rpoint,strpoint-1,str_thre)
                        temp_r = -deri_11[rpoint] + k
                        if temp_r <= end_thre:
                            #print('v1')
                            if -deri_11[rpoint+1] + BasicPart.slope(raw_data, lpoint, rpoint+1)[0] <= end_thre:
                                #print('V2')
                                if -deri_11[rpoint+2] + BasicPart.slope(raw_data, lpoint, rpoint+2)[0] <= end_thre:
                                    break
                        if rpoint>= endpoint:
                            rpoint= endpoint
                            k,lpoint= funtemp_l(raw_data,lpoint,rpoint,strpoint-1,str_thre)
                            lpiont=max(lpoint,strpoint)
                            break
                        if lpoint <=strpoint:
                            lpoint =strpoint
                            k,rpoint= funtemp_r(raw_data,lpoint,rpoint,endpoint+1,end_thre)
                            rpiont=min(rpoint,endpoint)
                            break
                    mediadata.loc[t,'bslend']=mediadata.loc[t,'end']=rpoint
                    mediadata.loc[t,'bslstr']=mediadata.loc[t,'str']=lpoint
                    mediadata.loc[t, 'slp'],mediadata.loc[t, 'dft']= BasicPart.slope(raw_data, lpoint, rpoint, step)
        media.append(mediadata.copy())
        last=timelist[2*i+1]+1
    if last<=apex.iloc[-1]:
        mediadata=ApexData.loc[last:,:].copy()
        mediadata.index=mediadata.iloc[:,0]
        media.append(mediadata)


        '''
        if len(ApexData.loc[timelist[max(2*i-1,0)]+1:timelist[2*i],:]) !=0 :
            media.append(ApexData.loc[timelist[max(2*i-1,0)]+1:timelist[2*i],:])
        #print(media)
        mediadata=ApexData.loc[timelist[2*i]:timelist[2*i+1],:]
        cluster = pe.search_cluster(mediadata, 10, 11)
        for each in cluster:
            print(len(each))
            if len(each) > 1:
                for t in range(len(each)):
                    each.iloc[t, 3]=each.iloc[t, 1]
                    each.iloc[t, 4]=each.iloc[t, 2]
            media.append(each)
    
    '''
    ApexData=pd.concat(media,axis=0)
    ApexData.index=index
    return ApexData
            
def intersecting_detect(raw_data, FusePeak,liftoff,touchdown,step=1.,):
    """
    检测是否有基线穿透原数据，有的话重新分割峰并定位基线
    :param raw_data: array_like(原数据)
    :param ApexData: DataFrame(未经过基线穿透检测的DataFrame)
    :return ApexData: DataFrame(更正过的DataFrame)
    """
    #print(FusePeak)
    st = FusePeak.iloc[:, 1]
    en = FusePeak.iloc[:, 2]
    bs = FusePeak.iloc[:, 3]
    be = FusePeak.iloc[:, 4]
    apex=FusePeak.iloc[:,7]
    cluster = pe.search_cluster(FusePeak, 10, 11)
    
    # print(cluster)
    times = 0
    droplist=[]
    for each in cluster:
        #print(each)
        
        if len(each) == 1:
            #print()
            #print(ApexData.iloc[each[0],:],bs[each[0]], be[each[0]])
            bsl = BasicPart.baseline(raw_data, FusePeak.iloc[each[0], 3],FusePeak.iloc[each[0], 4])
            interval = raw_data[FusePeak.iloc[each[0], 3]:FusePeak.iloc[each[0], 4] + 1]
            m_data = interval - bsl
            if len(m_data[m_data<0])>0:
                #t+=1
                bs_new=apex[each[0]]
                be_new=apex[each[0]]
                if len(raw_data[bs[each[0]]:apex[each[0]]])>0:
                    bs_new=np.argmin(raw_data[bs[each[0]]:apex[each[0]]])+bs[each[0]]
                if len(raw_data[apex[each[0]]:be[each[0]]])>0: 
                    be_new=np.argmin(raw_data[apex[each[0]]:be[each[0]]+1])+apex[each[0]]
                FusePeak.iloc[each[0], 1]=FusePeak.iloc[each[0], 3]=bs_new
                FusePeak.iloc[each[0], 2]=FusePeak.iloc[each[0], 4]=be_new
                slp_l, dft_l = BasicPart.slope(raw_data, bs_new, be_new, step)
                FusePeak.iloc[each[0], 8] = slp_l
                FusePeak.iloc[each[0], 9] = dft_l
                bsl = BasicPart.baseline(raw_data, FusePeak.iloc[each[0], 3],FusePeak.iloc[each[0], 4])
                interval = raw_data[FusePeak.iloc[each[0], 3]:FusePeak.iloc[each[0], 4] + 1]
                m_data = interval - bsl
                if len(m_data[m_data<0])>0:
                    strpoint=FusePeak.iloc[each[0], 1]
                    endpoint=FusePeak.iloc[each[0], 2]
                    deri_11=np.diff(raw_data, 1)
                    lpoint=FusePeak.iloc[each[0], 5]
                    rpoint=FusePeak.iloc[each[0], 6]
                    ini_slp, _ = BasicPart.slope(raw_data, lpoint, rpoint)
                    str_thre = np.tan((np.arctan(deri_11[lpoint]) - np.arctan(ini_slp)) * liftoff)
                    end_thre = np.tan((np.arctan(-deri_11[rpoint]) + np.arctan(ini_slp)) * touchdown)
                    con=0
                    while True:
                        con+=1
                        k,rpoint= funtemp_r(raw_data,lpoint,rpoint,endpoint+1,end_thre)
                        k,lpoint= funtemp_l(raw_data,lpoint,rpoint,strpoint-1,str_thre)
                        temp_r = -deri_11[rpoint] + k
                        if con>10:
                            break
                        if temp_r <= end_thre:
                            #print('v1')
                            if -deri_11[rpoint+1] + BasicPart.slope(raw_data, lpoint, rpoint+1)[0] <= end_thre:
                                #print('V2')
                                if -deri_11[rpoint+2] + BasicPart.slope(raw_data, lpoint, rpoint+2)[0] <= end_thre:
                                    break
                        if rpoint>= endpoint:
                            rpoint= endpoint
                            k,lpoint= funtemp_l(raw_data,lpoint,rpoint,strpoint-1,str_thre)
                            lpiont=max(lpoint,strpoint)
                            break
                        if lpoint <=strpoint:
                            lpoint =strpoint
                            k,rpoint= funtemp_r(raw_data,lpoint,rpoint,endpoint+1,end_thre)
                            rpiont=min(rpoint,endpoint)
                            break
                    FusePeak.iloc[each[0], 1]=FusePeak.iloc[each[0], 3]=lpoint
                    FusePeak.iloc[each[0], 2]=FusePeak.iloc[each[0], 4]=rpoint
                #print(ApexData,ApexData.iloc[each[0], 2])
                #print(ApexData.iloc[each[0],:])
        
        if len(each) > 1:
            cluster_start = bs[each[0]]
            cluster_end = be[each[0]]
            #print(cluster_start)
            bsl = BasicPart.baseline(raw_data, cluster_start, cluster_end)
            interval = raw_data[cluster_start: cluster_end + 1]
            m_data = interval - bsl

            # intersect_point = []
            # intersect_delta = []

            # for i in range(each[0], each[-1]):
            #     valley_point = en[i] - cluster_start

            #     #if m_data[valley_point] <= 0:
            if len(m_data[m_data<0]) >5:
                valley_point=np.argmin(m_data)+cluster_start
                if valley_point<FusePeak.iloc[each[0],7]:
                    #穿透点在峰簇前
                    FusePeak.iloc[each[0]:each[-1]+1,3]=valley_point
                    FusePeak.iloc[each[0] ,1]=valley_point
                    slp_n, dft_n = BasicPart.slope(raw_data, FusePeak.iloc[each[0] ,1],cluster_end, step)
                    FusePeak.iloc[each[0]:each[-1]+1, 8] = slp_n
                    FusePeak.iloc[each[0]:each[-1]+1, 9] = dft_n
                elif valley_point>FusePeak.iloc[each[-1],7]:
                    FusePeak.iloc[each[0]:each[-1]+1,4]=valley_point
                    FusePeak.iloc[each[-1] ,2]=valley_point
                    slp_n, dft_n = BasicPart.slope(raw_data, cluster_start,FusePeak.iloc[each[-1] ,2], step)
                    FusePeak.iloc[each[0]:each[-1]+1, 8] = slp_n
                    FusePeak.iloc[each[0]:each[-1]+1, 9] = dft_n
                else:
                    i=int(FusePeak[FusePeak.iloc[:, 7]<=valley_point].iloc[-1,0])
                    times += 1
                    #intersect_point.append(i)
                    #intersect_delta.append(m_data[valley_point])
                    #if len(intersect_point) == 1:  # 只有一个基线穿透点
                    #i = intersect_point[0]
                    #print(i,FusePeak.iloc[i, 4], int(FusePeak.iloc[i, 2]))
                    # print(1, i)
                    FusePeak.iloc[each[0]: i + 1, 4] =  int(FusePeak.iloc[i, 2])  # 更新穿透点之前峰族子峰的基线终点

                    # if i == each[-1]:
                    #     FusePeak.iloc[i + 1, 3] =  int(FusePeak.iloc[i , 2]) # 更新穿透点之后峰族子峰的基线起点
                    # else:
                    FusePeak.iloc[i + 1: each[-1]+1, 3] =  int(FusePeak.iloc[i ,2])
                    # print(en[i])
                    #print(i,FusePeak.iloc[i, 4], int(FusePeak.iloc[i, 2]))
                    FusePeak.iloc[i + 1, 10] = 'B'  # 更新穿透点后一处子峰的峰起始类型
                    FusePeak.iloc[i, 11] = 'B'  # 更新穿透点前一处子峰的峰结束类型
                    slp_l, dft_l = BasicPart.slope(raw_data, cluster_start,FusePeak.iloc[i, 2], step)
                    slp_r, dft_r = BasicPart.slope(raw_data, FusePeak.iloc[i+1, 1], cluster_end, step)
                    FusePeak.iloc[each[0]: i + 1, 8] = slp_l  # 更新穿透点之前的峰族子峰的基线斜率和漂移
                    FusePeak.iloc[each[0]: i + 1, 9] = dft_l
                    FusePeak.iloc[i + 1: each[-1]+1, 8] = slp_r  # 更新穿透点之后的峰族子峰的基线斜率和漂移
                    FusePeak.iloc[i + 1: each[-1]+1, 9] = dft_r
                
                # if len(intersect_point) > 1:  # 多个基线穿透点
                #     #print('t')
                #     # print('n', intersect_point)
                #     temp = np.argmin(intersect_delta)
                #     i = intersect_point[temp]
                #     #print(i,FusePeak.iloc[i, 4], int(FusePeak.iloc[i, 2]))
                #     # print('low', i)
                #     FusePeak.iloc[each[0]: i + 1, 4] = int(FusePeak.iloc[i , 2])
                #     #ApexData.iloc[i + 1: each[-1], 3] = en[i]
                #     if i + 1 == each[-1]:
                #         FusePeak.iloc[i + 1, 3] = int(FusePeak.iloc[i , 2]) # 更新穿透点之后峰族子峰的基线起点
                #     else:
                #         FusePeak.iloc[i + 1: each[-1]+1, 3] = int(FusePeak.iloc[i , 2])
                #     #print(each)
                #     #print(i,FusePeak.iloc[i, 4], int(FusePeak.iloc[i, 2]))
                #     FusePeak.iloc[i + 1, 10] = 'B'
                #     FusePeak.iloc[i, 11] = 'B'
                #     slp_l, dft_l = BasicPart.slope(raw_data, cluster_start, en[i], step)
                #     slp_r, dft_r = BasicPart.slope(raw_data, en[i], cluster_end, step)

                #     FusePeak.iloc[each[0]: i + 1, 8] = slp_l
                #     FusePeak.iloc[each[0]: i + 1, 9] = dft_l
                #     FusePeak.iloc[i + 1: each[-1]+1, 8] = slp_r
                #     FusePeak.iloc[i + 1: each[-1]+1, 9] = dft_r
    # #print(ApexData)
    # if len(droplist)!=0:
    #     #print(len(FusePeak))
    #     FusePeak=FusePeak.drop(droplist)
    #     #print(len(FusePeak))
    #     FusePeak.index=FusePeak['name']=np.arange(FusePeak.shape[0])
    return FusePeak, times


def reten_time(a, d2, name, s, e, infl, infr, apex, slp, dft, step=1.):
    """
    函数：用基线对原数据进行修正并得到真正的Apex顶点，对其进行抛物线拟合算出真实的保留时间
    :param a: array_like(用于计算保留时间的原数据)
    :param name: int(峰编号，用于标记，函数内不做处理)
    :param s: int(峰起点)
    :param e: int(峰终点)
    :param apex: int(二阶导局部最大值点即apex顶点)
    :param slp: float(基线斜率)
    :param dft: float(基线漂移)
    :param step: float, optional(数据点时间间隔，默认值为1)
    :return name : int(峰名称)
    :return rt : float(计算得到的保留时间)
    :return ph : float(计算得到的峰高[排除基线影响])
    :return area : float(积分面积)
    """
    # 计算峰起落点间距
    length = e - s
    # 时间数据
    t_data = np.asarray(range(s, e + 1)) * step

    # 基线以及基线修正后的原数据(index为[0: e - s + 1])
    bsl_data = t_data * slp + dft
    m_data = np.asarray(a[s: e + 1]) - bsl_data

    # 峰区域上apex点的index
    apex=apex
    pre_peak = apex - s
    #area = integrate.simps(m_data, t_data)
    #gussian_fit_data=m_data[m_data>=0.5*m_data[pre_peak]]
    '''
    if len(gussian_fit_data)>10:
        fit_sigma=len(gussian_fit_data)*step/(2.354)
        fit_miu=pre_peak*step
        fit_A=m_data[pre_peak]
        y = gussian_fit_data#m_datanp.asarray(m_data[s: e])
        x = np.arange(len(gussian_fit_data))#np.asarray(range(s, e+1)) * step#-s*step
        #print(s*step,x,'\n',fit_miu)
        p= BasicPart.gaussian_fit(x,y,[fit_A,fit_miu,fit_sigma])
        rt = p[1]+s*step
        ph = p[0]
        z=BasicPart.gaussian(t_data-s*step,p[0],p[1],p[2])#*m_data[pre_peak]
        #print(p,m_data,'\n',z)
        #area = integrate.simps(z, t_data)
        return name, rt, ph, area, None
    '''
    # 计算积分面积
    area = integrate.simps(m_data, t_data)

    # 计算拐点宽度
    inf_len = infr - infl

    # Apex点位于边界时，直接采用apex顶点
    if apex == s or apex == e:
        rt = apex * step
        ph = a[apex] - bsl_data[pre_peak]

        return name, rt, ph, area, 'margin'

    # 拐点宽小于四时，对原数据三点拟合
    if inf_len <= 4:
        y = np.asarray(m_data[pre_peak - 1: pre_peak + 2])
        x = np.asarray(range(apex - 1, apex + 2)) * step
        z = np.polyfit(x, y, 2)

        # 计算保留时间、峰高、面积
        rt = (-z[1] / (2. * z[0]))  # 顶点处一阶导数为0
        if s > (rt / step) or e < (rt / step):
            rt = apex * step
            ph = a[apex] - bsl_data[pre_peak]
            return name, rt, ph, area, '3, out of range'
        else:
            ph = z[0] * rt * rt + z[1] * rt + z[2]
            return name, rt, ph, area, '3'

    # 拐点宽大于四时，对原数据五点拟合
    if inf_len > 4:
        if apex - 2 < s or apex + 3 > e:  # apex点两侧不足5个点时
            y = np.asarray(m_data[pre_peak - 1: pre_peak + 2])
            x = np.asarray(range(apex - 1, apex + 2)) * step
            z = np.polyfit(x, y, 2)
            # 计算保留时间、峰高、面积
            rt = (-z[1] / (2. * z[0]))  # 顶点处一阶导数为0
            if s > (rt / step) or e < (rt / step):
                rt = apex * step
                ph = a[apex] - bsl_data[pre_peak]
                return name, rt, ph, area, '3n, out of range'
            else:
                ph = z[0] * rt * rt + z[1] * rt + z[2]
                return name, rt, ph, area, '3'
        else:  # 正常情况
            y = np.asarray(m_data[pre_peak - 2: pre_peak + 3])
            x = np.asarray(range(apex - 2, apex + 3)) * step
            z = np.polyfit(x, y, 2)
            # 计算保留时间、峰高、面积
            rt = (-z[1] / (2. * z[0]))  # 顶点处一阶导数为0
            if s > (rt / step) or e < (rt / step):
                peak = np.argmax(m_data) + s
                rt = peak * step
                ph = a[peak] - bsl_data[peak - s]
                return name, rt, ph, area, 'out of range'
            else:
                ph = z[0] * rt * rt + z[1] * rt + z[2]
                return name, rt, ph, area, None


def area_ratio(area):
    """
    计算面积百分比
    :param area: array_like(积分面积)
    :return ratio: array_like(面积百分比)
    """
    # print('area', area)
    total = np.sum(area)
    ratio = []
    for i in range(len(area)):
        if area[i] is None:  # 峰不足5个点时认为峰%面积为0
            ratio.append(0)
        else:
            #temp = np.around((area[i] / total) * 100, 2)  # 返回峰%面积，保留两位小数
            temp = (area[i] / total) * 100
            ratio.append(temp)

    return ratio


def baseline(data, x1, x2):
    """
    生成基线数据
    :param data: array_like(原数据)
    :param x1: int(基线起点)
    :param x2: int(基线终点)
    :return bsl: array_like()
    """
    y1 = data[x1]
    y2 = data[x2]
    n = x2 - x1
    delta = y2 - y1
    step = delta / n
    bsl = np.arange(y1, y2, step)
    bsl = np.append(bsl, y2)
    return bsl

def final_step(ApexData, d2, raw_data, raw_step=1.):
    """
    在ApexData的基础上计算出保留时间，峰高，峰面积等信息，合成并返回包含所有信息的DataFrame
    :param ApexData: DataFrame
    :param raw_data: array_like(原始数据)
    :param raw_step: float, optional(数据点时间间隔)
    :return FullData: DataFrame
    """
    '''
    media=[]
    length = len(ApexData)
    #NewData=ApexData.copy()
    for i in range(1,length-1):
        #media.append(ApexData.iloc[i,:])
        if (ApexData.iloc[i+1, 1]-ApexData.iloc[i, 2]<3) and (ApexData.iloc[i, 2]-ApexData.iloc[i, 7]<=5):
            ApexData.iloc[i+1, 1]=ApexData.iloc[i, 1]
            ApexData.iloc[i+1, 3]=ApexData.iloc[i, 3]
            ApexData.iloc[i+1, 11] =ApexData.iloc[i, 11]
            ApexData.iloc[i+1,7]=np.argmax(raw_data[ApexData.iloc[i+1, 1]:ApexData.iloc[i+1, 2]])+ApexData.iloc[i+1, 1]
            #NewData=NewData.drop(i)
            media.append(i)
        elif (ApexData.iloc[i, 1]-ApexData.iloc[i-1, 2]<3) and (ApexData.iloc[i, 7]-ApexData.iloc[i, 1]<=5):
            #print(ApexData.iloc[i-1, 2],ApexData.iloc[i, 2])
            ApexData.iloc[i-1, 2]=ApexData.iloc[i, 2]
            ApexData.iloc[i-1, 4]=ApexData.iloc[i, 4]
            ApexData.iloc[i-1, 10] =ApexData.iloc[i, 10]
            ApexData.iloc[i-1,7]=np.argmax(raw_data[ApexData.iloc[i-1, 1]:ApexData.iloc[i-1, 2]])+ApexData.iloc[i-1, 1]
            #media[-1]=ApexData.iloc[i-1,:]
            media.append(i)
            #print(len(media),i)
            #print(len(media),i)
        
    #ApexData=pd.concat(media,axis=1).T
    for i in media:
        ApexData=ApexData.drop(i)
    ApexData.index=range(len(ApexData))
    ApexData['name']=range(len(ApexData))
    #ApexData=NewData.copy()
    #print(ApexData)
    #'''
    #print(ApexData)
    cluster = pe.search_cluster(ApexData, 10, 11)
    length = len(ApexData)
    s = ApexData.iloc[:, 1]
    e = ApexData.iloc[:, 2]
    bs = ApexData.iloc[:, 3]
    be = ApexData.iloc[:, 4]
    infl = ApexData.iloc[:, 5]
    infr = ApexData.iloc[:, 6]
    apex = ApexData.iloc[:, 7]
    slp = ApexData.iloc[:, 8]
    dft = ApexData.iloc[:, 9]
    detail = []
    for i in range(length):
        detail.append(reten_time(raw_data, d2, i, s[i], e[i], infl[i], infr[i], apex[i], slp[i], dft[i], raw_step))

    Detail = pd.DataFrame(detail, columns=['name', 'rete', 'p_high', 'area', 'error'])
    area = Detail.iloc[:, 3]
    Detail['a_ratio'] = None
    a_ratio = area_ratio(area)
    Detail.iloc[:, 5] = a_ratio
    FullData = pd.merge(ApexData, Detail, on='name')
    return FullData


def output(FullData, raw_data, step=1.):
    """
    整理数据，排除过小峰，整理表头，输出最终结果DataFrame
    :param FullData: DataFrame(输入数据框)
    :param step: float, optional(数据时间间隔，缺省值为1.0)
    :return Output: DataFrame(输出结果)
    """
    name = FullData.iloc[:, 0]
    s = FullData.iloc[:, 1]
    e = FullData.iloc[:, 2]
    bsl = FullData.iloc[:, 3]
    bsr = FullData.iloc[:, 4]
    slp = FullData.iloc[:, 8]
    dft = FullData.iloc[:, 9]
    type1 = FullData.iloc[:, 10]
    type2 = FullData.iloc[:, 11]
    bsh=FullData.iloc[:, 12]
    beh = FullData.iloc[:, 13]
    rt = FullData.iloc[:, 14]
    ph = FullData.iloc[:, 15]
    area = FullData.iloc[:, 16]

    ar = FullData.iloc[:, 18]

    num = len(FullData)

    # 整理输出数据
    out = []
    for i in range(num):
        # 与empower统一单位
        #ren = np.around(rt[i], 3)
        ren=rt[i]
        #m_area = np.round(area[i] * 60000000)
        m_area = area[i]
        #m_ph = np.around(ph[i] * 1000000)  # 保留三位小数
        m_ph = ph[i]
        start_time =s[i]
        end_time = e[i]
        bs_time = bsl[i]
        be_time = bsr[i]
        # 合并积分类型
        p_type = type1[i] + type2[i]
        # 科学计数法
        """
        slp_temp, dft_temp = BasicPart.slope(raw_data, bsl[i], bsr[i], step)
        slope = format(slp_temp, '.3e')
        drift = format(dft_temp, '.3e')
        """
        slope = format(slp[i], '.3e')
        drift = format(dft[i], '.3e')
        points = e[i] - s[i]+1
        bs_h=bsh[i]
        be_h=beh[i]
        width =(points-1)*step

        out.append([name[i], ren, m_area, ar[i], m_ph, p_type, start_time, end_time, bs_time, be_time,
                    slope, drift, points,width,bs_h,be_h])

    OutPut = pd.DataFrame(out, columns=['Peak Number', 'Retention Time', 'Area', '% Area', 'Height', 'Int Type', 'Start Point', 'End Point',
                                        'Baseline Start', 'Baseline End', 'Slope', 'float','Points Across Peak', 'Width','Baseline Start Height','Baseline End Height'])
    #OutPut = pd.DataFrame(out, columns=['峰编号', '保留时间', '面积', '%面积', '峰高', '积分类型', '开始时间', '结束时间',
    #                                    '基线开始', '基线结束', '斜率', '漂移', '峰内点数'])
    '''
    final=[]
    for i in range(1,OutPut.shape[0]-1):
        if OutPut['End Time'][i]==OutPut['Start Time'][i+1] and OutPut['End Time'][i]-OutPut['Retention Time'][i]<0.01 :
            OutPut['Start Time'][i+1]=OutPut['Start Time'][i]
            OutPut['%Area'][i+1]+=OutPut['%Area'][i]
            OutPut['Area'][i+1]+=OutPut['Area'][i]
            OutPut['Area'][i+1]+=OutPut['pointnumber'][i]
            '''
    return OutPut
