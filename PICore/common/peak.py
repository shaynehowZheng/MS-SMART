from scipy import integrate
import numpy as np

"""
画图建议：
所有点都是索引值，需要在时间轴上找出来对应时间点
基线参数是根据时间轴计算的，可直接画在图上
"""
class Baseline():
    """
    输出格式：
    :param: startPoint: int 基线起点索引
    :param: endPoint: int 基线终点索引
    :param: curveType: string 基线曲线类型：linear/exp/gaussian
    :param: curvePoints: []float 基线曲线点集合，直线为首尾两个端点，曲线为全部点的高度集合
    """
    def __init__(self, startpoint = None, endpoint = None, rawData = None, inteval = None,bslType =None,startY = None,endY = None):
        if bslType =='auto':
            try:
                a = rawData[endpoint]
                b = rawData[startpoint]
                slope = (rawData[endpoint] - rawData[startpoint]) / ((endpoint - startpoint)*inteval)
                drift = rawData[startpoint]-slope * startpoint*inteval
            except:
                slope = None
                drift = None
        else:
            #正向水平的话是以第一个点
            if bslType =='Horizon':
                if startY !=None and endY !=None:
                    drift = endY
                    slope = 0
                else:
                    try:
                        slope = 0
                        drift = rawData[startpoint]
                    except:
                        slope = None
                        drift = None
            else:
                # 反向水平在则是最后一个点
                if startY !=None and endY !=None:
                    drift = endY
                    slope = 0
                else:
                    try:
                        slope = 0
                        drift = rawData[endpoint]

                    except:
                        slope = None
                        drift = None
        self.startPoint = startpoint
        self.endPoint = endpoint
        self.slope = slope
        self.drift = drift
        self.bslType = bslType
        self.endY = endY
        self.startY = startY


    def generatebslData(self,rawData,inteval):
        s, e = self.startPoint, self.endPoint
        k, b = self.slope, self.drift
        # bsl_data = np.linspace(rawData[self.startPoint], rawData[self.endPoint], e - s + 1)  #当峰是融合峰时就会出错 因为融合峰的终点不是基线的终点
        bsl_data = np.linspace(k*s*inteval + b, k*e*inteval + b, e-s+1)
        return bsl_data

    def exportDict(self,inteval):
        return {"startIdx": self.startPoint,
                "endIdx": self.endPoint,
                "slope": self.slope,
                "drift":self.drift/((self.startPoint+1)*inteval/3600),
                "bslType":self.bslType,
                "endY":self.endY,
                "startY":self.startY
                }
class Peak():
    """
    输出格式：
    """
    def __init__(self, clusterID = None, peakID= None, peakType= None, 
                    startPoint= None, apexPoint= None, endPoint= None, 
                    retentionTime = None, area = None, height = None,
                    percentArea = None, percentHeight = None,baseline = None,startBLY = None,endBLY = None):
        self.ClusterID = clusterID
        self.PeakID = peakID
        self.ApexPoint = apexPoint
        self.RetentionTime = retentionTime
        self.Area = area
        self.PercentArea = percentArea
        self.PercentHeight = percentHeight
        self.Height = height
        self.PeakType = peakType
        self.Baseline = baseline
        self.startPoint = startPoint
        self.endPoint = endPoint
        self.endBLY = endBLY
        self.startBLY = startBLY




    def findBLY(self,bsl_data,startBslIdx):
        endBLY = bsl_data[self.endPoint-startBslIdx]
        startBY = bsl_data[self.startPoint-startBslIdx]
        self.startBLY = startBY
        self.endBLY = endBLY


    # 更新基线
    def updateBsl(self, baseline):
        self.Baseline = baseline
    # 积分、保留时间
    def peakCalculation(self, rawData, timeData,step):
        """

        :param rawData:
        :param timeData:
        :param step: 时间间隔
        :return:
        """
        # 峰面积积分
        s = self.startPoint
        e = self.endPoint
        bsln = self.Baseline
        k, b = bsln.slope, bsln.drift
        bsl_data = np.linspace(k*s*step + b, k*e*step + b, e-s+1)
        # bsl_data = np.linspace(rawData[self.startPoint], rawData[self.endPoint], e - s + 1)

        t_data = np.asarray(timeData[s: e + 1])
        m_data = np.asarray(rawData[s: e + 1]) - bsl_data   
        self.Area = abs(integrate.simps(m_data, t_data))
        # 保留时间
        apex = self.ApexPoint
        #直接根据坐标得到的结果
        # rt = timeData[apex]
        # ph = rawData[apex]
        # self.RetentionTime = rt
        # self.Height = ph

        # 原始的计算
        pre_peak = apex - s
        # _, _, inf_len = infl(rawData[s: e + 1],apex) #求拐点
        # 由于拐点计算有问题，此处直接用峰的x轴宽度1/2作为峰宽使用
        inf_len = (e - s) / 2
        # 以后需要重构此处计算
        # 顶点位于边界时，直接采用apex顶点
        if apex == s or apex == e:
            rt = apex * step
            ph = rawData[apex] - bsl_data[pre_peak]
        # 拐点宽小于四时，对原数据三点拟合
        elif inf_len <= 4:
            y = np.asarray(m_data[pre_peak - 1: pre_peak + 2])
            x = np.asarray(range(apex - 1, apex + 2)) * step
            z = np.polyfit(x, y, 2)
            # 计算保留时间、峰高
            rt = (-z[1] / (2. * z[0]))  # 顶点处一阶导数为0
            if s > (rt / step) or e < (rt / step):
                rt = apex * step
                ph = rawData[apex] - bsl_data[pre_peak]
            else:
                ph = z[0] * rt * rt + z[1] * rt + z[2]
        # 拐点宽大于四时，对原数据五点拟合
        else: 
            if apex - 2 < s or apex + 3 > e:  # apex点两侧不足5个点时
                y = np.asarray(m_data[pre_peak - 1: pre_peak + 2])
                x = np.asarray(range(apex - 1, apex + 2)) * step
                z = np.polyfit(x, y, 2)
                # 计算保留时间、峰高
                rt = (-z[1] / (2. * z[0]))  # 顶点处一阶导数为0
                if s > (rt / step) or e < (rt / step):
                    rt = apex * step
                    ph = rawData[apex] - bsl_data[pre_peak]
                else:
                    ph = z[0] * rt * rt + z[1] * rt + z[2]
            else:  # 正常情况
                y = np.asarray(m_data[pre_peak - 2: pre_peak + 3])
                x = np.asarray(range(apex - 2, apex + 3)) * step
                z = np.polyfit(x, y, 2)
                # 计算保留时间、峰高
                rt = (-z[1] / (2. * z[0]))  # 顶点处一阶导数为0
                if s > (rt / step) or e < (rt / step):
                    peak = np.argmax(m_data) + s
                    rt = peak * step
                    ph = rawData[peak] - bsl_data[peak - s]
                else:
                    ph = z[0] * rt * rt + z[1] * rt + z[2]
        self.RetentionTime = rt
        self.Height = ph
    def exportDict(self):
        # 转为字典
        return {"PeakID":self.PeakID,
                            "ApexIdx":self.ApexPoint,
                            "RetentionTime":self.RetentionTime,
                            "Area":self.Area,
                            "PercentArea":self.PercentArea,
                            "PercentHeight": self.PercentHeight,
                            "Height":self.Height,
                            "PDType":self.PeakType,
                            "StartIdx":self.startPoint,
                            "EndIdx":self.endPoint,
                            "endBLY":self.endBLY,
                            "startBLY":self.startBLY
                }


class PeakCluster():
    def __init__(self, clusterID = None):
        self.ClusterID = clusterID
        self.Peaks = []
        self.Baseline = None
    def setBaseline(self, baseline):
        self.Baseline = baseline
    def updateAllPeaks(self):
        pass

class integralResult():
    def __init__(self, clusters = []):
        self.Clusters = clusters
    def ExportResDict(self,inteval):
        if len(self.Clusters) == 0:
            return None
        else:
            res = []
            for cluster in self.Clusters:
                curr = {"ClusterID":cluster.ClusterID, 
                        "Peaks":[],
                        "Baseline":cluster.Baseline.exportDict(inteval)}
                for p in cluster.Peaks:
                    curr["Peaks"].append(p.exportDict())
                res.append(curr)
        return res
        #json.dumps(res)
    def UpdatePercentArea(self):
        sumarea = 0
        sumheight = 0
        for cluster in self.Clusters:
            for peak in cluster.Peaks:
                sumarea += peak.Area
                sumheight += peak.Height
        for cluster in self.Clusters:
            for peak in cluster.Peaks:
                peak.PercentArea = peak.Area/sumarea*100
                peak.PercentHeight = peak.Height / sumheight * 100
    # def ImportJsonStr(self,impstr):
    #     pass