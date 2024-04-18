import numpy as np
from scipy import optimize as op
from scipy import signal

def infl(input_data, x, d_2=None):
    """
    函数：求拐点及拐点间隔
    * 仅考虑正峰情况
    :param a: array_like(输入需要求拐点的数组)
    :param x: int(Apex峰顶点的index)
    :param d_2: array_like, optional(二阶导数，缺省时自动计算)
    :return i : int(左侧拐点)
    :return j : int(右侧拐点)
    :return j - i: int(拐点距离)
    """
    n = len(input_data)
    if d_2 is None:
        #求二阶导
        d_2 = -np.diff(input_data, 2)
    rp=x-1
    lp=x-1
    for i in range(x - 1, 1, -1):
        if i != 1:#在数据结尾处强制终止
            if d_2[i] == 0 or d_2[i] * d_2[min(i+1,n-2)] < 0:#找左端点
                #加入判断 
                if d_2[max(i - 1,1)] * d_2[min(i + 2,x - 1)] > 0:
                    continue
                lp = i
                break
            #if d_2[max(i - 2,1)] * d_2[min(i + 3,x - 1)] > 0:
                #continue
                #剩下全部缩进
            #判断点二阶导是否为0
    for j in range(x, n - 2, 1):
        if j == n - 3:#在数据结尾处强制终止
            return lp + 1, j + 1, j - lp
        if d_2[j] == 0 or d_2[max(j - 1,1)] * d_2[j] < 0:#右端点
            if d_2[max(j - 1,x)] * d_2[min(j + 2,n - 2)] > 0:
                continue
            rp=j
            break
            #if d_2[max(j - 5,x)] * d_2[min(j + 5,n - 2)] > 0:
                #continue
            #d_2[j - 2] * d_2[j + 1]>0:
            #continue
            #剩下全部缩进
            #注意n的取值范围
    #print(lp + 1, rp + 1, rp-lp)
    return lp + 1, rp + 1, rp-lp
    #return x, x, 0


def threshold(raw_data, right_num, lift_num):
    """
    函数：阈值计算
    * 阈值为噪声区域标准差四倍
    :param a: array_like(需要计算阈值的数组)
    :return thre: float(通过输入噪声区域计算得到的阈值)
    """
    thredata=np.array(raw_data[int(right_num): int(lift_num)])
    thredata=signal.resample(thredata,106)
    thre = np.std(thredata) * 4.
    return thre


def six_sigma(input_data):
    """
    函数：6σ方法计算噪音
    :param a: array_like(需要计算噪声的数组)
    :return ns: float(给定区域上输入数组的噪声)
    """
    l = int(input('请输入计算噪声所需区域的起点：'))
    r = int(input('请输入计算噪声所需区域的终点：'))
    temp_x = np.asarray(range(l, r))
    temp_y = np.asarray(input_data[l: r])
    num = r - l
    k, b = np.polyfit(temp_x, temp_y, 1)#拟合线性
    y = k * temp_x + b
    delta = (y - temp_y) ** 2.#偏差计算
    s = sum(delta)
    ns = 6. * np.sqrt(s / num)
    return ns


def peakwidth(input_data, x=None):
    """
    函数：峰宽计算
    :param a: array_like(需要计算阈值的数组)
    :param x: int(需要计算峰宽的峰顶点的index,为'None'时计算自动峰宽)
    :return pw_index: int(峰宽值取整)
    :return pw: float(峰宽值，未乘时间间隔)
    """
    if x is None:
        deri_2 = -np.diff(input_data, 2)#二阶导
        x = np.argmax(deri_2)#最大值索引
    il, ir, temp_index = infl(input_data, x)#找拐点
    # print('infl width =', temp_index)
    # print('l =', il, 'r =', ir)
    pw = temp_index * 4.89549 / 2.
    pw_index = int(np.round(pw))
    return pw_index, pw


def filter_parameter(input_data):
    """
    函数：滤波系数计算(必须为奇数(为偶数则+1)且不小于7)
    :param a: array_like(需要计算滤波系数的数组)
    :return flt_par: int(计算得到的滤波系数)
    """
    auto_pw, _ = peakwidth(input_data)
    flt_par = int(np.round(auto_pw / 1.5))
    if flt_par % 2 == 0:
        flt_par += 1
    #if flt_par <= 3:
        #flt_par = 3
    if flt_par <= 7:
        flt_par = 7

    # print('filter parameter =', flt_par)
    return flt_par


def slope(input_data, strponit, endpoint, step=1.):
    """
    函数：计算某曲线上给定两点连线的斜率以及漂移 √
    :param a: array_like(原数据数组)
    :param i: int(计算斜率的起点)
    :param j: int(计算斜率的终点)
    :param step: float, optional(数据点间的时间间隔，默认值为1)
    :return slp: float(计算得到的线段斜率)
    :return dft: float(计算得到的线段漂移[线段所在直线在y轴上的截距])
    """
    if strponit == endpoint:
        return 0, 0

    elif strponit < endpoint:
        x = [strponit * step, endpoint * step]
        y = [input_data[strponit], input_data[endpoint]]

    else:
        x = [endpoint * step, strponit * step]
        y = [input_data[endpoint], input_data[strponit]]
    slp = (y[1] - y[0]) / (x[1] - x[0])
    dft = y[0] - slp * x[0]

    return slp, dft


def baseline(data, bs, be):
    """
    生成基线数据
    :param data: array_like(原数据)
    :param x1: int(基线起点)
    :param x2: int(基线终点)
    :return bsl: array_like()
    """
    y1 = data[bs]
    y2 = data[be]
    num = be - bs + 1
    """
    n = x2 - x1
    delta = y2 - y1
    step = delta / n
    bsl = np.arange(y1, y2, step)

    if bsl[-1] != y2:
        #print('-1', bsl[-1], 'y2', y2, 'judge', bsl[-1] != y2)
        bsl = np.append(bsl, y2)
    #print('x1', x1, 'x2', x2, 'n', n, 'contain', len(bsl))
    #print(bsl)
    return bsl
    """
    bsl = np.linspace(y1, y2, num)

    return bsl


# def data_input(input_path, skiprow=None):
#     """
#     函数：读取数据 √
#     :param a_path: string(存放数据的文件的路径)
#     :param a: int, optional(读取有标题的数据文件时用于跳过标题行，a为需要跳过的行数，缺省情况下不跳过)
#     :return RawData: DataFrame(直接从文件中读取的原始数据，包含时间和强度信息)
#     :return raw_data: array_like(原始数据中的强度列)
#     :return raw_step: float(通过原始时间数据算出的时间间隔)
#     """
#     if skiprow is None:
#         RawData = pd.read_csv(input_path, sep='\s+', header=None,dtype='float64')
#         RawData.columns=['time', 'intensity']

#     else:
#         temp = list(range(skiprow))
#         RawData = pd.read_csv(input_path, sep='\s+', header=[0], skiprows=temp,dtype='float64')
#         RawData.columns = ['time', 'intensity']
#     print(RawData)
#     raw_data = np.asarray(RawData.iloc[:, 1])
#     raw_step = RawData.iloc[1, 0] - RawData.iloc[0, 0]
#     #raw_data=np.hstack((raw_data[0:7],raw_data))
#     #Hz=1/raw_step
#     b, skiprow=signal.butter(8, 0.9, 'lowpass')  
#     raw_data=signal.filtfilt(b, skiprow, raw_data)

#     return RawData, raw_data, raw_step

def data_input_ob(data):
    """
    函数：读取数据 √
    :param a_path: string(存放数据的文件的路径)
    :param a: int, optional(读取有标题的数据文件时用于跳过标题行，a为需要跳过的行数，缺省情况下不跳过)
    :return RawData: DataFrame(直接从文件中读取的原始数据，包含时间和强度信息)
    :return raw_data: array_like(原始数据中的强度列)
    :return raw_step: float(通过原始时间数据算出的时间间隔)
    """
    RawData = data
    RawData.columns=['time', 'intensity']
    #print(RawData)
    raw_data = np.asarray(RawData.iloc[:, 1])
    raw_step = (RawData.iloc[-1, 0] - RawData.iloc[0, 0])/(len(RawData)-1)
    #raw_data=np.hstack((raw_data[0:1],raw_data))
    #Hz=1/raw_step
    b, skiprow=signal.butter(8, 0.9, 'lowpass')  
    raw_data=signal.filtfilt(b, skiprow, raw_data)
    return RawData, raw_data, raw_step


def gaussian(x,a,b,c):
    """
    高斯函数公式(与scipy.signal中的不同，用于高斯切削等步骤)
    :param x: array_like(自变量数组)
    :param param: array_like(参数数组，包含三个参数，param[0]为H0，param[1]为x0，param[2]为σ)
    :return gau: array_like(根据数组和参数得到的高斯函数数据点)
    """
    gau = a * np.exp(-np.power(x - b, 2.) / (2 * np.power(c, 2.)))
    return gau


def gaussian_fit(x,y,p0):
    """
    高斯拟合
    :param a: array_like(进行拟合所用的数据)
    :param param: array_like(参数数组，包含三个参数，param[0]为H0，param[1]为x0，param[2]为σ)
    :return g_fit: array_like(拟合得到的数据)
    """
    #num = len(y)
    #x = np.asarray(range(num))
    #y = y
    p, g_fit = op.curve_fit(gaussian, x, y, p0)
    return p ,g_fit

# def setBslExpr(curvetype, coefs):
#     """
#     函数：基线/切削曲线参数转为dict形式进行描述
#     :param functype: string(待转换曲线类型, "linear"/"exp"/"gaussian")
#     :param coefs: float[](参数组)
#     :return res: dict(组装好的曲线参数形式)
#     """
#     if functype == "linear":
#         # 直线 y = k * x + b
#         res = {"type": "linear", "k": coefs[0], "b": coefs[1]}
#     elif functype == "exp":
#         # 指数曲线 y = m * exp(n * x)
#         res = {"type": "exp", "m": coefs[0], "n": coefs[1]}
#     elif functype == "gaussian":
#         # 高斯曲线 y = n(mu, sigma)
#         res = {"type": "gaussian", "mu": coefs[0], "sigma": coefs[1], "b": coefs[2]}
#     else:
#         res = None
#     return res

# def fitCurve(curvetype, rawData, timeData, startPoint, endPoint):
#     """
#     直接对指定区间、指定类型的曲线进行参数拟合，返回组装好的dict形式曲线参数
#     """
#     Y = rawData[startPoint:endPoint+1]
#     X = timeData[startPoint:endPoint+1]
#     if curvetype == "linear":
#         coefs =  
#     elif curvetype == "exp":
#         coefs =
#     elif curvetype == "gaussian":
#         coefs = 
#     else:
#         return None
#     return setBslExpr(curvetype, coefs)