import numpy as np


# ###最小面积的过滤事件 这是过滤完 后基线不变的版本
def Areafilter(pre_result, startTime, endTime, paramValue, inteval):
    startinx = int(startTime // inteval)
    endinx = int(endTime // inteval)

    for index in range(len(pre_result) - 1, -1, -1):

        for i in range(len(pre_result[index]["Peaks"]) - 1, -1, -1):
            if pre_result[index]["Peaks"][i]["StartIdx"] >= startinx and pre_result[index]["Peaks"][i][
                "EndIdx"] <= endinx:
                if pre_result[index]["Peaks"][i]["Area"] < paramValue:
                    del pre_result[index]["Peaks"][i]

        if len(pre_result[index]["Peaks"]) == 0:
            del pre_result[index]

    return pre_result




def Heightfilter(pre_result, startTime, endTime, flag, paramValue, inteval):
    startinx = int(startTime // inteval)
    endinx = int(endTime // inteval)

    if flag == "SetMinHeight":

        for index in range(len(pre_result) - 1, -1, -1):

            for i in range(len(pre_result[index]["Peaks"]) - 1, -1, -1):
                if pre_result[index]["Peaks"][i]["StartIdx"] >= startinx and pre_result[index]["Peaks"][i][
                    "EndIdx"] <= endinx:
                    if pre_result[index]["Peaks"][i]["Height"] < paramValue:
                        del pre_result[index]["Peaks"][i]

            if len(pre_result[index]["Peaks"]) == 0:
                del pre_result[index]
    else:

        for index in range(len(pre_result) - 1, -1, -1):

            for i in range(len(pre_result[index]["Peaks"]) - 1, -1, -1):
                if pre_result[index]["Peaks"][i]["StartIdx"] >= startinx and pre_result[index]["Peaks"][i][
                    "EndIdx"] <= endinx:
                    if pre_result[index]["Peaks"][i]["Height"] > paramValue:
                        del pre_result[index]["Peaks"][i]

            if len(pre_result[index]["Peaks"]) == 0:
                del pre_result[index]

    return pre_result


def Widthfilter(pre_result, startTime, endTime, paramValue, inteval):
    startinx = int(startTime // inteval)
    endinx = int(endTime // inteval)

    for index in range(len(pre_result) - 1, -1, -1):

        for i in range(len(pre_result[index]["Peaks"]) - 1, -1, -1):
            if pre_result[index]["Peaks"][i]["StartIdx"] >= startinx and pre_result[index]["Peaks"][i][
                "EndIdx"] <= endinx:
                if pre_result[index]["Peaks"][i]["EndIdx"] - pre_result[index]["Peaks"][i]["StartIdx"] > int(
                        paramValue // inteval):
                    del pre_result[index]["Peaks"][i]

        if len(pre_result[index]["Peaks"]) == 0:
            del pre_result[index]

    return pre_result


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


def VtoV(res, startinx, endinx, rawData):

    temp = []  # 存储删除的峰信息
    for index in range(len(res) - 1, -1, -1):
        if res[index]["bsl"][0] < endinx and res[index]["bsl"][-1] > startinx:  # 融合峰基线与指定的位置有交集
            if len(res[index]["peak"]) != 1:
                for i in range(len(res[index]["peak"]) - 1, -1, -1):  # 先挑出谷到谷事件中的峰

                    if res[index]["peak"][i][0] >= startinx and res[index]["peak"][i][-1] <= endinx:  # 起点和终点都落在指定区间中

                        raw_list = np.asarray(rawData[res[index]["peak"][i][0]:res[index]["peak"][i][-1]])
                        bsl_list = generatebasline(rawData, res[index]["peak"][i][0], res[index]["peak"][i][-1] - 1)
                        dif = raw_list - bsl_list < 0
                        bslstart = res[index]["peak"][i][0]
                        bslend = res[index]["peak"][i][-1]
                        while dif.any():  # 原数据有部分在基线下方
                            if dif[0:len(dif) // 2].any():  # 基线的前半部分有在原数据下方的点
                                bslstart += 1
                                raw_list = np.asarray(rawData[bslstart:res[index]["peak"][i][-1]])
                                bsl_list = generatebasline(rawData, bslstart, res[index]["peak"][i][-1] - 1)
                                dif = raw_list - bsl_list < 0
                            elif dif[len(dif) // 2:-1].any():
                                bslend -= 1
                                raw_list = np.asarray(rawData[res[index]["peak"][i][0]:bslend])
                                bsl_list = generatebasline(rawData, res[index]["peak"][i][0], bslend - 1)
                                dif = raw_list - bsl_list < 0

                        if bslstart >= res[index]["peak"][i][1]:
                            bslstart = res[index]["peak"][i][0]
                        if bslend <= res[index]["peak"][i][1]:
                            bslend = res[index]["peak"][i][-1]

                        curr = {"peak": [[bslstart, res[index]["peak"][i][1], bslend]],
                                "bsl": [bslstart, bslend],
                                "bslType":"auto"}
                        temp.append(curr)  # 临时list存储被删除的峰
                        del res[index]["peak"][i]

                if len(res[index]["peak"]) > 1:  # 对一个峰簇中剩下的峰进行处理
                    curr = {"peak": [], "bsl": None,"bslType":"auto"}
                    curr["peak"].append(res[index]["peak"][-1])
                    firststart = res[index]["peak"][-1][0]  # 最后一个峰的起点
                    del res[index]["peak"][-1]
                    for i in range(len(res[index]["peak"]) - 1, -1, -1):
                        if firststart == res[index]["peak"][i][-1]:  # 判断后一个峰的起点和前一个终点是否相等
                            curr["peak"].append(res[index]["peak"][i])
                        else:
                            curr["peak"].sort()
                            # globallow = curr["peak"][0][0] + np.argmin(rawData[curr["peak"][0][0]:curr["peak"][-1][-1]])
                            raw_list = np.asarray(rawData[curr["peak"][0][0]:curr["peak"][-1][-1]])
                            bsl_list = generatebasline(rawData, curr["peak"][0][0], curr["peak"][-1][-1] - 1)
                            globallow = curr["peak"][0][0] + np.argmin(raw_list - bsl_list)
                            if rawData[globallow] < bsl_list[globallow - curr["peak"][0][0]]:
                                previouspart = {"peak": [], "bsl": None,"bslType":"auto"}
                                behindpart = {"peak": [], "bsl": None,"bslType":"auto"}
                                for j in range(len(curr["peak"])):
                                    if curr["peak"][j][1] >= curr["peak"][0][0] and curr["peak"][j][1] <= globallow:
                                        previouspart["peak"].append(curr["peak"][j])

                                    elif curr["peak"][j][1] >= globallow and curr["peak"][j][1] <= curr["peak"][-1][-1]:
                                        behindpart["peak"].append(curr["peak"][j])
                                if len(previouspart["peak"]) != 0:
                                    previouspart["peak"][-1][-1] = globallow
                                    previouspart["bsl"] = [curr["peak"][0][0], globallow]
                                    temp.append(previouspart)
                                if len(behindpart["peak"]) != 0:
                                    behindpart["peak"][0][0] = globallow
                                    behindpart["bsl"] = [globallow, curr["peak"][-1][-1]]
                                    temp.append(behindpart)

                            else:
                                curr["bsl"] = [curr["peak"][0][0], curr["peak"][-1][-1]]
                                temp.append(curr)

                            curr = {"peak": [], "bsl": None,"bslType":"auto"}
                            curr["peak"].append(res[index]["peak"][i])

                        firststart = res[index]["peak"][i][0]
                        del res[index]["peak"][i]

                    # 谷到谷事件，可能选取的一个峰簇中中间的某部分，因为采取的从后到前的遍历
                    # 方法，所以下面的代码处理的是谷到谷前面的融合峰
                    curr["peak"].sort()
                    raw_list = np.asarray(rawData[curr["peak"][0][0]:curr["peak"][-1][-1]])
                    bsl_list = generatebasline(rawData, curr["peak"][0][0], curr["peak"][-1][-1] - 1)
                    globallow = curr["peak"][0][0] + np.argmin(raw_list - bsl_list)
                    if rawData[globallow] < bsl_list[globallow - curr["peak"][0][0]]:
                        previouspart = {"peak": [], "bsl": None,"bslType":"auto"}
                        behindpart = {"peak": [], "bsl": None,"bslType":"auto"}
                        for j in range(len(curr["peak"])):
                            if curr["peak"][j][1] >= curr["peak"][0][0] and curr["peak"][j][1] <= globallow:
                                previouspart["peak"].append(curr["peak"][j])

                            elif curr["peak"][j][1] >= globallow and curr["peak"][j][1] <= curr["peak"][-1][-1]:
                                behindpart["peak"].append(curr["peak"][j])
                        if len(previouspart["peak"]) != 0:
                            previouspart["peak"][-1][-1] = globallow
                            previouspart["bsl"] = [curr["peak"][0][0], globallow]
                            temp.append(previouspart)
                        if len(behindpart["peak"]) != 0:
                            behindpart["peak"][0][0] = globallow
                            behindpart["bsl"] = [globallow, curr["peak"][-1][-1]]
                            temp.append(behindpart)
                    else:
                        curr["bsl"] = [curr["peak"][0][0], curr["peak"][-1][-1]]
                        temp.append(curr)

                if len(res[index]["peak"]) == 1:
                    curr = {"peak": [res[index]["peak"][0]],
                            "bsl": [res[index]["peak"][0][0], res[index]["peak"][0][-1]],
                            "bslType":"auto"}
                    temp.append(curr)

                del res[index]

    # 将删除的峰添加回去
    for i in range(len(temp)):
        res.append(temp[i])

    res = sorted(res, key=lambda x: x["bsl"][0])  # 按照基线的起始点排序

    return res