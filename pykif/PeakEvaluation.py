import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
import re

# ----------------------------------------------------------------------
# 1. 定义高斯函数模型
# ----------------------------------------------------------------------
def gaussian(x, A, mu, sigma):
    """
    高斯函数（正态分布）模型
    A: 峰高
    mu: 峰中心（保留时间）
    sigma: 峰宽的标准差
    """
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# ----------------------------------------------------------------------
# 2. 定义拟合和评估函数
# ----------------------------------------------------------------------
def fit_and_evaluate_peak(row):
    """
    对单行色谱峰数据进行高斯拟合，并计算R^2和残差统计。

    参数:
        row (pd.Series): 包含 'x' 和 'y' 列的单行数据。

    返回:
        pd.Series: 包含拟合参数 (A, mu, sigma), R^2, AIC, BIC 和残差统计的Series。
    """
    try:
        # 提取并清理 x (时间) 和 y (信号强度) 数据
        # 从字符串表示中提取数字列表
        x_str = row['x'].strip('[]').replace('\n', ' ')
        y_str = row['y'].strip('[]').replace('\n', ' ')
        
        # 使用正则表达式匹配浮点数（包括科学计数法）
        x_data = np.array([float(val) for val in re.findall(r'[-+]?\d*\.\d+|\d+', x_str)])
        y_data = np.array([float(val) for val in re.findall(r'[-+]?\d*\.\d+|\d+', y_str)])

        if len(x_data) < 3 or len(x_data) != len(y_data):
            return pd.Series({
                'fit_A': np.nan, 'fit_mu': np.nan, 'fit_sigma': np.nan,
                'R_squared': np.nan, 'residual_mean': np.nan, 
                'residual_std': np.nan, 'AIC': np.nan, 'BIC': np.nan,
                'fit_error': 'Data_Mismatch_or_Too_Short'
            })
            
        # ------------------- 确定初始参数 (p0) -------------------
        max_y_index = np.argmax(y_data)
        A_p0 = y_data[max_y_index]  # 初始峰高 = 实际最大高度
        mu_p0 = x_data[max_y_index] # 初始峰中心 = 最大高度对应的x值
        
        # 估算 sigma: 使用峰高一半处的峰宽 (FWHM)
        # FWHM = 2 * sqrt(2 * ln(2)) * sigma ≈ 2.355 * sigma
        # sigma ≈ FWHM / 2.355
        half_max = A_p0 / 2
        # 找到峰高一半处的左右边界索引
        left_idx = np.where(y_data[:max_y_index] < half_max)[0][-1] if np.any(y_data[:max_y_index] < half_max) else 0
        right_idx = np.where(y_data[max_y_index:] < half_max)[0][0] + max_y_index if np.any(y_data[max_y_index:] < half_max) else len(x_data) - 1
        
        # 如果左右边界有效，估算FWHM和sigma
        if right_idx > left_idx:
            FWHM = x_data[right_idx] - x_data[left_idx]
            sigma_p0 = FWHM / 2.355
        else:
            # 简单估算，例如使用峰宽 (width) 列的值作为基础
            if 'width' in row and not np.isnan(row['width']):
                 sigma_p0 = row['width'] / 4 # 粗略估计
            else:
                 # 最后的默认粗略估计
                 sigma_p0 = (x_data[-1] - x_data[0]) / 10 
                 
        p0 = [A_p0, mu_p0, sigma_p0]
        
        # ------------------- 执行拟合 -------------------
        # bounds 限制参数在合理的物理范围内，防止拟合发散
        popt, pcov = curve_fit(gaussian, x_data, y_data, p0=p0, 
                               bounds=([0, x_data.min(), 0], [np.inf, x_data.max(), np.inf]))
        
        A_fit, mu_fit, sigma_fit = popt
        
        # ------------------- 拟合评估 -------------------
        # 1. R-squared (决定系数)
        y_fit = gaussian(x_data, *popt)
        ss_total = np.sum((y_data - np.mean(y_data))**2)
        ss_residual = np.sum((y_data - y_fit)**2)
        R_squared = 1 - (ss_residual / ss_total) if ss_total != 0 else np.nan
        
        # 2. 残差统计
        residuals = y_data - y_fit
        residual_mean = np.mean(residuals)
        residual_std = np.std(residuals)

        # 3. 信息准则 (AIC, BIC) - 评估模型复杂度与拟合优度
        n = len(y_data) # 数据点数量
        k = len(popt)  # 参数数量 (3: A, mu, sigma)
        
        # 计算 log-likelihood (假设残差服从正态分布)
        # s_sq = ss_residual / (n - k) # 样本方差
        # log_likelihood = -n/2 * np.log(2*np.pi) - n/2 * np.log(s_sq) - 1/(2*s_sq) * ss_residual
        # 简化版 AIC/BIC 计算 (基于最小二乘残差平方和)
        AIC = n * np.log(ss_residual/n) + 2 * k
        BIC = n * np.log(ss_residual/n) + k * np.log(n)
        
        return pd.Series({
            'fit_A': A_fit,
            'fit_mu': mu_fit,
            'fit_sigma': sigma_fit,
            'R_squared': R_squared,
            'residual_mean': residual_mean,
            'residual_std': residual_std,
            'AIC': AIC,
            'BIC': BIC,
            'fit_error': 'None'
        })
    
    except RuntimeError:
        return pd.Series({
            'fit_A': np.nan, 'fit_mu': np.nan, 'fit_sigma': np.nan,
            'R_squared': np.nan, 'residual_mean': np.nan, 
            'residual_std': np.nan, 'AIC': np.nan, 'BIC': np.nan,
            'fit_error': 'RuntimeError_No_Fit'
        })
    except Exception as e:
        return pd.Series({
            'fit_A': np.nan, 'fit_mu': np.nan, 'fit_sigma': np.nan,
            'R_squared': np.nan, 'residual_mean': np.nan, 
            'residual_std': np.nan, 'AIC': np.nan, 'BIC': np.nan,
            'fit_error': f'General_Error: {e}'
        })

# ----------------------------------------------------------------------
# 3. 主程序
# ----------------------------------------------------------------------
def process_chromatography_data(file_path):
    """
    读取色谱数据文件，进行高斯拟合和评估，并将结果保存回原路径。

    参数:
        file_path (str): CSV 文件的完整路径。
    """
    try:
        # 检查文件是否存在
        if not os.path.exists(file_path):
            print(f"错误: 文件未找到: {file_path}")
            return

        # 1. 读取数据 (跳过第一行，因为它是标题，但不是CSV标准格式)
        # 注意: 您的数据格式不太规范，需要特殊处理
        # 假设第一行是实际的列名，第二行是数据，但数据中的 x/y 是列表字符串
        df = pd.read_csv(file_path, skipinitialspace=True, engine='python', on_bad_lines='skip')
        

        # 清理列名中的序号，并确保使用正确的列名
        df.columns = [col.split('.')[-1] if '.' in col else col for col in df.columns]

        # 2. 对每一行应用拟合和评估函数
        print(f"开始处理 {len(df)} 个色谱峰...")
        
        # 应用 fit_and_evaluate_peak 函数
        fit_results = df.apply(fit_and_evaluate_peak, axis=1)
        
        # 3. 将拟合结果合并回原始 DataFrame
        df_final = pd.concat([df, fit_results], axis=1)
        
        # 4. 准备保存路径
        directory = os.path.dirname(file_path)
        base_name = os.path.basename(file_path).replace('.csv', '')
        output_file = os.path.join(directory, f'{base_name}_GaussianFitResults.csv')
        
        # 5. 保存结果
        df_final.to_csv(output_file, index=False)
        
        print(f"\n处理完成！结果已保存到：\n{output_file}")
        print("-" * 50)
        print("主要结果摘要:")
        print(f"平均 R²: {df_final['R_squared'].mean():.4f}")
        print(f"拟合失败/错误峰数: {df_final['fit_error'].astype(str).str.count('Error|No_Fit|Mismatch').sum()}")
        print("-" * 50)

    except Exception as e:
        print(f"主程序发生错误: {e}")

# ----------------------------------------------------------------------
# 4. 运行示例 (使用您提供的虚拟路径)
# ----------------------------------------------------------------------

# 虚拟文件路径（请根据您的实际文件位置修改！）
# 警告：此路径在您的系统中可能不存在，请创建它或将其更改为实际文件路径。
data_directory = r"results\20250702-YY-G2-step2-N-1\peaks_energy=unknown_ppm4=10_ppm5=10_ppm6=10_top=5.csv"
file_name = "your_data_file.csv" # 假设您已将数据保存为这个名字


# 执行主函数
process_chromatography_data(data_directory)