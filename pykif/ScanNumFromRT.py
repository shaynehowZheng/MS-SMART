import pandas as pd
import os
from tqdm import tqdm
# 引入 ctypes 以支持 GetPrecursorMassForScanNum 的参数类型 (c_double, byref)
from ctypes import c_double, byref 

# 假设 MSFileReader 已经正确导入或定义，且包含 ScanNumFromRT 和 GetPrecursorMassForScanNum 方法
# 如果您使用的是 pymsfilereader 或类似库，请确保这些方法和 ctypes 参数的使用是匹配的。
from MSFileReader import MSFileReader


# --- 配置路径 ---
INPUT_CSV_PATH = r"C:\work\msclassifier\logs\20250627BanMoNA\multi_predictions_Glycyrrhiza uralensis.csv"
RAW_FILE_BASE_DIR = r"C:\work\msclassifier\data\Glycyrrhiza uralensis"
OUTPUT_CSV_PATH = INPUT_CSV_PATH.replace('.csv', '_accurate_scannum.csv')

# --- 常量配置 ---
RT_COLUMN_PREFIX = 'Rt_' 
MZ_COLUMN_PREFIX = 'Mz_' # 新增：用于构造 Mz 列前缀
NEW_SCANNUM_COLUMN = 'Accurate_ScanNum'
PPM_TOLERANCE = 20.0 # 20 ppm 容忍度
RT_STEP = 0.001 # 每次调整的 Rt 步长 (分钟)
MAX_RT_TRIES = 100 # 最大尝试次数 (控制搜索范围)
MS_ORDER_TO_CHECK = 2 # 假设我们要查找的是 MS2 谱图的前体离子

def calculate_ppm_error(mz_theoretical, mz_observed):
    """计算两个 m/z 值之间的 ppm 误差"""
    if mz_theoretical == 0:
        return float('inf')
    return abs(mz_observed - mz_theoretical) * 1e6 / mz_theoretical

def find_scannum_by_mz_walk(raw_reader, initial_rt, theoretical_mz):
    """
    以初始 RT 为起点，微调 RT 并在 RAW 文件中搜索匹配 m/z 的扫描号。
    
    Args:
        raw_reader (MSFileReader): 已加载的 RAW 文件读取器实例。
        initial_rt (float): 初始保留时间。
        theoretical_mz (float): 表格中的前体离子 m/z。
        
    Returns:
        int: 找到的准确扫描号，如果未找到则返回 None。
    """
    best_scan_num = None
    min_ppm_error = float('inf')
    
    # 定义搜索步长和方向
    search_steps = [0] + \
                   [i * RT_STEP for i in range(1, MAX_RT_TRIES)] + \
                   [-i * RT_STEP for i in range(1, MAX_RT_TRIES)]
    search_steps = sorted(list(set(search_steps)), key=abs) # 保证先尝试接近 0 的步长
    
    for step in search_steps:
        current_rt = initial_rt + step
        
        try:
            # 1. 根据当前 RT 获取扫描号
            current_scan_num = raw_reader.ScanNumFromRT(current_rt)
            
            if current_scan_num is None or current_scan_num <= 0:
                continue

            # 2. 从 RAW 文件中获取该扫描号的前体离子 m/z
            # 假设我们只关心 MS2 扫描的前体离子 (MSOrder=2)
            observed_mz = raw_reader.GetPrecursorMassForScanNum(
                current_scan_num, MS_ORDER_TO_CHECK
            )
            
            # 3. 计算误差并确认
            ppm_error = calculate_ppm_error(theoretical_mz, observed_mz)
            
            if ppm_error <= PPM_TOLERANCE:
                # 找到准确匹配，立即返回
                return current_scan_num
            
            # 记录最小误差，以防找不到完美匹配
            if ppm_error < min_ppm_error:
                min_ppm_error = ppm_error
                best_scan_num = current_scan_num

        except Exception as e:
            # 忽略 GetPrecursorMassForScanNum 抛出的错误，可能只是当前扫描不是 MS2
            # print(f"RT {current_rt:.3f} 查找失败: {e}")
            continue

    # 如果搜索完成仍未找到满足 PPM_TOLERANCE 的扫描号，返回 None
    return None

def get_scannum_from_rt_by_raw(input_csv_path, base_raw_dir):
    """
    根据 Filename、Rt 和 Mz 进行精确的扫描号查找和验证。
    """
    print(f"--- 1. 读取数据文件: {input_csv_path}")
    try:
        df = pd.read_csv(input_csv_path)
    except FileNotFoundError:
        print(f"错误: 找不到文件 {input_csv_path}")
        return
    
    if 'Filename' not in df.columns:
        print("错误: CSV 文件中缺少 'Filename' 列。")
        return
    
    # 确保目标列存在，并初始化为 None
    filename_col_index = df.columns.get_loc('Filename')
    df.insert(filename_col_index + 1, NEW_SCANNUM_COLUMN, None)
    
    # 存储已加载的 MSFileReader 实例，避免重复加载大型 RAW 文件
    raw_reader_cache = {}
    
    print(f"--- 2. 按 Filename 分组并处理 ---")
    unique_filenames = df['Filename'].dropna().unique()

    for filename in tqdm(unique_filenames, desc="处理 RAW 文件"):
        
        # 🌟 动态构造 Rt 列名和 Mz 列名
        rt_col_name = RT_COLUMN_PREFIX + filename
        mz_col_name = MZ_COLUMN_PREFIX + filename # 新增：Mz 列名
        
        if rt_col_name not in df.columns or mz_col_name not in df.columns:
            print(f"警告: 缺少 Rt 列 ('{rt_col_name}') 或 Mz 列 ('{mz_col_name}')。跳过该文件相关的数据。")
            continue
            
        raw_file_reader = None
        file_path = None
        try:
            # 3. 构造文件路径并加载 MSFileReader
            raw_file_name = filename + '.raw'
            file_path = os.path.join(base_raw_dir, raw_file_name)
            
            if not os.path.exists(file_path):
                print(f"警告: RAW 文件未找到: {file_path}。跳过该文件。")
                continue
            
            # 从缓存中获取 MSFileReader 实例
            if file_path not in raw_reader_cache:
                raw_reader_cache[file_path] = MSFileReader(file_path)
            raw_file_reader = raw_reader_cache[file_path]
            
            # 获取当前 Filename 对应的所有行索引
            indices = df.index[df['Filename'] == filename].tolist()
            
            for index in indices:
                initial_rt = df.loc[index, rt_col_name]
                theoretical_mz = df.loc[index, mz_col_name]
                
                if pd.notna(initial_rt) and pd.notna(theoretical_mz):
                    # 4. 调用新的精确查找函数
                    found_scan_num = find_scannum_by_mz_walk(
                        raw_file_reader, float(initial_rt), float(theoretical_mz)
                    )
                    
                    # 5. 保存扫描号
                    if found_scan_num is not None:
                        df.loc[index, NEW_SCANNUM_COLUMN] = found_scan_num
                    else:
                        # 如果多次尝试后仍未找到匹配的扫描号
                        df.loc[index, NEW_SCANNUM_COLUMN] = "Not Confirmed"
                
        except Exception as e:
            print(f"处理文件 {filename} 时发生致命错误: {e}")
            # 清理缓存以防加载了损坏的文件
            if file_path in raw_reader_cache and raw_file_reader:
                # 尝试关闭文件
                try: raw_file_reader.Close() 
                except: pass
                del raw_reader_cache[file_path]
            continue
            
    # --- 6. 保存最终结果 ---
    print(f"--- 6. 保存最终结果到: {OUTPUT_CSV_PATH}")
    df.to_csv(OUTPUT_CSV_PATH, index=False, encoding='utf-8-sig')
    print("✅ 任务完成。")


if __name__ == '__main__':
    # 注意：在实际运行前，请确保您的 MSFileReader 类已正确定义，
    # 且能够正确处理 ctypes 的 c_double 和 byref 参数。
    search_steps = [0] + \
                   [i * RT_STEP for i in range(1, MAX_RT_TRIES)] + \
                   [-i * RT_STEP for i in range(1, MAX_RT_TRIES)]
    search_steps = sorted(list(set(search_steps)), key=abs) # 保证先尝试接近 0 的步长
    get_scannum_from_rt_by_raw(
        input_csv_path=INPUT_CSV_PATH,
        base_raw_dir=RAW_FILE_BASE_DIR
    )
