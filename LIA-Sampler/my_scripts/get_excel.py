import os
import pandas as pd

def extract_percentage(filepath):
    """从文件中读取百分数"""
    with open(filepath, 'r') as file:
        content = file.read().strip()
    return content

def scan_directory(directory):
    """扫描目录并提取所有文件的百分数"""
    data = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".out"):  # 假设文件是文本文件，你可以根据需要更改
                filepath = os.path.join(root, filename)
                percentage = extract_percentage(filepath)
                filename_without_ext = os.path.splitext(filename)[0]
                data.append([filename_without_ext, percentage])
    return data

def save_to_excel(data, output_file):
    """将数据保存到Excel文件"""
    df = pd.DataFrame(data, columns=['Filename', 'Percentage'])
    df.to_excel(output_file, index=False)

# 使用示例
directory = 'calc_res'  # 替换为你的文件夹路径
output_file = 'output.xlsx'  # 输出Excel文件名
data = scan_directory(directory)
save_to_excel(data, output_file)
