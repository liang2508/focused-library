# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:51:45 2023

@author: HY
"""

import pandas as pd
import matplotlib.pyplot as plt

# 读取CSV文件
data = pd.read_csv('./NIK/new_molecules_file_allmethod_nik_SA_QED.csv')
#file='./docking/docking_score_asinex40000.csv'
#data = pd.read_csv(file)
# 提取第二列和第三列数据
category_column = data.iloc[:, 8]
column_data = data.iloc[:, 7]

# 定义区间范围和步长
bins = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#bins = [-12.0,-11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
# 统计各区间的数据数量，并按照类别进行分组
hist_data, bin_edges = pd.cut(column_data, bins=bins, include_lowest=True, right=False, retbins=True)
counts = hist_data.value_counts().sort_index()

# 按照类别进行分组，统计各类别各区间的数据数量
grouped_data = data.groupby(['method', pd.cut(data.iloc[:, 7], bins=bins, include_lowest=True, right=False)]).size().unstack()

# 打印结果
print(counts)
print("\nGrouped Data:")
print(grouped_data)
