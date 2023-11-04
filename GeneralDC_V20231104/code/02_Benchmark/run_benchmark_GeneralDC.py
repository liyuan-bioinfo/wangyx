# 基于MouseSpleen的数据, 模拟不同的情况
# 使用command exe的形式执行不同的benchmark方法
# Input的数据为: 
#   1, 基于single-cell proteomics train好的model或者single-cell proteomics H5AD格式文件; 
#   2, 基于Bulk protemics(需要用single-cell 进行对齐)得到的H5AD格式文件. 

# Output的数据类型为: 
#   1, cell-type abundance, ID为spot, 产生 Dataset * Method个数的结果
#   2, 每个数据集 每个方法都会产生一个结果文件



import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import pandas as pd
import scanpy as sc
import numpy as np


import glob
import configparser
import warnings


import benchmark_scripts.BenchmarkMethods as DeconvolutionSpot

warnings.filterwarnings("ignore")

celltype_key = "celltype"

save_dir = "/aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/write/01_Benchmark_Deconvolution/"

Methods = ['GeneralDC']#'SpatialDC','Tangram',DestVI,Cell2location,Stereoscope

    
sc_adata_path = "/aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/data/Lung_celltype9_sig_counts.h5ad"
# spatial_adata_path = "/aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/data/Lung_pseudobulk_counts.h5ad"
Bulk_adata_path = "/aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/simulated/Lung_sim_Similar_0.h5ad"

Bulk_data_name = Bulk_adata_path.split('/')[-1].replace(".h5ad","")

# Run de-convonlution methods
for method in Methods:

    print(f"{Bulk_data_name} {method}")
    temp_output_dir = os.path.join(save_dir,method)
    if not os.path.exists(temp_output_dir):
        os.makedirs(temp_output_dir)

    temp_output_path = os.path.join(temp_output_dir, Bulk_data_name) # prefix
    test = DeconvolutionSpot.Deconvolutions(RNA_h5ad = sc_adata_path, Bulk_h5ad = Bulk_adata_path, 
                                            celltype_key = celltype_key, 
                                            output_path = temp_output_path)    
    test.Dencon(method)

