# GeneralDC_V20231104 文件说明
## @author Li-Yuan
## @time 20231104

-code: 展示我们当前的方法的源代码以及对比方法
	-01_Simulation：
		"Simulation.ipynb" 用于模拟数据（可1，参考如何基于single-cell 的数据，进行下采样生成显著性矩阵）；2，基于single-cell的数据，构建pseduo-bulk
		"visualization_sig_matrix.R 用于将下采样后基于CIBERSROTx在线工具生成的显著性矩阵进行可视化。
		
	-02_Benchmark：
		"run_benchmark_GeneralDC.py" 通过该方法直接调用我们的模型
		"run_benchmark_methods.py" 通过该方法调用其他对比方法
		"benchmark_scripts" 包含benchmark方法的python源码，也包括直接调用的Tangram以及我们的General-DC
-data
	"Lung_celltype9_sig_matrix.txt" 基于single-cell 进行下采用+CIBERSORTx的应用后，构建的显著性矩阵。
	"Lung_sim_Normal_0.h5ad" 基于single-cell 模拟的Normal状态下的bulk transcriptome data, 即细胞比例随机分布
	"Lung_sim_Rare_0.h5ad" 基于single-cell 模拟的Rare状态下的bulk transcriptome data，即细胞比例分布不均匀
	"Lung_sim_Similar_0.h5ad" 基于single-cell 模拟的Similar状态下的bulk transcriptome data，即细胞比例分布相似

-output
	"GroundTruth" 模拟bulk data时，产生的细胞比例，用于比较不同方法的预测结果。
	"Tangram" 对比方法预测的结果
	"GerneralDc" 我们方法预测的结果
	
-visualization
	xxx.R 分别用于读取不同状态下的结果文件，进行方法的性能评估以及可视化。主要是计算误差以及相似性。
	
- other
	环境的配置：推荐用conda管理环境，配置PyTorch，Scanpy, Numpy, Pandas等，根据官网教程即可。其他主要是对比方法依赖的包。R语言安装基本的ggplot2, dplyr即可。
	模型的调用：
		以我们的方法为例：python run_benchmark_GeneralDC.py # 该方法通过1，加载Sig. Matrix, Bulk data的位置；2，命令行调用GerneralDC脚本的形式运行。
		对比方法是调用sc-RNA data的，这个文件比较大，需要另外下载。