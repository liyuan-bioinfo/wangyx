
# desc: Benchmark methods
# benchmark methods: GeneralDC, Tangram, .... , et al. 
# history: 
#   Version_20231104 [Li-Yuan]: create this script for benchmark different De-Convolution methods.




import sys
import os

class Deconvolutions:
    """

        # usage
        >>> import Benchmarking.Deconvolution as Deconvolution
        >>> test = Deconvolution.Deconvolutions(RNA_h5ad, Bulk_file, Bulk_h5ad, celltype_key, python_path, output_path)
        >>> Methods = ['Tangram', 'GeneralDC']
        >>> Result = test.Dencon(Methods)
        
        Parameters
        -------
        
        RNA_h5ad : str
        scRNA-seq data file with h5ad format.
            
        Bulk_file : str
        Bulk data count file.
        
        Bulk_h5ad : str
        Bulk data file with h5ad format.
            
        celltype_key : str
        celltype annotataion title in scRNA-seq data h5ad file
        
        celltype_file : str
        celltype annotataion file
        
        python_path : str
        which python path used for Benchmark methods
        
        output_path : str
        Outfile path
        
    """



    def __init__(self, RNA_h5ad = None, Bulk_h5ad = None, celltype_key = None, python_path = None, output_path = None):

        
        self.RNA_h5ad = RNA_h5ad
        self.Bulk_h5ad = Bulk_h5ad
        self.celltype_key = celltype_key
        self.python_path = python_path
        self.output_path = output_path
    
    def Dencon(self, need_tools):
        print(need_tools)
 
        if "Stereoscope" in need_tools:
            RNA_h5ad = self.RNA_h5ad
            Bulk_h5ad = self.Bulk_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python /aaa/zihanwu/yyyli2/project1_Bulk_deconv/02_benchmark/benchmark_scripts/scripts/Stereoscope_pipeline.py ' + RNA_h5ad + ' ' + Bulk_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "Tangram" in need_tools:
            RNA_h5ad = self.RNA_h5ad
            Bulk_h5ad = self.Bulk_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python /aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/write/01_Benchmark_Deconvolution/benchmark_scripts/scripts/Tangram_pipeline.py ' + RNA_h5ad + ' ' + Bulk_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "GeneralDC" == need_tools:
            RNA_h5ad = self.RNA_h5ad
            Bulk_h5ad = self.Bulk_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python /aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/write/01_Benchmark_Deconvolution/benchmark_scripts/scripts/GeneralDC_pipeline.py ' + RNA_h5ad + ' ' + Bulk_h5ad + ' ' + celltype_key + ' ' + output_path)


        if "DestVI" in need_tools:
            RNA_h5ad = self.RNA_h5ad
            Bulk_h5ad = self.Bulk_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python /aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/write/01_Benchmark_Deconvolution/benchmark_scripts/scripts/DestVI_pipeline.py ' + RNA_h5ad + ' ' + Bulk_h5ad + ' ' + celltype_key + ' ' + output_path)




