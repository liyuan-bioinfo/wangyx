#'@desc 用于Lung transcriptomics 数据集的结果分析及可视化
#'@update_history 20231104: bulid script
#'

library(dplyr)
library(ggplot2)


rm(list=ls())
setwd("/aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/visualization/Lung")
# QUANTIFIED 
# combine output as table
output_names = c("/aaa/zihanwu/yyyli2/projectx_General_Deconv/01_datasets/write/01_Benchmark_Deconvolution/")
methods = c("GeneralDC","Tangram","GroundTruth")

output_df = data.frame()

for (method in methods){
  for (output_name in output_names){ # different condition
    temp_files = list.files(paste0(output_name,method),pattern = "*.Normal_0.csv",full.names = T)
    temp_filesname = list.files(paste0(output_name,method),pattern = "*.Normal_0.csv",full.names = F)
    
    for (i in 1:length(temp_files)){
      
      temp_full_filename = temp_files[i]
      temp_filename = temp_filesname[i]
      
      temp_file_df = read.csv(temp_full_filename,header = T,row.names = 1)  
      temp_file_df$SampleID = row.names(temp_file_df)
      temp_file_df$filename = temp_filename
      temp_file_df$method = method
      row.names(temp_file_df) = NULL
      
      output_df = rbind(output_df, temp_file_df)            
    }
  }
}

output_df$filename = NULL

# write.csv(output_df,file=paste0("raw/Merged_Prediction_",Sys.Date(),".csv"))

pred_df_gather = output_df %>% tidyr::gather(key="celltype",value="pred_prec",-c("SampleID","method"))

gd_df_gather = pred_df_gather %>% dplyr::filter(method == "GroundTruth")
gd_df_gather$ID = paste0(gd_df_gather$SampleID,"_",gd_df_gather$celltype)
gd_df_gather$gd_prec = gd_df_gather$pred_prec
gd_df_gather = gd_df_gather %>% dplyr::select(ID, gd_prec)

pred_df_gather = pred_df_gather %>% dplyr::filter(method != "GroundTruth")
pred_df_gather$ID = paste0(pred_df_gather$SampleID,"_",pred_df_gather$celltype)

merged_pred_df = merge(pred_df_gather, gd_df_gather,by="ID",all.x = F,all.y = F) %>% unique()

# write.csv(merged_pred_df,file = paste0("noise/RMSE_noise_Method1_",Sys.Date(),".csv"))
# --------------------------------------------------------------------------------------------------
# RMSE
rmse = c()
for (i in 1:dim(merged_pred_df)[1]){
  x = as.numeric(merged_pred_df$pred_prec[i])
  y = as.numeric(merged_pred_df$gd_prec[i])
  rmse = c(rmse, sqrt(mean((x - y)^2)))
}

merged_pred_df$rmse = rmse

write.csv(merged_pred_df,file = paste0("raw/Pseudo_bulk_Lung_Metric_Normal_",Sys.Date(),".csv"))


# write.csv(merged_pred_df,file = paste0("raw/Application_HumanTonsil2023_RMSE_Method5_",Sys.Date(),".csv"))

p2_rmse_1=merged_pred_df %>%#%>% dplyr::filter(celltype=="B") %>%
  ggplot(aes(x=method,y=rmse,fill=method)) + 
  geom_boxplot() + 
  # geom_jitter(width = 0.5,alpha = 0.2) +
  # facet_wrap(~celltype,ncol = 3) +
  theme_bw()+
  theme(axis.text.x = element_blank(),plot.background = element_blank(),
        legend.background = element_blank(),panel.grid = element_blank(),
        text = element_text(size = 16,face = "bold"),axis.ticks.x.bottom = element_blank()
  )+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25),label=seq(0,1,0.25),expand = c(0,0.02))+
  labs(x = "", y = "RMSE")

pdf(file = paste0("Boxplot_RMSE_Pseudo_bulk_Lung_Normal_",Sys.Date(),".pdf"),width = 8,height = 8)
print(p2_rmse_1)
dev.off()


# Other Metric
### CCC, PCC, SPCC
pcc = c()
spcc = c()
ccc = c()
out_df = data.frame()

    for (c in unique(merged_pred_df$celltype)){
      for (m in unique(merged_pred_df$method)){
        temp_df = merged_pred_df %>% 
          filter(method==m) %>% filter(celltype==c) 
        x = as.numeric(temp_df$pred_prec)
        y = as.numeric(temp_df$gd_prec)
        # pcc = c(pcc, rep(cor(x,y),length(x)))
        # spcc = c(spcc, rep( cor(x,y,method = "spearman"),length(x)) )
        # ccc = c(ccc, rep(DescTools::CCC(x, y)$rho.c$est,length(x)))
        pcc = rep(cor(x,y),length(x))
        spcc = rep( cor(x,y,method = "spearman"),length(x)) 
        ccc = rep(DescTools::CCC(x, y)$rho.c$est,length(x))
        temp_df$pcc = pcc
        temp_df$spcc = spcc
        temp_df$ccc = ccc
        out_df = rbind(out_df, temp_df)
        
        
      }
    }
  
metric_df = out_df %>% dplyr::select(method, celltype,pcc,spcc,ccc) %>% unique()
write.csv(metric_df,file = paste0("Pseudo_bulk_Lung_Metric_Method5_",Sys.Date(),".csv"))
#----------------------------------------------
# PLOT, ALL
# PCC
# metric_df$noise = factor(metric_df$noise, levels= names(table(metric_df$noise)))
metric_df$celltype = factor(metric_df$celltype, levels= names(table(metric_df$celltype)))
metric_df$method = factor(metric_df$method, levels= c("GeneralDC", "Tangram"))


p2_pcc_1=metric_df %>%#%>% dplyr::filter(celltype=="B") %>%
  ggplot(aes(x=method,y=pcc,fill=method)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.5,alpha = 0.2) +
  # facet_wrap(~celltype,ncol = 3) +
  theme_bw()+
  theme(axis.text.x = element_blank(),plot.background = element_blank(),
        legend.background = element_blank(),panel.grid = element_blank(),
        text = element_text(size = 16,face = "bold"),axis.ticks.x.bottom = element_blank()
  )+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25),label=seq(0,1,0.25),expand = c(0,0.02))+
  labs(x = "", y = "PCC")

pdf(file = paste0("Boxplot_PCC_Pseudo_bulk_Lung_Normal_",Sys.Date(),".pdf"),width = 4,height = 4)
print(p2_pcc_1)
dev.off()

p2_ccc_1=metric_df %>% 
  ggplot(aes(x=method,y=ccc,fill=method)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.5,alpha = 0.2) +
  # facet_wrap(~celltype,ncol = 2) +
  theme_bw()+
  theme(axis.text.x = element_blank(),plot.background = element_blank(),
        legend.background = element_blank(),panel.grid = element_blank(),
        text = element_text(size = 16,face = "bold"),axis.ticks.x.bottom = element_blank()
  )+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25),label=seq(0,1,0.25),expand = c(0,0.02))+
  labs(x = "", y = "CCC")

pdf(file = paste0("Boxplot_CCC_Pseudo_bulk_Lung_Normal_",Sys.Date(),".pdf"),width = 4,height = 4)
print(p2_ccc_1)
dev.off()




