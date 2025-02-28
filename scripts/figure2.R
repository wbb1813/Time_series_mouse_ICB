library(slingshot)
library(ggplot2)
library(Seurat)
library(reshape)
library(ggpubr)
library('cutpointr')
#library(immunarch)
## ------- Inputs and parameters -------
cd8_t_obj=readRDS('/Volumes/Binbin_Kun/binbin/silvio/results/seurat/tcell/cd8/res0.2/annotated_seurat_res0.2.rds')

tcr_meta=read.delim('../data/sc_tcr_meta.txt')
teff_subtype_average_clone_size=read.delim('../data/teff_subtype_average_clone_size.txt')
teff_subtype_extend_clone_fre=read.delim('../data/teff_subtype_extend_clone_fre.txt')

sample_average_clone_size=read.delim('../data/sample_average_clone_size.txt')

bulk_tcr_div_inv_simp=read.delim('../data/bulk_tcr_div_inv_simp.txt')
bulk_bcr_div_inv_simp=read.delim('../data/bulk_bcr_div_inv_simp.txt')

cd4_tcr_meta=read.delim('../data/sc_cd4_tcr_meta.txt')

patient_meta=read.delim('../data/LiBIO_scRNA_sample_Summary_metrics.txt')

## Outputs 
outdir='../results/figure2'
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ------- F2A: UMAP of CD8 T cells -------
## Functions 
tcr_umap=function(df,prefix){
  tcr_meta_clone1=df[which(df$frequency==1),]
  tcr_meta_clone_g1=df[which(df$frequency>1),]
  sp2<-ggplot() + 
    geom_point(data = tcr_meta_clone1, aes(x=umap_1, y=umap_2, color=frequency),size=0.2)+
    geom_point(data = tcr_meta_clone_g1, aes(x=umap_1, y=umap_2, color=frequency),size=0.2)+
    #scale_color_gradient(low="#AED6F1", high="#FC4E07")+theme_classic()
    scale_color_gradient(low="#66c2a5", high="#FC4E07")+theme_classic()
  sp2
  
  ggsave(file.path(outdir,paste0(prefix,'_clonetype.pdf')),sp2,width = 4.5,height = 3.5)
  return(sp2)
}

## UMAP of CD8 T cells of all samples 

cd8_all_tcr_umap=tcr_umap(df=tcr_meta,prefix='cd8t_all_samples')

# UMAP of CD8 T cells of responders, non-responder
responder_samples=patient_meta$Sample.Name[which(patient_meta$group=='responder')]
nonresponder_samples=patient_meta$Sample.Name[which(patient_meta$group%in%c('non-responder','stable/escape'))]

responder_tcr_meta=tcr_meta[which(tcr_meta$sample_id2%in%responder_samples),]
nonresponder_tcr_meta=tcr_meta[which(tcr_meta$sample_id2%in%nonresponder_samples),]

cd8_res_tcr_umap=tcr_umap(df=responder_tcr_meta,prefix='cd8t_responders')
cd8_non_tcr_umap=tcr_umap(df=nonresponder_tcr_meta,prefix='cd8t_nonresponder')

## ------- Distribution of expanded clones ------- 
tcr_meta_efftcell=tcr_meta[which(tcr_meta$Manually_curation=='Effector memory CD8'),]
tcr_meta_othertcell=tcr_meta[which(tcr_meta$cell_id%in%setdiff(tcr_meta$cell_id,tcr_meta_efftcell$cell_id)),]

tcr_meta_efftcell_te=tcr_meta_efftcell[which(tcr_meta_efftcell$frequency>=2),]
tcr_meta_efftcell_ne=tcr_meta_efftcell[which(tcr_meta_efftcell$frequency==1),]

tcr_meta_othertcell_te=tcr_meta_othertcell[which(tcr_meta_othertcell$frequency>=2),]
tcr_meta_othertcell_ne=tcr_meta_othertcell[which(tcr_meta_othertcell$frequency==1),]

eff_other_2b2=data.frame(te=c(nrow(tcr_meta_efftcell_te),nrow(tcr_meta_othertcell_te)),ne=c(nrow(tcr_meta_efftcell_ne),nrow(tcr_meta_othertcell_ne)),row.names = c('effect','other'))
test_res=fisher.test(eff_other_2b2,alternative = 'greater')
format.pval(test_res$p.value,digits = 20)
# barplot 
tmp_df=as.data.frame(t(eff_other_2b2))
tmp_df$group=rownames(tmp_df)

df_eff=tmp_df[,c('effect','group')]
colnames(df_eff)[1]='cell_number'
df_eff$Freq=df_eff$cell_number/sum(df_eff$cell_number)
df_eff$Group2='Tem'

df_other=tmp_df[,c('other','group')]
colnames(df_other)[1]='cell_number'
df_other$Freq=df_other$cell_number/sum(df_other$cell_number)
df_other$Group2='Other T'

df=rbind(df_eff,df_other)
df$group[which(df$group=='te')]='Expanded'
df$group[which(df$group=='ne')]='Non-expanded'
df$group=factor(df$group,c('Non-expanded','Expanded'))

p=ggplot(data=df, aes(x=Group2, y=Freq, fill=group)) +
  geom_bar(stat="identity",width = 0.7)+
  theme_classic()+xlab('')+
  geom_text(x=1, y=1.05, label=paste0('P: ',test_res$p.value))+
  geom_text(x=2, y=1.05, label=paste0('Odds: ',signif(test_res$estimate,4)))+
  #scale_fill_brewer(palette="Dark2")+"
  scale_fill_manual(values=c("#66c2a5","#fc8d62"))
  #expand_limits(y = c(0, 1.05))+theme(legend.title = element_blank() )
p
ggsave(file.path(outdir,paste0('cd8t','_clonetype_distribution.pdf')),p,width = 4,height = 4)

## ======= Clonal expansion index (average clonotype size) across time points  =======
## Functions
# boxplot compare responders with non-responders 
twogroup_boxplot=function(df,prefix, y_value,alternative='less'){
  df$timepoint=factor(df$timepoint,c('1','2','3','4'))
  df$Response=factor(df$Response,levels = c('Responders','Non-responders'))
  
  p <- ggboxplot(df, x = "timepoint", y = y_value,
                 color = "Response", palette = "jco",
                 add = "jitter")
  p = p + stat_compare_means(aes(group = Response),label = "p.signif",method = 'wilcox.test',method.args = list(alternative = alternative))+
    xlab('Timepoints')
  p
  ggsave(file.path(outdir,paste0(prefix,'_tcr.pdf')),p,width = 6,height = 4)
  return(p)
}

## CD8 effector memory T cell 
teff_clone_size_box=twogroup_boxplot(df=teff_subtype_average_clone_size,y_value = 'mean_clone_size',prefix='sc_teff_tcr_mean_clone_size')

teff_clone_fre_box=twogroup_boxplot(df=teff_subtype_extend_clone_fre,y_value = 'mean_extend_fre',prefix='sc_teff_tcr_mean_extend_fre')

## ------- Test if expanded cells enriched in responders, frequency of expanded clonotypes -------
## merge all time points, Fisher test
all_te=tcr_meta[which(tcr_meta$frequency>=2),]
all_ne=tcr_meta[which(tcr_meta$frequency==1),]
responder_te=responder_tcr_meta[which(responder_tcr_meta$frequency>=2),]
responder_ne=responder_tcr_meta[which(responder_tcr_meta$frequency==1),]
nonresponder_te=nonresponder_tcr_meta[which(nonresponder_tcr_meta$frequency>=2),]
nonresponder_ne=nonresponder_tcr_meta[which(nonresponder_tcr_meta$frequency==1),]

responder_2b2=data.frame(te=c(nrow(responder_te),nrow(nonresponder_te)),ne=c(nrow(responder_ne),nrow(nonresponder_ne)),row.names = c('responders','non-responders'))
clono_size_test=fisher.test(responder_2b2,alternative = 'greater')

# barplot 
tmp_df=as.data.frame(t(responder_2b2))
tmp_df$expanded=rownames(tmp_df)

df_responder=tmp_df[,c('responders','expanded')]
df_nonresponder=tmp_df[,c('non-responders','expanded')]

colnames(df_responder)[1]=c('cell_number')
colnames(df_nonresponder)[1]=c('cell_number')

df_responder$Freq=df_responder$cell_number/sum(df_responder$cell_number)
df_nonresponder$Freq=df_nonresponder$cell_number/sum(df_nonresponder$cell_number)

df_responder$Group='Responder'
df_nonresponder$Group='Non-Responder'

df=rbind(df_responder,df_nonresponder)

tmp_df=df[which(df$expanded=='te'),]
p=ggplot(tmp_df, aes(x=Group, y=Freq, fill=Group))+
  geom_bar(stat="identity", color="black",width = 0.7)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_classic()+
  #geom_text(x=1.5, y=0.055, label="Pvalue: 5.43e-07")+
  geom_text(x=1, y=0.055, label=paste0('P: ',signif(clono_size_test$p.value,3)))+
  geom_text(x=2, y=0.055, label=paste0('Odds: ',signif(clono_size_test$estimate,3)))+
  xlab('')+ylab('Frequncy of expanded clonetype')+theme(legend.position = 'none')
p
ggsave(file.path(outdir,paste0('cd8t','_clonetype_expantion.pdf')),p,width = 3,height = 4)

## ------- Bulk TCR -------
bulk_tcr_simp=twogroup_boxplot(df=bulk_tcr_div_inv_simp,y_value = 'Value',prefix='bulk_tcr_simp')

## ------- Bulk BCR -------
bulk_bcr_simp=twogroup_boxplot(df=bulk_bcr_div_inv_simp,y_value = 'Value',prefix='bulk_bcr_simp')

## ------- UMAP of CD4 T cells -------
cd4_all_tcr_umap=tcr_umap(df=cd4_tcr_meta,prefix='cd4t_all_samples')

# UMAP of CD4 T cells of responders, non-responder
responder_samples=patient_meta$Sample.Name[which(patient_meta$group=='responder')]
nonresponder_samples=patient_meta$Sample.Name[which(patient_meta$group%in%c('non-responder','stable/escape'))]

responder_tcr_meta=cd4_tcr_meta[which(cd4_tcr_meta$sample_id2%in%responder_samples),]
nonresponder_tcr_meta=cd4_tcr_meta[which(cd4_tcr_meta$sample_id2%in%nonresponder_samples),]

cd4_res_tcr_umap=tcr_umap(df=responder_tcr_meta,prefix='cd4t_responders')
cd4_non_tcr_umap=tcr_umap(df=nonresponder_tcr_meta,prefix='cd4t_nonresponder')

## ------- Distribution of expanded CD4 clones ------- 
th1_tcr_meta=cd4_tcr_meta[which(cd4_tcr_meta$Manually_curation=='Th1'),]
othercd4_tcr_meta=cd4_tcr_meta[which(cd4_tcr_meta$cell_id%in%setdiff(cd4_tcr_meta$cell_id,th1_tcr_meta$cell_id)),]

th1_tcr_te=th1_tcr_meta[which(th1_tcr_meta$frequency>=2),]
th1_tcr_ne=th1_tcr_meta[which(th1_tcr_meta$frequency==1),]

othercd4_tcr_te=othercd4_tcr_meta[which(othercd4_tcr_meta$frequency>=2),]
othercd4_tcr_ne=othercd4_tcr_meta[which(othercd4_tcr_meta$frequency==1),]

th1_other_2b2=data.frame(te=c(nrow(th1_tcr_te),nrow(othercd4_tcr_te)),ne=c(nrow(th1_tcr_ne),nrow(othercd4_tcr_ne)),row.names = c('th1','other'))
test_res=fisher.test(th1_other_2b2,alternative = 'greater')

# barplot 
tmp_df=as.data.frame(t(th1_other_2b2))
tmp_df$group=rownames(tmp_df)

df_th1=tmp_df[,c('th1','group')]
colnames(df_th1)[1]='cell_number'
df_th1$Freq=df_th1$cell_number/sum(df_th1$cell_number)
df_th1$Group2='Th1'

df_other=tmp_df[,c('other','group')]
colnames(df_other)[1]='cell_number'
df_other$Freq=df_other$cell_number/sum(df_other$cell_number)
df_other$Group2='Other T'

df=rbind(df_th1,df_other)
df$group[which(df$group=='te')]='Expanded'
df$group[which(df$group=='ne')]='Non-expanded'
df$group=factor(df$group,c('Non-expanded','Expanded'))

cd4_tcr_dis=ggplot(data=df, aes(x=Group2, y=Freq, fill=group)) +
  geom_bar(stat="identity",width = 0.7)+
  theme_classic()+xlab('')+
  geom_text(x=1, y=1.05, label=paste0('P: ',test_res$p.value))+
  geom_text(x=2, y=1.05, label=paste0('Odds: ',signif(test_res$estimate,4)))+
  #scale_fill_brewer(palette="Dark2")+
  #scale_fill_manual(values=c("#1B9E77","#E7298A"))+
  scale_fill_manual(values=c("#66c2a5","#FC4E07"))+
  expand_limits(y = c(0, 1.05))+theme(legend.title = element_blank() )
cd4_tcr_dis
ggsave(file.path(outdir,paste0('cd4t','_clonetype_distribution.pdf')),cd4_tcr_dis,width = 4,height = 4)

## ------- Calculate AUC with T and B cell expansion index -------
## Functions
cal_auc_odd=function(df,score_col,response_col,group_col,timepoint_col){
  df=df[,c(score_col,response_col,group_col,timepoint_col)]
  colnames(df)=c('score','Response','group','timepoint')
  
  pred_res=data.frame()
  for (i in unique(df$group)){
    for (j in unique(df$timepoint)){
      tmp_df=df[which(df$group==i&df$timepoint==j),]
      if (length(unique(tmp_df$Response))==2){
        opt_cut <- cutpointr(tmp_df, score, Response, direction = ">=", pos_class = "1",
                             neg_class = "0", method = maximize_metric, metric = sum_sens_spec)
        opt_cut$group=i
        opt_cut$timepoint=j
        pred_res=rbind(pred_res,opt_cut)
      }
    }
  }
  return(pred_res)
}

## Calculate AUC with single cell Tem TCR data 
df=teff_subtype_average_clone_size
df$treatment='aPD1'
df$Response=ifelse(df$Response=='Responders','1',"0")

sc_auc=cal_auc_odd(df,score_col='mean_clone_size',response_col='Response',group_col='treatment',timepoint_col='timepoint')

## Calculate AUC with single cell Tem TCR data 
df=sample_average_clone_size
df$treatment='aPD1'
df$Response=ifelse(df$group=='responder','1',"0")

sc_allt_auc=cal_auc_odd(df,score_col='mean_clone_size',response_col='Response',group_col='treatment',timepoint_col='timepoint')

## Calculate AUC with bulk TCR and BCR data 
# TCR
df=bulk_tcr_div_inv_simp
df$treatment='aPD1'
df$Response=ifelse(df$Response=='Responders','1',"0")

bulk_tcr_auc=cal_auc_odd(df,score_col='Value',response_col='Response',group_col='treatment',timepoint_col='timepoint')

# BCR
df=bulk_bcr_div_inv_simp
df$treatment='aPD1'
df$Response=ifelse(df$Response=='Responders','1',"0")

bulk_bcr_auc=cal_auc_odd(df,score_col='Value',response_col='Response',group_col='treatment',timepoint_col='timepoint')

## plot AUC
sc_auc$Group='sc_Tem_TCR'
sc_allt_auc$Group='sc_all_Tcell_TCR'
bulk_tcr_auc$Group='bulk_TCR'
bulk_bcr_auc$Group='bulk_BCR'

df_plot=rbind(sc_auc,bulk_tcr_auc,bulk_bcr_auc)


p = ggplot(data=df_plot, aes(x=timepoint, y=AUC, color=Group,fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.9),width = 0.8,color="black")+
  geom_text(aes(label=signif(AUC,2)), vjust=-1, color="black",
            position = position_dodge(0.9), size=2.5)+
  scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette="Dark2")+
  theme_classic()
p
ggsave(file.path(outdir,paste0('tcr','_auc_mouse.pdf')),p,width = 6,height = 3.5)

## Add sc All T cells 
df_plot=rbind(sc_auc,sc_allt_auc,bulk_tcr_auc,bulk_bcr_auc)

p = ggplot(data=df_plot, aes(x=timepoint, y=AUC, color=Group,fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.9),width = 0.8,color="black")+
  geom_text(aes(label=signif(AUC,2)), vjust=-1, color="black",
            position = position_dodge(0.9), size=2.5)+
  scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette="Dark2")+
  theme_classic()
p
ggsave(file.path(outdir,paste0('tcr','all_Tcell_auc_mouse.pdf')),p,width = 7,height = 3.5)

## ------- ROC curve for each timepoint -------
# Function
library(pROC)
library(RColorBrewer)
roc_plot=function(df_score){
  
  df_plot=df_score[,c('filtered_t2_B cell','filtered_t2_Effector memory CD8','mean_score','Response')]
  colnames(df_plot)=c('Sig_B','Sig_Teff','Comb_sig','Response')
  
  # Create a list to store ROC curves
  roc_list <- list(
    Sig_B = pROC::roc(df_plot$Response, df_plot$Sig_B),
    Sig_Teff = pROC::roc(df_plot$Response, df_plot$Sig_Teff),
    Comb_sig = pROC::roc(df_plot$Response, df_plot$Comb_sig)
  )
  
  # Calculate AUC and update names
  auc_values <- sapply(roc_list, auc)
  names(roc_list) <- paste0(names(roc_list), " (AUC: ", round(auc_values, 2), ")")
  
  # Define custom colors using a color palette
  color_palette <- brewer.pal(n = length(roc_list), name = "Set1")
  custom_colors <- setNames(color_palette, names(roc_list))
  
  # Plot with ggroc and custom colors
  p = ggroc(roc_list) +
    ggtitle("ROC Curves with AUC Values") +
    theme_classic2() +
    labs(color = "Model") +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  return(p)
}

# 
comb_t="filtered_t2_Effector memory CD8"
comb_b="filtered_t2_B cell"

# Luoma dataset
sc_patient_icb=readRDS('../results/ICB_ML/datasets/patient_sc_dat.rds')
sc_luoma_pbmc_score=sig_score_func(expr_df=as.matrix(sc_patient_icb$Luoma_pbmc$Pseudobulk),meta=sc_patient_icb$Luoma_pbmc$meta,sig_ls=clean_sig_ls$human_sig_ls,comb_score = T,comb_t = comb_t,comb_b = comb_b)

df_score=teff_subtype_average_clone_size
bulk_tcr_div_inv_simp
bulk_bcr_div_inv_simp

## bulk TCR and BCR
tmp_tcr=bulk_tcr_div_inv_simp[,c('Sample','Value','timepoint','Response')]
tmp_bcr=bulk_bcr_div_inv_simp[,c('Sample','Value')]
colnames(tmp_tcr)[2]='TCR_clonality'
colnames(tmp_bcr)[2]='BCR_clonality'

df_score=merge(tmp_tcr,tmp_bcr,by='Sample')

for (i in unique(df_score$timepoint)){
  tmp_df_score=df_score[which(df_score$timepoint==i),]
  # Create a list to store ROC curves
  roc_list <- list(
    TCR_clonality = pROC::roc(tmp_df_score$Response, tmp_df_score$TCR_clonality),
    BCR_clonality = pROC::roc(tmp_df_score$Response, tmp_df_score$BCR_clonality)
  )
  # Calculate AUC and update names
  auc_values <- sapply(roc_list, auc)
  names(roc_list) <- paste0(names(roc_list), " (AUC: ", round(auc_values, 2), ")")
  
  # Define custom colors using a color palette
  # color_palette <- brewer.pal(n = length(roc_list), name = "Set2")
  #color_palette = c("#66C2A5","#FC8D62")
  color_palette = c("#E41A1C","#377EB8") 
  custom_colors <- setNames(color_palette, names(roc_list))
  # Plot with ggroc and custom colors
  p = ggroc(roc_list) +
    ggtitle("Bulk RNAseq TCR and BCR") +
    theme_classic() +
    labs(color = "Model") +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  p
  ggsave(file.path(outdir,paste0('ROC_AUC_T',i,'.pdf')),width = 5,height = 3)
}

## Single cell TCR of Tem
df_score=teff_subtype_average_clone_size
for (i in unique(df_score$timepoint)){
  tmp_df_score=df_score[which(df_score$timepoint==i),]
  # Create a list to store ROC curves
  roc_list <- list(
    TCR_clonality = pROC::roc(tmp_df_score$Response, tmp_df_score$mean_clone_size)
  )
  # Calculate AUC and update names
  auc_values <- sapply(roc_list, auc)
  names(roc_list) <- paste0(names(roc_list), " (AUC: ", round(auc_values, 2), ")")
  
  # Define custom colors using a color palette
  #color_palette <- brewer.pal(n = length(roc_list), name = "Set1")
  color_palette = "#4DAF4A"
  custom_colors <- setNames(color_palette, names(roc_list))
  
  # Plot with ggroc and custom colors
  p = ggroc(roc_list) +
    ggtitle("Single cell TCR") +
    theme_classic() +
    labs(color = "Model") +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  p
  ggsave(file.path(outdir,paste0('sc_ROC_AUC_T',i,'.pdf')),width = 5,height = 3)
}

## Single cell TCR of Tem and All T cells
tmp_df=teff_subtype_average_clone_size[,c('sample_id','mean_clone_size')]
colnames(tmp_df)[2]='Tem_clone_size'
colnames(sample_average_clone_size)[1]='Tcell_clone_size'
df_score=merge(tmp_df,sample_average_clone_size,by='sample_id')
df_score$Response=ifelse(df_score$group=='responder','1','0')

for (i in unique(df_score$timepoint)){
  tmp_df_score=df_score[which(df_score$timepoint==i),]
  # Create a list to store ROC curves
  roc_list <- list(
    Tem_clonality = pROC::roc(tmp_df_score$Response, tmp_df_score$Tem_clone_size),
    Tcell_clonality = pROC::roc(tmp_df_score$Response, tmp_df_score$Tcell_clone_size)
  )
  # Calculate AUC and update names
  auc_values <- sapply(roc_list, auc)
  names(roc_list) <- paste0(names(roc_list), " (AUC: ", round(auc_values, 2), ")")
  
  # Define custom colors using a color palette
  # color_palette <- brewer.pal(n = length(roc_list), name = "Set2")
  #color_palette = c("#66C2A5","#FC8D62")
  color_palette = c("#4DAF4A","#E78AC3") 
  custom_colors <- setNames(color_palette, names(roc_list))
  # Plot with ggroc and custom colors
  p = ggroc(roc_list) +
    ggtitle("Bulk RNAseq TCR and BCR") +
    theme_classic() +
    labs(color = "Model") +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  p
  ggsave(file.path(outdir,paste0('ROC_AUC_Tem_Tcell',i,'.pdf')),width = 5,height = 3)
}

## ------- Expression markers of transit T cell markers -------
features=c('Cd101','Havcr2','Tcf7','Cx3cr1','Tbx21','Gzmb','Mki67')
p=RidgePlot(cd8_t_obj$RNA, features = features, ncol = 2,group.by = 'Manually_curation',log = T)
ggsave(file.path(outdir,paste0("T_transit_makers_ridge.pdf")), p,  width=16, height=20)

DotPlot(cd8_t_obj$RNA, features = features,group.by = 'Manually_curation') + RotatedAxis()
ggsave(file.path(outdir,paste0("T_transit_makers_dot.pdf")),  width=8, height=4)

FeaturePlot(cd8_t_obj$RNA, features = c('Havcr2','Tcf7','Cx3cr1','Tbx21','Gzmb','Mki67'),ncol = 3)
ggsave(file.path(outdir,paste0("T_transit_makers_feature.pdf")),  width=12, height=8)

VlnPlot(cd8_t_obj$RNA, features = features,group.by = 'Manually_curation')
ggsave(file.path(outdir,paste0("T_transit_makers_vln.pdf")),  width=12, height=12)

