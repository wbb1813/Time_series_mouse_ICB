library(ggplot2)
library(Seurat)
library(ggpubr)
library(trend)
library(Kendall)
library(reshape)

## ------- Inputs and parameters -------
#allcell_abj=readRDS('/vf/users/Binbin_Kun/binbin/silvio/results/seurat/res0.6/annotated_seurat_res0.6.rds')
allcell_abj=readRDS('/Volumes/Binbin_Kun/binbin/silvio/results/seurat/res0.6/annotated_seurat_res0.6.rds')
cell_abund=read.delim('../data/cell_abundance_timepoint.txt')
cell_abund_sample=read.delim('../data/cell_abundance_sample_timepoint.txt')
subtype_abund_sample=read.delim('../data/subtype_abund_sample_timepoint.txt')

outdir='../results/figure2'

if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ------- F1A: UMAP of all cells -------
p1 = DimPlot(allcell_abj$RNA, reduction = "umap", label = F, repel = F, group.by = 'merged_tcell',
            cols = c('Erythroblast' = '#D95F02', "CD8 T" = '#E7298A',"CD4 T" = '#B22222', "B cell"  = '#7570B3',
                     "Monocyte"  = '#1B9E77',"NK cell"  = '#66A61E',"Neutrophil"  = '#E6AB02',
                     "Dendritic cell"  = '#A6761D',"Basophil"  = '#666666'),raster=FALSE) + ggtitle("")
p1
ggsave(file.path(outdir,'f1a_all_cell_umap.png'), p1,  width=6, height=4)

## ------- F1B: abundance of major cell types across timepoints (merge samples in each timepoint) -------
cols = c('Erythroblast' = '#D95F02', "CD8 T" = '#E7298A',"CD4 T" = '#B22222', "B cell"  = '#7570B3',
         "Monocyte"  = '#1B9E77',"NK cell"  = '#66A61E',"Neutrophil"  = '#E6AB02',
         "Dendritic cell"  = '#A6761D',"Basophil"  = '#666666')

p2 = ggplot(data=cell_abund, aes(x=timepoint, y=Relative_frequency, fill=cell_type)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = cols)+
  theme_classic()+ylab('Relative frequency')+xlab('Timepoint')+
  theme(legend.title=element_blank())
p2

ggsave(file.path(outdir,'f1b_all_cell_abundace_timepoint.pdf'),p2,width = 3.5,height = 4)

## ------- F1C: box plot of CD8, CD4, B cell, and Neutrophil cell abundance at different time points -------
## Functions 
box_line_plot=function(tmp_df,alternative='less',prefix){
  my_comparisons <- list( c("1", "2"),c("2", "3"),c("3", "4"),c("1", "3"),c("1", "4"))
  p=ggplot(tmp_df, aes(x = timepoint, y = Relative_frequency)) + 
    # geom_boxplot(aes(fill = Group)) +
    geom_boxplot(aes(fill = timepoint),alpha=0.5, width=0.5) +
    geom_line(aes(group = Patient),color='gray45') + 
    #geom_point(size = 2, aes(color = Group))+ 
    geom_point(size = 2, color = 'black')+ 
    theme_classic()+
    theme(legend.position = "none")+
    #stat_compare_means(comparisons = my_comparisons,paired = TRUE)
    #stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test',method.args = list(alternative = "greater"))+
    #stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')+ 
    xlab('Timepoints')+ylab('Relative frequency')
  if (alternative=='less'){
    p = p + stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test',method.args = list(alternative = "less"),label = "p.signif")
  } else if (alternative=='greater'){
    p = p + stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test',method.args = list(alternative = "greater"),label = "p.signif")
  }
  ggsave(file.path(outdir,paste0(prefix,'.pdf')),p,width = 3.5,height = 4.5)
  return(p)
}

# Mann-Kendall test for monotonic trend 
one_tail_mk=function(df){
  test_result <- Kendall(x=as.numeric(as.character(df$timepoint)), y=df$Relative_frequency)
  
  # Extract p-value and tau
  tau <- test_result$tau
  p_value <- test_result$sl
  
  # Adjust p-value for one-tailed test
  if (tau > 0) {
    p_one_tailed <- p_value / 2  # For testing an increasing trend
  } else {
    p_one_tailed <- 1 - (p_value / 2)  # For testing a decreasing trend
  }
  
  return(p_one_tailed)
}


## Figures
cell_abund_sample=na.omit(cell_abund_sample)
cell_abund_sample$timepoint=factor(cell_abund_sample$timepoint,levels = c(1,2,3,4))

# CD8 cell 
tmp_df=cell_abund_sample[which(cell_abund_sample$cell_type=='CD8 T'),]
tmp_df=na.omit(tmp_df)
p3_1=box_line_plot(tmp_df,alternative='less',prefix='CD8_cell_all_sample')

# CD4 cell 
tmp_df=cell_abund_sample[which(cell_abund_sample$cell_type=='CD4 T'),]
tmp_df=na.omit(tmp_df)
p3_2=box_line_plot(tmp_df,alternative='less',prefix='CD4_cell_all_sample')

# B cell 
tmp_df=cell_abund_sample[which(cell_abund_sample$cell_type=='B cell'),]
tmp_df=na.omit(tmp_df)
p3_3=box_line_plot(tmp_df,alternative='less',prefix='b_cell_all_sample')

# Neutrophil cell 
tmp_df=cell_abund_sample[which(cell_abund_sample$cell_type=='Neutrophil'),]
tmp_df=na.omit(tmp_df)
p3_4=box_line_plot(tmp_df,alternative='greater',prefix='Neutrophil_all_sample')

## ------- F1D: box plot of Teff, Th1 cell abundance between responders and non-responders -------
## Functions 
# Box plot
diff_abundance_box=function(cond1,cond2,prefix,alternative = "less"){
  df=rbind(cond1,cond2)
  df$Response=factor(df$Response,levels = c('Responders','Non-responders'))
  p <- ggboxplot(df, x = "timepoint", y = "Relative_frequency",
                 color = "Response", palette = "jco",
                 add = "jitter")
  # p = p + stat_compare_means(aes(group = Response),label = "p.format",method = 'wilcox.test',method.args = list(alternative = "less"))
  p = p + stat_compare_means(aes(group = Response),method = 'wilcox.test',method.args = list(alternative = alternative),label = "p.signif")
  
  # p = p + stat_compare_means(aes(group = Response),label = "p.format",method = 'wilcox.test')
  p
  ggsave(file.path(outdir,paste0(prefix,'.pdf')),p,width = 6,height = 4)
  return(p)
}

# Line plot
diff_abundance_line=function(cond1,cond2,prefix,alternative = "less"){
  df=rbind(cond1,cond2)
  df$Response=factor(df$Response,levels = c('Responders','Non-responders'))
  
  p=ggline(df, x = "timepoint", y = "Relative_frequency", add = "mean_se",
         color = "Response", palette = "jco")+
    #stat_compare_means(aes(group = Response), label = "p.signif")
    stat_compare_means(aes(group = Response),method = 'wilcox.test',method.args = list(alternative = alternative),label = "p.signif")
  
  ggsave(file.path(outdir,paste0(prefix,'.pdf')),p,width = 6,height = 4)
  return(p)
}


## Teff abundance 
subtype_abund_sample=na.omit(subtype_abund_sample)
subtype_abund_sample$timepoint=factor(subtype_abund_sample$timepoint,levels = c(1,2,3,4))

# response
df_res=subtype_abund_sample[which(subtype_abund_sample$Response=='Responders'&subtype_abund_sample$cell_type=='Effector memory CD8'),]
df_res=na.omit(df_res)
p4_1=box_line_plot(df_res,alternative='less',prefix='teff_response')

mk_res=one_tail_mk(df_res)

p4_1=p4_1+geom_text(x=3.5, y=0.0535, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('teff_response_mk','.pdf')),p4_1,width = 3.5,height = 4.5)

# non-response
df_non=subtype_abund_sample[which(subtype_abund_sample$Response=='Non-responders'&subtype_abund_sample$cell_type=='Effector memory CD8'),]
df_non=na.omit(df_non)
p4_2=box_line_plot(df_non,alternative='less',prefix='teff_non_response')

mk_res=one_tail_mk(df_non)

p4_2=p4_2+geom_text(x=3.5, y=0.052, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('teff_non_response_mk','.pdf')),p4_2,width = 3.5,height = 4.5)

# responder vs non-responder
p4_3=diff_abundance_box(cond1=df_res,cond2=df_non,prefix='teff_res_non',alternative = "less")
p4_4=diff_abundance_line(cond1=df_res,cond2=df_non,prefix='teff_res_non_line',alternative = "less")

## Th1 CD4 abundance
# response
df_res=subtype_abund_sample[which(subtype_abund_sample$Response=='Responders'&subtype_abund_sample$cell_type=='Th1'),]
df_res=na.omit(df_res)
sp1=box_line_plot(df_res,alternative='less',prefix='th1_response')

mk_res=one_tail_mk(df_res)

sp1=sp1+geom_text(x=3.5, y=0.0535, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('th1_response_mk','.pdf')),sp1,width = 3.5,height = 4.5)

# non-response
df_non=subtype_abund_sample[which(subtype_abund_sample$Response=='Non-responders'&subtype_abund_sample$cell_type=='Th1'),]
df_non=na.omit(df_non)
sp1_non=box_line_plot(df_non,alternative='less',prefix='th1_non_response')

mk_res=one_tail_mk(df_non)

sp1_non=sp1_non+geom_text(x=3.5, y=0.052, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('th1_non_response_mk','.pdf')),sp1_non,width = 3.5,height = 4.5)

# responder vs non-responder
sp1_3=diff_abundance_box(cond1=df_res,cond2=df_non,prefix='th1_res_non',alternative = "less")
sp1_4=diff_abundance_line(cond1=df_res,cond2=df_non,prefix='th1_res_non_line',alternative = "less")

## Treg abundance
# response
df_res=subtype_abund_sample[which(subtype_abund_sample$Response=='Responders'&subtype_abund_sample$cell_type=='Treg'),]
df_res=na.omit(df_res)
sp_treg=box_line_plot(df_res,alternative='less',prefix='treg_response')

mk_res=one_tail_mk(df_res)

sp_treg=sp_treg+geom_text(x=3.5, y=0.02, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('treg_response_mk','.pdf')),sp_treg,width = 3.5,height = 4.5)

# non-response
df_non=subtype_abund_sample[which(subtype_abund_sample$Response=='Non-responders'&subtype_abund_sample$cell_type=='Treg'),]
df_non=na.omit(df_non)
sp_treg_non=box_line_plot(df_non,alternative='less',prefix='treg_non_response')

mk_res=one_tail_mk(df_non)

sp_treg_non=sp_treg_non+geom_text(x=3.5, y=0.02, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('treg_non_response_mk','.pdf')),sp_treg_non,width = 3.5,height = 4.5)

# responder vs non-responder
sp_treg_3=diff_abundance_box(cond1=df_res,cond2=df_non,prefix='treg_res_non',alternative = "less")
sp_treg_4=diff_abundance_line(cond1=df_res,cond2=df_non,prefix='treg_res_non_line',alternative = "less")

## B cell abundance 
# response
df_res=subtype_abund_sample[which(subtype_abund_sample$Response=='Responders'&subtype_abund_sample$cell_type=='B cell'),]
df_res=na.omit(df_res)
p5_1=box_line_plot(df_res,alternative='less',prefix='b_response')

mk_res=one_tail_mk(df_res)
p5_1=p5_1+geom_text(x=3.5, y=0.052, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('b_response_mk','.pdf')),p5_1,width = 3.5,height = 4.5)

# non-response
df_non=subtype_abund_sample[which(subtype_abund_sample$Response=='Non-responders'&subtype_abund_sample$cell_type=='B cell'),]
df_non=na.omit(df_non)
p5_2=box_line_plot(df_non,alternative='less',prefix='b_non_response')

mk_res=one_tail_mk(df_non)
p5_2=p5_2+geom_text(x=3.5, y=0.052, label=paste0("Mann-Kendall test p-value: ", signif(mk_res[1],2)))
ggsave(file.path(outdir,paste0('b_non_response_mk','.pdf')),p5_2,width = 3.5,height = 4.5)

# responder vs non-responder
p5_3=diff_abundance_box(cond1=df_res,cond2=df_non,prefix='b_res_non',alternative = "less")
p5_4=diff_abundance_line(cond1=df_res,cond2=df_non,prefix='b_res_non_line',alternative = "less")

# ## ------- Arrange figure 1 -------
# library("cowplot")
# p=ggdraw() +
#   draw_plot(p1, x = 0, y = .7, width = .3, height = .25) +
#   draw_plot(p2, x = 0.3, y = .7, width = .2, height = .25) +
#   draw_plot(p3_1, x = 0.5, y = .7, width = .125, height = 0.2) +
#   draw_plot(p3_2, x = 0.625, y = .7, width = .125, height = 0.2) +
#   draw_plot(p3_3, x = 0.75, y = .7, width = .125, height = 0.2) +
#   draw_plot(p3_4, x = 0.875, y = .7, width = .125, height = 0.2) +
#   draw_plot(p4_1, x = 0, y = .4, width = .125, height = 0.2) +
#   draw_plot(p4_2, x = 0.2, y = .4, width = .125, height = 0.2) +
#   draw_plot(p4_3, x = 0.4, y = .4, width = .6, height = 0.25) +
#   draw_plot(p5_1, x = 0, y = .1, width = .125, height = 0.2) +
#   draw_plot(p5_2, x = 0.2, y = .1, width = .125, height = 0.2) +
#   draw_plot(p5_3, x = 0.4, y = .1, width = .6, height = 0.25) +
#   draw_plot_label(label = c("A", "B", "C",'D','E','F','G','H','I'), size = 15,
#                   x = c(0, 0.25, 0.4,0,0.2,0.4,0,0.2,0.4), y = c(0.98, 0.98, 0.98,0.65,0.65,0.65,0.35,0.35,0.35))
# p
# ggsave('figure1.pdf',p,width = 17,height = 23)
# 
# f1_1=ggarrange(p1, p2, p3_1, p3_2, p3_3, p3_4 + rremove("x.text"), 
#           labels = c("A", "B", "C"),widths=c(2,1.5,1,1,1,1),
#           ncol = 6, nrow = 1)
# f1_2=ggarrange(p4_1, p4_2, p4_3 + rremove("x.text"), 
#                labels = c("D", "E", "F"),widths=c(2,2,6),
#                ncol = 3, nrow = 1)
# f1_3=ggarrange(p5_1, p5_2, p5_3 + rremove("x.text"), 
#                labels = c("G", "H", "I"),widths=c(2,2,6),
#                ncol = 3, nrow = 1)
# p = ggarrange(f1_1,f1_2,f1_3,ncol = 1, nrow = 3)
# p
# ggsave('figure1.pdf',p,width = 17,height = 23)
# 
# 
# 
