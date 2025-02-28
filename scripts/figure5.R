library(cutpointr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(plyr)
library(ModelMetrics)
library(glmnet)
library(cutpointr)
library(survival)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(pROC)
library(scales)
library(rstatix)
library(abind)
library(circlize)
library(survival)
library(survminer)

## ------- Inputs and parameters -------
inter_score=readRDS('../data/patient_interaction_score.rds')

outdir='../results/figure5'

if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}
## -------  Calculate AUC R vs. NR across each timepoint RDI/RAI in mouse sc data -------
auc_input_rdi=readRDS('../data/mouse_RDI_AUC.rds')
auc_input_rdi_234=auc_input_rdi[which(auc_input_rdi$feature!='merge_feature'),]
ggplot(data=auc_input_rdi_234, aes(x=timepoint, y=auc, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), color='black')+
  geom_text(aes(label=round(auc,2)), vjust=-0.3, size=3.5, position = position_dodge(0.9))+ 
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  theme_pubr() +
  theme(legend.position = "right")+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  facet_wrap(~feature, scales='fixed',nrow = 1)+
  ylab('AUC')

ggsave(file.path(outdir,'mouse_sc_auc_rdi_timepoint.pdf'),width = 7,height = 3.5)

## ------- AUC in patient cohorts -------
df_auc=read.delim('../data/RDI_AUC_patients.txt')
p=ggplot(data=df_auc, aes(x=Group, y=AUC)) +
  geom_bar(stat="identity", fill="#88BECD",color='black')+
  geom_text(aes(label=signif(AUC,2)), vjust=-0.5, color="black", size=3.5)+
  theme_classic() + theme(axis.text.x = element_text(hjust = 1,angle = 60))+xlab('')
p
ggsave(file.path(outdir,paste0('mouse_IRIS_interaction_bulk_patient_AUC.pdf')),p,width = 5,height = 4.5)

## ------- Interaction plot -------
## Function 
plot_interactions <- function(DD, cell_color_mapping = NULL, annotations = NULL, annotations_color_map = NULL, cell_separation_angle = 20, annotation_track_height = 4, annotation_track_margin_mm = 3) {
  Netw <- data.frame(Source = sapply(DD[,1], function(x) {
    aa <- strsplit(x,split = "_")[[1]]
    return(paste(aa[1],":",aa[3],sep = ""))
  }),
  Target = sapply(DD[,1], function(x) {
    aa <- strsplit(x,split = "_")[[1]]
    return(paste(aa[2],":",aa[4],sep = ""))
  }),
  weights = sapply(DD[,2], function(x) {
    return(x)
  }),
  stringsAsFactors = F
  )
  cc <- sapply(union(Netw$Source, Netw$Target), function(x) strsplit(x, split = ":")[[1]][1])
  dd <- data.frame(from =  Netw$Source, to = Netw$Target, value =  Netw$weights, stringsAsFactors = F)
  grid.col = rep("grey", length(cc))
  if(is.null(cell_color_mapping)){
    set.seed(1234) 
    cl <- colors(distinct = TRUE)
    mycl = sample(cl, length(unique(cc)))
    for(i in which(names(cc) %in% Netw$Source)){
      grid.col[i] <- mycl[which(unique(cc) == strsplit(names(cc)[i],split = ":")[[1]][1])]
    }
    names(grid.col) <- names(cc)
  }
  
  else{
    mycl = cell_color_mapping
    for(i in which(names(cc) %in% Netw$Source)){
      grid.col[i] <- mycl[strsplit(names(cc)[i],split = ":")[[1]][1]]
    }
    names(grid.col) <- names(cc)
  }
  
  chordDiagram(dd,preAllocateTracks = list(list(
    track.height = mm_h(4),
    track.margin = c(mm_h(1), 0)
  ),list(
    track.height = mm_h(annotation_track_height),
    track.margin = c(mm_h(annotation_track_margin_mm), 0)
  ), list(
    track.height = mm_h(4),
    track.margin = c(mm_h(3), mm_h(3))
  )),annotationTrack=NULL, 
  grid.col = grid.col, 
  order = names(cc), 
  group = cc, 
  big.gap = cell_separation_angle, 
  directional = 1, 
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow")
  circos.trackPlotRegion(track.index = 3, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    lr = strsplit(CELL_META$sector.index,split = ":")[[1]][2]
    circos.text(mean(xlim),mean(ylim),lr, col = "black", cex = 0.5, facing = 'reverse.clockwise', niceFacing = TRUE)
  }, bg.border = 0)
  
  if(is.null(cell_color_mapping)){
    for(ct in unique(cc)){
      highlight.sector(names(cc)[cc == ct], track.index = 1, col = mycl[which(unique(cc) == ct)], 
                       text = ct, cex = 0.8, text.col = "black", niceFacing = TRUE, text.vjust = -1.5)
    }
  }
  
  else{
    for(ct in unique(cc)){
      highlight.sector(names(cc)[cc == ct], track.index = 1, col = mycl[ct], 
                       text = ct, cex = 0.8, text.col = "black", niceFacing = TRUE, text.vjust = -1.5)
    }
  }
  
  
  if(!is.null(annotations)){
    if(is.null(annotations_color_map))
      stop("annotations_color_map missing!")
    for(aa in unique(annotations)){
      ss <- Netw$Target[annotations == aa]
      highlight.sector(ss, track.index = 2, col = annotations_color_map[aa], 
                       text = "", cex = 0.8, text.col = "black", niceFacing = TRUE, text.vjust = -1.5)
    }
  }
  
  
  circos.clear()
  
}


## ------- Plot interaction -------
com_intr=readRDS('../data/patient_RDI.rds')

## Make CELL LR DF
Feature_df = com_intr %>% as.data.frame(.) %>% set_colnames('interaction') %>%
  data.frame(.,str_split_fixed(.$interaction, '\\_', 4)) %>%
  set_colnames(c('interaction','Lcell','Rcell','L','R')) %>% mutate(GenePair=paste(L,R,sep='_')) #%>% mutate(annotation=plyr::mapvalues(GenePair, from=LR_anno$GenePair, to=LR_anno$annotations))

## ------- Figure 3A -------
DD = Feature_df %>% dplyr::select(names=1) %>% mutate(weights=0.25, annotations = 'UNK') #DELETE
DD$annotations = ifelse(is.na(DD$annotations), 'UNK', DD$annotations)
int_ann = DD$annotations
names(int_ann) <- DD$names
acm <- c("white","darkorange","aquamarine2","darkorchid","deeppink","deepskyblue")
names(acm) <- c('UNK', 'chemotaxis', 'cell-adhesion/LTEM', 'activating/co-stimulatory', 'pro-inflammatory', 'checkpoint/inhibitory')

## save figure
pdf(file.path(outdir,'interactions.pdf'),width = 6,height = 6)
plot_interactions(DD[,1:2], cell_color_mapping = NULL, annotations = DD$annotations, annotations_color_map = acm, cell_separation_angle = 2, annotation_track_height = 2, annotation_track_margin_mm = 4)
dev.off()

## ------- Survival analysis -------
## Functions 
hr_plot=function(df,surv_time,surv_status,sex_info,fig_title,prefix){
  colnames(df)[which(colnames(df)==surv_time)]='time'
  colnames(df)[which(colnames(df)==surv_status)]='status'
  df$time=as.numeric(df$time)
  df$status=as.numeric(df$status)
  df$Age=as.numeric(df$Age)
  df$Sex=as.factor(df$Sex)
  
  if ('HPV' %in% colnames(df)){
    if (isTRUE(sex_info)){
      cox_model <- coxph(Surv(time, status) ~ Age + Sex + score + HPV, data = df)
    }else{
      cox_model <- coxph(Surv(time, status) ~ Age + score + HPV, data = df)
    }
    
  }else{
    if (isTRUE(sex_info)){
      cox_model <- coxph(Surv(time, status) ~ Age + Sex + score, data = df)
    }else{
      cox_model <- coxph(Surv(time, status) ~ Age + score, data = df)
    }
  }
  
  
  # Extract coefficients and confidence intervals
  cox_summary <- summary(cox_model)
  hr_data <- data.frame(
    term = rownames(cox_summary$coefficients),
    hazard_ratio = (cox_summary$coefficients[, "coef"]),
    lower_ci = log(cox_summary$conf.int[, "lower .95"]),
    upper_ci = log(cox_summary$conf.int[, "upper .95"]),
    p_value = cox_summary$coefficients[,'Pr(>|z|)']
  )
  
  hr_data <- hr_data %>% mutate(
    p_value_label = ifelse(!is.na(p_value) & p_value < 0.05, 
                           sprintf("%.3f", p_value), sprintf("%.3f", p_value))
  )
  
  # Plot
  p=ggplot(hr_data, aes(x = hazard_ratio, y = term)) +
    geom_point(shape = 15) +  # Square point to match original
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dotted") +  # Reference line at HR = 1
    labs(
      title = fig_title,
      x = "Hazard Ratio",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_text(aes(label = p_value_label, x = upper_ci), hjust = -0.3)  # P-value text
  
  ggsave(file.path(outdir,paste0(prefix,'.pdf')),p,width = 4,height = 6)
  return(p)
}

## Foy
Foy_score=inter_score$Foy_score
df=Foy_score
df$HPV_status[which(df$HPV_status==0)]='Neg'
df$HPV_status[which(df$HPV_status=='1')]='Pos'
df$HPV_status[which(df$HPV_status=='2')]='Unknown'
colnames(df)[which(colnames(df)=='HPV_status')]='HPV'

p_foy_os=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=T,fig_title='Foy et al. OS','foy_os_cox')
p_foy_pfs=hr_plot(df = df,surv_time='PFS_days',surv_status='PFS_status',sex_info=T,fig_title='Foy et al. PFS','foy_pfs_cox')

## INSPIRE
INSPIRE_score=inter_score$INSPIRE_score
df=INSPIRE_score
p_INSPIRE_os=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=F,fig_title='INSPIRE OS','INSPIRE_os_cox')
p_INSPIRE_pfs=hr_plot(df = df,surv_time='PFS_days',surv_status='PFS_status',sex_info=F,fig_title='INSPIRE PFS','INSPIRE_pfs_cox')

## ------- KM curve -------
## Functions 
km_plot=function(df,surv_time='OS_days',surv_status='OS_status',prefix){
  df$score_group=ifelse(df$score > mean(df$score),'High','Low')
  
  colnames(df)[which(colnames(df)==surv_time)]='time'
  colnames(df)[which(colnames(df)==surv_status)]='status'
  df$time=as.numeric(df$time)
  df$status=as.numeric(df$status)
  
  # Fit the Kaplan-Meier model
  km_fit <- survfit(Surv(time, status) ~ score_group, data = df)
  
  # Plot the Kaplan-Meier curves
  p_km=ggsurvplot(km_fit, data = df, pval = TRUE, conf.int = F, 
                  risk.table = TRUE, ggtheme = theme_minimal(),palette = c("darkred", "#2E9FDF"))
  
  ggsave(file.path(outdir,paste(prefix,'_km_curve.pdf',sep = "_")),p_km$plot,width = 4,height = 3)
  return(p_km)
}


## Foy
km_foy_os=km_plot(df=Foy_score,surv_time='OS_days',surv_status='OS_status',prefix='foy_os')
km_foy_pfs=km_plot(df=Foy_score,surv_time='PFS_days',surv_status='PFS_status',prefix='foy_pfs')

## Foy
km_INSPIRE_os=km_plot(df=INSPIRE_score,surv_time='OS_days',surv_status='OS_status',prefix='INSPIRE_os')
km_INSPIRE_pfs=km_plot(df=INSPIRE_score,surv_time='PFS_days',surv_status='PFS_status',prefix='INSPIRE_pfs')




