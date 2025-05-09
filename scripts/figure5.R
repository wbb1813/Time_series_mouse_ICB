library(ggplot2)
library(survival)
library(survminer)
library(dplyr)

outdir='../results/figure5'
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}
## ------- Inputs and parameters -------
bulk_patient_score=readRDS('../data/bulk_patient_score.rds')
tcga_hnsc=read.delim('../data/tcga_hnsc_mean_tb_score.txt')

df_comb_filter_pbmc=read.delim('../data/Teff_B_comb_AUC_PBMC.txt')
df_comb_filter_tumor=read.delim('../data/Teff_B_comb_AUC_tumor.txt')
df_compare_auc=read.delim('../data/auc_compare.txt')

bulk_fixed_cutoff_odds_sum=read.delim('../data/fix_cutoff_bulk_test.txt')
sc_fixed_cutoff_odds_sum=read.delim('../data/fix_cutoff_sc_test.txt')
#fixed_cutoff=readRDS('../data/fix_cutoff.rds')

## -------F5A: PBMC AUC barplot -------
## PBMC
df_comb_filter_pbmc$data_type='PBMC single-cell'
df_comb_filter_pbmc$sig_name=factor(df_comb_filter_pbmc$sig_name,levels = c('t2_Effector memory CD8','t2_B cell','mean_score'))
p=ggplot(data=df_comb_filter_pbmc, aes(x=id, y=AUC,fill=sig_name)) +
  geom_bar(stat="identity", position=position_dodge(),color='black')+
  geom_text(aes(label=AUC), position = position_dodge(0.9), vjust=-0.5, color="black", size=1.75)+
  theme_classic() + theme(axis.text.x = element_text(hjust = 1,angle = 60))+xlab('')+
  facet_grid(~data_type,scales = 'free_x',space = "free")+scale_fill_brewer(palette="Set2")

p
ggsave(file.path(outdir,'Teff_B_comb_filter_AUC_pbmc.pdf'),p,width = 5.5,height = 4.5)

## Tumor
df_comb_filter_tumor$data_type=ifelse(df_comb_filter_tumor$data_type=='Bulk','Tumor bulk','Tumor single-cell')
df_comb_filter_tumor$sig_name=factor(df_comb_filter_tumor$sig_name,levels = c('t2_Effector memory CD8','t2_B cell','mean_score'))
df_comb_filter_tumor=df_comb_filter_tumor[which(df_comb_filter_tumor$Cohort!='Aoki_2021'),] ## This is a esophageal cancer dataset

p=ggplot(data=df_comb_filter_tumor, aes(x=id, y=AUC,fill=sig_name)) +
  geom_bar(stat="identity", position=position_dodge(),color='black')+
  geom_text(aes(label=AUC), position = position_dodge(0.9), vjust=-0.5, color="black", size=1.75)+
  theme_classic() + theme(axis.text.x = element_text(hjust = 1,angle = 60))+xlab('')+
  facet_grid(~data_type,scales = 'free_x',space = "free")+scale_fill_brewer(palette="Set2")

p
ggsave(file.path(outdir,'Teff_B_comb_filter_AUC_tumor.pdf'),p,width = 10,height = 4.5)

## ------- Comparasion between combine score and other public signatuer score -------
# AUC
colors <- c("#d62728", "#ff7f0e", "#2ca02c","#1f77b4" , "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#aec7e8", "#ffbb78", "#98df8a")
sig_id=c('mean_score','teff_ifng','POPLAR_teff','proliferation','tgfb','stroma_emt','chemo_12','CD38','CXCL9','CD274','MEX3B')
df_compare_auc=df_compare_auc[which(df_compare_auc$id!='Luoma_tumor_sc_Combo_Post'),]
df_compare_auc=df_compare_auc[which(df_compare_auc$Cohort!='Aoki_2021'),]
df_compare_auc$sig_name=factor(df_compare_auc$sig_name,levels = sig_id)
p = ggplot(df_compare_auc, aes(x=sig_name, y=AUC)) +
  geom_boxplot(color="#1f77b4",alpha=0.3,outlier.shape = NA) +
  # Box plot with dot plot
  geom_jitter(aes(colour = id,shape=id), position=position_jitter(0.2),size=2)+
  theme_classic()+theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 45))+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=seq(0,15))+
  xlab('')+
  geom_hline(yintercept = 0.5,linetype='dashed',color="#d62728")
p

p = p + stat_compare_means(method = "wilcox.test",paired = T,ref.group = "mean_score",label='p.signif',method.args = list(alternative = "less")) # other groups compare to ref group, the alternative here should be "less"
p
ggsave(file.path(outdir,'auc_compare.pdf'),p,width = 7,height = 3.5)

## ------- Identify fix threshold ------
# # Plot the odds ratios against cutoffs
# pdf(file.path(outdir,'fixed_threshold_max_odds.pdf'),width = 4,height = 3)
# plot(fixed_cutoff$cutoffs, fixed_cutoff$odds_ratios, type = "l", xlab = "Cutoff", ylab = "Odds Ratio", main = "Maximizing Odds Ratio")
# abline(v = fixed_cutoff$best_cutoff, col = "red", lty = 2)
# dev.off()

## ------- Test fixed threshold in testing cohorts -------
## single cell
p=ggplot(sc_fixed_cutoff_odds_sum, aes(x = id, y = odds_ratio, fill = Type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Adds a black border to bars
  scale_fill_manual(values = c("Training" = "#D55E00", "Testing" = "#009E73")) +  # Dark orange and teal+
  geom_text(aes(label=signif(odds_ratio,2)), vjust=-0.3, size=3.5)+
  theme_classic() +
  labs(
    x = "Cohort",
    y = "Odds Ratio",
    fill = "Cohort Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave(file.path(outdir,'odds_ratio_train_test_sc.pdf'),p,width = 3.5,height = 3.5)
## Bulk
bulk_fixed_cutoff_odds_sum=bulk_fixed_cutoff_odds_sum[which(bulk_fixed_cutoff_odds_sum$cohort!='Aoki_2021'),]
p=ggplot(bulk_fixed_cutoff_odds_sum, aes(x = id, y = odds_ratio, fill = Type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Adds a black border to bars
  scale_fill_manual(values = c("Training" = "#D55E00", "Testing" = "#009E73")) +  # Dark orange and teal+
  geom_text(aes(label=signif(odds_ratio,2)), vjust=-0.3, size=3.5)+
  theme_classic() +
  labs(
    x = "Cohort",
    y = "Odds Ratio",
    fill = "Cohort Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
ggsave(file.path(outdir,'odds_ratio_train_test_bulk.pdf'),p,width = 3.5,height = 3.5)

## ------- Survival analysis -------
## Functions 
hr_plot=function(df,surv_time,surv_status,sex_info,fig_title,xlim=3.5,width = 3,height = 5,prefix){
  colnames(df)[which(colnames(df)==surv_time)]='time'
  colnames(df)[which(colnames(df)==surv_status)]='status'
  df$time=as.numeric(df$time)
  df$status=as.numeric(df$status)
  df$Age=as.numeric(df$Age)
  df$Sex=as.factor(df$Sex)
  
  if ('HPV' %in% colnames(df)){
    if (isTRUE(sex_info)){
      cox_model <- coxph(Surv(time, status) ~ Age + Sex + mean_score + HPV, data = df)
    }else{
      cox_model <- coxph(Surv(time, status) ~ Age + mean_score + HPV, data = df)
    }
    
  }else{
    if (isTRUE(sex_info)){
      cox_model <- coxph(Surv(time, status) ~ Age + Sex + mean_score, data = df)
    }else{
      cox_model <- coxph(Surv(time, status) ~ Age + mean_score, data = df)
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
    theme_minimal() +expand_limits(x = xlim)+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_text(aes(label = p_value_label, x = upper_ci), hjust = -0.3)  # P-value text
  
  ggsave(file.path(outdir,paste0(prefix,'.pdf')),p,width = width,height = height)
  return(hr_data)
}

## TCGA
df = tcga_hnsc
df$HPV[which(df$HPV=='negative')]='Neg'
df$HPV[which(df$HPV=='positive')]='Pos'
df=df[which(df$HPV!='indeterminate'),]
p_tcga=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=T,fig_title='TCGA HNSC',xlim = 15,width = 2,height = 4,prefix = 'tcga_hnsc_cox')


## Foy
df=bulk_patient_score$Foy1_2022
df$HPV_status[which(df$HPV_status==0)]='Neg'
df$HPV_status[which(df$HPV_status=='1')]='Pos'
df$HPV_status[which(df$HPV_status=='2')]='Unknown'
#df$HPV_status <- relevel(df$HPV_status, ref = "Unknown")
colnames(df)[which(colnames(df)=='HPV_status')]='HPV'

p_foy_os=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=T,fig_title='Foy et al. OS',xlim = 25,width = 2.3,height = 4,prefix = 'foy_os_cox')
p_foy_pfs=hr_plot(df = df,surv_time='PFS_days',surv_status='PFS_status',sex_info=T,fig_title='Foy et al. PFS',xlim = 25,width = 2.3,height = 4,prefix = 'foy_pfs_cox')

## INSPIRE
df=bulk_patient_score$INSPIRE
p_INSPIRE_os=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=F,fig_title='INSPIRE OS',xlim = 30,width = 2.3,height = 2,prefix ='INSPIRE_os_cox')
p_INSPIRE_pfs=hr_plot(df = df,surv_time='PFS_days',surv_status='PFS_status',sex_info=F,fig_title='INSPIRE PFS',xlim = 30,width = 2.3,height = 2,prefix ='INSPIRE_pfs_cox')

## Figure with merged datasets 
p_tcga$Cohort='TCGA HNSC OS'
p_foy_os$Cohort='Foy et al. OS'
p_foy_pfs$Cohort='Foy et al. PFS'
p_INSPIRE_os$Cohort='INSPIRE OS'
p_INSPIRE_pfs$Cohort='INSPIRE PFS'

df_sur=rbind(p_tcga,p_foy_os,p_foy_pfs,p_INSPIRE_os,p_INSPIRE_pfs)
df_sur$term[which(df_sur$term=='Sexmale')]='SexM'
df_sur$p_label='ns'
df_sur$p_label[which(df_sur$p_value_label<=0.05)]='*'
df_sur$p_label[which(df_sur$p_value_label<=0.01)]='**'
df_sur$p_label[which(df_sur$p_value_label<=0.001)]='***'
df_sur$Cohort=factor(df_sur$Cohort,c("TCGA HNSC OS","Foy et al. OS","Foy et al. PFS","INSPIRE OS","INSPIRE PFS"))
p=ggplot(df_sur, aes(x = hazard_ratio, y = term)) +
  geom_point(shape = 15) +  # Square point to match original
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dotted") +  # Reference line at HR = 1
  labs(
    title = '',
    x = "Hazard Ratio",
    y = NULL
  ) +
  theme_minimal() + expand_limits(x = 20)+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  geom_text(aes(label = p_label, x = upper_ci), hjust = -0.3) +  # P-value text
  facet_grid(~Cohort,scales='free_x')
p
ggsave(file.path(outdir,paste0('cox_hr','.pdf')),p,width = 6,height = 3)
## ------- KM curve -------
## Functions 
km_plot=function(df,surv_time='OS_days',surv_status='OS_status',prefix){
  df$score_group=ifelse(df$mean_score > mean(df$mean_score),'High','Low')
  
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

## TCGA
km_tcga_os=km_plot(df=tcga_hnsc,surv_time='time',surv_status='status',prefix='tcga_os')

km_tcga_neg_os=km_plot(df=tcga_hnsc[which(tcga_hnsc$HPV=='negative'),],surv_time='time',surv_status='status',prefix='tcga_os_neg')
km_tcga_pos_os=km_plot(df=tcga_hnsc[which(tcga_hnsc$HPV=='positive'),],surv_time='time',surv_status='status',prefix='tcga_os_pos')

## Foy
km_foy_os=km_plot(df=bulk_patient_score$Foy1_2022,surv_time='OS_days',surv_status='OS_status',prefix='foy_os')
km_foy_pfs=km_plot(df=bulk_patient_score$Foy1_2022,surv_time='PFS_days',surv_status='PFS_status',prefix='foy_pfs')

## Foy
km_INSPIRE_os=km_plot(df=bulk_patient_score$INSPIRE,surv_time='OS_days',surv_status='OS_status',prefix='INSPIRE_os')
km_INSPIRE_pfs=km_plot(df=bulk_patient_score$INSPIRE,surv_time='PFS_days',surv_status='PFS_status',prefix='INSPIRE_pfs')

