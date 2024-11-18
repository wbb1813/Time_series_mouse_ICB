library(ggplot2)
library(survival)
library(survminer)
library(dplyr)
outdir='../results/figure4'
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}
## ------- Inputs and parameters -------
bulk_patient_score=readRDS('/data/Binbin_Kun/binbin/silvio/results/ICB_odds_fix_cutoff/bulk_patient_score.rds')
tcga_hnsc=read.delim('/data/Binbin_Kun/binbin/silvio/results/ICB_response_signature_genes/t2_response_tb/survival/TCGA_HNSC/tcga_hnsc_mean_tb_score.txt')

df_comb_filter=read.delim('../data/auc_teff_b_cell_patients.txt')
fixed_cutoff_odds_sum=read.delim('../data/fixed_cutoff_odds_sum.txt')
fixed_cutoff=readRDS('../data/fix_cutoff.rds')

## -------F4B: AUC barplot -------=
p=ggplot(data=df_comb_filter, aes(x=id, y=AUC,fill=sig_name)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=AUC), position = position_dodge(0.9), vjust=1.6, color="white", size=2.5)+
  theme_classic() + theme(axis.text.x = element_text(hjust = 1,angle = 60))+xlab('')+
  facet_grid(~data_type,scales = 'free_x',space = "free")+scale_fill_brewer(palette="Dark2")

p
ggsave(file.path(outdir,'Teff_B_comb_filter_AUC.pdf'),p,width = 12,height = 5)

## ------- Identify fix threshold ------
# Plot the odds ratios against cutoffs
pdf(file.path(outdir,'fixed_threshold_max_odds.pdf'),width = 4,height = 3)
plot(fixed_cutoff$cutoffs, fixed_cutoff$odds_ratios, type = "l", xlab = "Cutoff", ylab = "Odds Ratio", main = "Maximizing Odds Ratio")
abline(v = fixed_cutoff$best_cutoff, col = "red", lty = 2)
dev.off()

## ------- Test fixed threshold in testing cohorts -------
p=ggplot(fixed_cutoff_odds_sum, aes(x = id, y = odds_ratio, fill = Type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Adds a black border to bars
  scale_fill_manual(values = c("Training" = "#D55E00", "Testing" = "#009E73")) +  # Dark orange and teal+
  geom_text(aes(label=signif(odds_ratio,2)), vjust=-0.3, size=3.5)+
  theme_minimal() +
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
ggsave(file.path(outdir,'odds_ratio_train_test.pdf'),p,width = 4,height = 3.5)

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
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_text(aes(label = p_value_label, x = upper_ci), hjust = -0.3)  # P-value text
  
  return(p)
  ggsave(file.path(outdir,paste0(prefix,'.pdf')),p,width = 4,height = 6)
}

## TCGA
df = tcga_hnsc
df$HPV[which(df$HPV=='negative')]='Neg'
df$HPV[which(df$HPV=='positive')]='Pos'
df=df[which(df$HPV!='indeterminate'),]
p_tcga=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=T,fig_title='TCGA HNSC','tcga_hnsc_cox')

## Foy
df=bulk_patient_score$Foy1_2022
df$HPV_status[which(df$HPV_status==0)]='Neg'
df$HPV_status[which(df$HPV_status=='1')]='Pos'
df$HPV_status[which(df$HPV_status=='2')]='Unknown'
#df$HPV_status <- relevel(df$HPV_status, ref = "Unknown")
colnames(df)[which(colnames(df)=='HPV_status')]='HPV'

p_foy_os=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=T,fig_title='Foy et al. OS','foy_os_cox')
p_foy_pfs=hr_plot(df = df,surv_time='PFS_days',surv_status='PFS_status',sex_info=T,fig_title='Foy et al. PFS','foy_pfs_cox')

## INSPIRE
df=bulk_patient_score$INSPIRE
p_INSPIRE_os=hr_plot(df = df,surv_time='OS_days',surv_status='OS_status',sex_info=F,fig_title='INSPIRE OS','INSPIRE_os_cox')
p_INSPIRE_pfs=hr_plot(df = df,surv_time='PFS_days',surv_status='PFS_status',sex_info=F,fig_title='INSPIRE PFS','INSPIRE_pfs_cox')

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
  return(p_km)
  ggsave(file.path(outdir,paste(prefix,'_km_curve.pdf',sep = "_")),p_km$plot,width = 4,height = 3)
  
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

