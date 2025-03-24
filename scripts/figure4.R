library('devtools')
# devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)

## ------- Inputs and parameters -------
outdir='../results/figure4'
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ------- heatmap anno with GO terms -------
hm_input=readRDS('../data/bulk_cluster_enrichment.rds')

## responder
anno_res=hm_input$res$anno_res
res_cm=hm_input$res$res_cm

ann_res_df=data.frame(id=c(rep('C1',3),rep('C3',3),rep('C6',4)),
                      term=c('mRNA processing','DNA replication','mitochondrial gene expression',#C1
                             'immunoglobulin production','B cell mediated immunity','protein-RNA complex assembly',#C3
                             'immune response','myeloid leukocyte activation','neutrophil activation','chemokine production'#C6
                      ))

pdf(file.path(outdir,'response_anno.pdf'),height = 10,width = 10)
visCluster(object = res_cm,
           plot.type = "both",
           column_names_rot = 45,
           annoTerm.data = ann_res_df,
           line.side = "left")
dev.off()

## Non-responder
non_cm=hm_input$non_res$non_cm

anno_non_df=data.frame(id=c(rep('C1',3),rep('C3',3),rep('C4',3),rep('C5',3),rep('C7',3),rep('C8',4)),
                       term=c('rRNA processing','protein-RNA complex assembly','tRNA processing',#C1
                              'activation of immune response','response to interferon-alpha','humoral immune response',#C3
                              'mRNA processing','RNA splicing','chromatin remodeling',#C4
                              'regulation of innate immune response','leukocyte mediated cytotoxicity','cell killing',#C5
                              'immunoglobulin production','B cell mediated immunity','lymphocyte mediated immunity',#C7
                              'epithelial cell differentiation','response to wounding','cell chemotaxis','cell-cell junction organization'
                              
                       ))

pdf(file.path(outdir,'nonresponse_anno.pdf'),height = 10,width = 10)
visCluster(object = non_cm,
           plot.type = "both",
           column_names_rot = 45,
           annoTerm.data = anno_non_df,
           line.side = "left")
dev.off()

## ------- Barplot of ML AUC -------
single_tp_res=data.frame(timepoint=c('1','2','3','4'),AUC=c(0.62,0.79,0.7,0.79))
diff_tp_res=data.frame(timepoint=c('T2-T1','T3-T1','T3-T2','T4-T1','T4-T2','T4-T3'),AUC=c(0.52,0.55,0.74,0.69,0.82,0.69))

## single timepoint
p = ggplot(data=single_tp_res, aes(x=timepoint, y=AUC)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.9),width = 0.8,color="black",fill='#66c2a5')+
  geom_text(aes(label=signif(AUC,2)), vjust=-1, color="black",
            position = position_dodge(0.9), size=2.5)+
  #scale_color_brewer(palette="Dark2")+
  theme_classic()
p
ggsave(file.path(outdir,paste0('single_timepoint_auc.pdf')),p,width = 3.5,height = 3)

## diff timepoint
p = ggplot(data=diff_tp_res, aes(x=timepoint, y=AUC)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.9),width = 0.8,color="black",fill='#fc8d62')+
  geom_text(aes(label=signif(AUC,2)), vjust=-1, color="black",
            position = position_dodge(0.9), size=2.5)+
  #scale_color_brewer(palette="Dark2")+
  theme_classic()
p
ggsave(file.path(outdir,paste0('diff_timepoint_auc.pdf')),p,width = 5,height = 3)
