library(ClusterGVis)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)

## ------- Inputs and parameters -------
outdir='../results/figure3'
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
anno_non=read.delim(file.path(outdir,'non','cluster_all_pathway_res.tsv'),sep = ' ')

anno_non0=anno_non[which(anno_non$cluster==8&anno_non$p.adjust<=0.01),]

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
