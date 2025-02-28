library(ggplot2)
library(ggpubr)
library(reshape)

## ------- Inputs and parameters -------
mouse_survival_response=read.delim('../data/mouse_responder_survival.txt')
mouse_survival_non=read.delim('../data/mouse_non_responder_survival.txt')

outdir='../results/figure0'
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ------- Mouse survival curve of responders and non-responders ------
volumn_res=melt(mouse_survival_response,id.vars = c('Days'))
volumn_non=melt(mouse_survival_non,id.vars = c('Days'))

p=ggplot(data=volumn_res, aes(x=Days, y=value, group=variable)) +
  geom_line(color='#0073C2')+ylim(c(0,400))+theme_classic() +ylab('Tumor volumn')+
  geom_vline(xintercept = c(4,9,17,24), linetype = "dashed", color = "black")
p
ggsave(file.path(outdir,'responder_tumor_volumn.pdf'),width = 4,height = 3)

p=ggplot(data=volumn_non, aes(x=Days, y=value, group=variable)) +
  geom_line(color='#EFC000')+ylim(c(0,400))+theme_classic()+ylab('Tumor volumn')+
  geom_vline(xintercept = c(4,9,17,24), linetype = "dashed", color = "black")
p
ggsave(file.path(outdir,'nonresponder_tumor_volumn.pdf'),p,width = 4,height = 3)

## ------- Dataset summary -------
library(ggplot2)
library(dplyr)
library(reshape2)

# Load the data
df <- read.table("../data/datasets_summary.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE,row.names = 1)

# Extract numerical data
df=as.data.frame(t(df))
df$`Samples size`=as.numeric(df$`Samples size`)
df$responders=as.numeric(df$responders)
df$`non-responder`=as.numeric(df$`non-responder`)

df$Responder_ratio=df$responders/df$`Samples size`
df$non_Responder_ratio=df$`non-responder`/df$`Samples size`
df$cohort=rownames(df)

sample_order=c("Mouse sc","Mouse bulk","Allen et al. Pre","Luoma et al. PBMC Pre","Luoma et al. PBMC On/Post","Foy et al.Pre","INSPIRE Pre","Liu et al. On/Post","Liu et al. Pre","Obradovic et al. On/Post","Bill et al. Pre","Luoma et al. Tumor On/Post")
df$cohort=factor(df$cohort,sample_order)
# Melt for ggplot
num_data_melt <- melt(df[,c("Responder_ratio","non_Responder_ratio","cohort")], id.vars = "cohort")
num_data_melt$cohort=factor(num_data_melt$cohort,sample_order)
# Plot
ggplot(data=num_data_melt, aes(x=cohort, y=value, fill=variable)) +
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_brewer(palette="Set2")+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
ggsave(file.path(outdir,'responder_ratio.pdf'),width = 5,height = 3.5)

ggplot(data=df, aes(x=cohort, y=`Samples size`)) +
  geom_bar(stat="identity", fill="#8DA0CB")+
  theme_classic()+
  scale_fill_brewer(palette="Set2")+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggsave(file.path(outdir,'sample_size.pdf'),width = 4  ,height = 3.5)

library(ggplot2)
library(reshape2)

df_heat=df[,c('Availability','Speciman','Data type','sc_TCRseq','Timepoint','cohort')]
df_heat_melt=melt(df_heat,id.vars = c('cohort'))
df_heat_melt$value[which(df_heat_melt$value=='Singel cell')]='Single cell'

# Define a **professional & visually comfortable color palette**
colors_list <- c(
  # Specimen
  "PBMC" = "#2E86C1",       # Muted Blue
  "Tumor" = "#E74C3C",      # Soft Red
  
  # Availability
  "Public" = "#28B463",     # Calming Green
  "In-house" = "#8E44AD",    # Soft Purple
  
  # Data Type
  "Single cell" = "#F5B041", # Warm Orange
  "Bulk" = "#D35400",       # Earthy Brown
  
  # Timepoint (Distinguishing "Both", "Pre", and "Post")
  "Pre" = "#1ABC9C",        # Teal
  "Post" = "#5D6D7E",       # Soft Gray-Blue
  "Both" = "#9B59B6",       # Gentle Lavender
  
  # sc-TCRseq (Yes/No)
  "Yes" = "#FFC300",        # Soft Yellow-Gold
  "No" = "white"          # Dark Blue-Gray
)


# Plot with professional color scheme
ggplot(df_heat_melt, aes(x = cohort, y = variable, color = value)) +
  geom_point(size = 8, shape = 16) +  # Solid filled circles
  scale_color_manual(values = colors_list) +  # Apply the comfortable colors
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave(file.path(outdir,'heatmap.pdf'),width = 7,height = 6)
