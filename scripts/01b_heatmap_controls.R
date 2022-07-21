#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 21 Jan 2021
#                         Last updated: 18 Jul 2022
#
################################################################################

## CODE FOR: creating heatmap of ASVs found at high abundances in extraction kit 
#  control samples 

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load ASV Table of control samples generated 
#                 in script 01a_remove_contaminants.R               
################################################################################

load("data/01_controls_ASV_table.Rdata");

# load ASV taxonomy file
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T);  

# attach ASV genus label to the ASV table 
tax=tax[,c(1,8)];
rownames(tax)=tax$ASV; tax$ASV=NULL;
asv_tax=merge(controls,tax,by="row.names"); 


################################################################################
#             2.  prepare ASV table for plotting with pheatmap                
################################################################################

# clean up ASV taxonomy table
asv_tax$Row.names=NULL;
rownames(asv_tax)=asv_tax$asv_genus; 
asv_tax$asv_genus=NULL;
asv_tax=as.matrix(asv_tax);

# calculate ASV relative abundances 
asv_heat=prop.table(asv_tax,2);
print(colSums(asv_heat));

# keep ASVs >0.7% relative abundance across samples
asv_heat=as.data.frame(asv_heat);
asv_heat$AVG=rowMeans(asv_heat);
asv_heat=asv_heat[asv_heat$AVG>0.007,];
asv_heat=asv_heat[order(-asv_heat$AVG),];
asv_heat$AVG=NULL; 
paste("the number of ASVs that will show on heatmap is:", nrow(asv_heat), sep=" ");


################################################################################
#             3.  Plot heatmap of the most abundant ASVs                
################################################################################

# select color-palette to indicate ASV relative abundances
heat_color<-colorRampPalette(c("cornsilk","maroon4"), 
                             space = "rgb")(nrow(asv_heat));

# plot heatmap
pheatmap(
  mat               = asv_heat,
  color             = heat_color,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 9,
  main              = ""
);

# save manually by first manually adjusting R Plots' window to desired size
# then going to Plots > Export > Save as PDF
# save file in directory "figures", with filename= 01_controls_ASVheatmap.pdf
# 6.67 x 5.66

