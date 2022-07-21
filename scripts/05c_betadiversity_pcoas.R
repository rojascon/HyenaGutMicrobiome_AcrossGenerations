#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 10 Feb 2021
#                         Last updated: 27 May 2021
#
################################################################################

## CODE FOR: generating PCoA plots based on distances color-coded by:
#     A) matriline
#     B) age in years
#     C) prey abundance
#     D) calendar year

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load dissimilarity matrices generated in script 05a
#                         & formatted metadata table                 
################################################################################

load("data/02_sample_metadata_formatted.Rdata");
load("data/05_dissimilarity_distances.Rdata");

# remove samples without prey availability from metadata
meta2=metadf[metadf$preym>3,];

# rename matrilines (e.g. 1 to "M1"
meta2$mat=as.character(meta2$mat);
meta2$mat2=paste("M",meta2$mat,sep="");
meta2$mat2<-factor(meta2$mat2, levels=c("M1","M2","M3","M4"));


################################################################################
#             2. Make data frame of PCoA coordinates 
################################################################################

# calculate coordinates for PCoA
pcoa_dec=cmdscale(wuni.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,meta2,by="sampleID"); 

# calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);


################################################################################
#             2. Plot PCoA ordinations
################################################################################

# PCoA-color coded by matriline
my_col=c("#74c476","dodgerblue3","#fdcdac","mediumorchid3");

pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=mat2),
             size = 2.5,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Matriline")+
  scale_fill_manual(values=my_col)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa1);

# PCoA-color coded by age
pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=age_yrs),
             size = 2.5,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Age (yrs)")+
  scale_fill_gradient(low = "bisque", high = "#774411")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill = guide_colorbar(reverse = TRUE));plot(pcoa2);

# PCoA-color coded by prey density
pcoa3=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=preym),
             size = 2.5,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Prey Avail.")+
  scale_fill_gradient(low = "#ccece6", high = "#253494")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill = guide_colorbar(reverse = TRUE));plot(pcoa3);

# PCoA-color coded by sample year
pcoa4=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=year),
             size = 2.5,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Year")+
  scale_fill_gradient(low = "black", high = "#f0f0f0")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill = guide_colorbar(reverse = TRUE));plot(pcoa4);


################################################################################
#             3. Save plots
################################################################################
mypcas=arrangeGrob(pcoa1, pcoa2,
                    pcoa3,pcoa4, nrow=2); 

ggsave(filename="05_betadiv_pcoas.pdf",
       device="pdf",path="./figures",
       plot=mypcas,
       width=10,
       height=7,
       units="in",
       dpi=500);









