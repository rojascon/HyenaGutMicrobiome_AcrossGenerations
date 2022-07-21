#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 20 Jul 2022
#                         Last updated:  20 Jul 2022
#
################################################################################

## CODE FOR: A) determine whether KEGG or COG abundances vary with host
#               individual identity, matriline, or prey abundance
#           B) construct PCoAs based on KEGG or COG abundances

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1. Load KEGG and COG gene counts 
#             (see script 02_tidy_genefunctions.R"               
################################################################################

load("data/metagenomes/02_KEGG_abundances.Rdata");
load("data/metagenomes/02_COG_abundances.Rdata");

# load sample metadata
load("data/metagenomes/00_metagenomes_metadata.Rdata");

# format gene data for vegdist
rownames(keggabun)=keggabun$gf; keggabun$gf=NULL;
rownames(cogabun)=cogabun$gf; cogabun$gf=NULL;

keggabun=keggabun[,order(colnames(keggabun))];
cogabun=cogabun[,order(colnames(cogabun))];

# transpose for vegdist
keggabun=as.data.frame(t(keggabun));
cogabun=as.data.frame(t(cogabun));


################################################################################
#             2. Generate Bray-Curtis distances and run PERMANOVA               
################################################################################
# get bray-curtis
kegg.dist=vegdist(keggabun, method="bray");
cog.dist=vegdist(cogabun, method="bray");

adonis2(kegg.dist ~    
          mat+
          hyenaID2+
          preym,
        data=meta, by="terms",
        method = met[i],     
        permutations = 999);


adonis2(cog.dist ~    
          mat+
          hyenaID2+
          preym,
        data=meta, by="terms",
        method = met[i],     
        permutations = 999);


################################################################################
#             3. PCoA of KEGG distances ~ hyena identity            
################################################################################

# calculate coordinates for PCoA
pcoa_dec=cmdscale(kegg.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,meta,by="sampleID"); 
pcoa_met$hyenaID2=factor(pcoa_met$hyenaID2,
                         levels=c("D1","G1","D3","G3"));

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

# PCoA-color coded by hyena identity
my_col=c("#253494","#7fcdbb","tan3","navajowhite");

pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=hyenaID2),
             size = 2.5,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Hyena ID")+
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

################################################################################
#             4. PCoA of COG distances ~ prey abundance         
################################################################################

# calculate coordinates for PCoA
pcoa_dec=cmdscale(cog.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,meta,by="sampleID"); 

# recode large prey densities >500 = 500 for plotting purposes
pcoa_met$preym[pcoa_met$preym>=500]=500;

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);


################################################################################
#             5. PCoA of COG distances ~ prey abundance         
################################################################################

pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=preym),
             size = 2.5,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Hyena ID")+
  scale_fill_viridis_c(begin=0.18)+
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
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa2);


################################################################################
#             6. save figures 
################################################################################

pcas=arrangeGrob(pcoa1,pcoa2, nrow=1); 

ggsave(filename="03_gene_PCoAs.pdf",
       device="pdf",path="./figures/metagenomes",
       plot=pcas,
       width=9.8,
       height=4,
       units="in",
       dpi=500);

