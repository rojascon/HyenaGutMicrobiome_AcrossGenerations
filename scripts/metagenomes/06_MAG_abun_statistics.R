#################################################################################
#
#               Gut Microbiome Across 3 Generations of wild spotted hyenas
#                      
#              Rojas et al 2021.Gut Microbiome stability and function across
#                   three generations of wild spotted hyenas
#
#                               By: Connie Rojas
#                               Created: 20 Mar 2022
#                            Last updated: 21 July 2022
#
################################################################################

#           "METAGENOME ASSEMBLED GENOMES - MAGS"
##CODE FOR: A) PERMANOVA statistics on MAG abundances across samples
#           B) PCoA ordinations based on MAG abundances

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#         1. Load MAG abundances as calculated by CoverM software 
#                         and sample metadata
################################################################################

coverm=read.table("data/metagenomes/00_CoverM_MAG_abun.tsv", sep="\t", header=T);
load("data/metagenomes/00_metagenomes_metadata.Rdata");

# tidy sample names
coverm=coverm[coverm$Genome!="unmapped",];
rownames(coverm)=coverm$Genome; coverm$Genome=NULL;
colnames(coverm)=gsub("_notmapped.1.fq.Relative.Abundance....", "", colnames(coverm));


################################################################################
#           2. Generate Bray-Curtis distances based on 
#           MAG abundances across samples and run PERMANOVA
################################################################################

#transpose taxa abundance table so taxa are columns
coverm=coverm[,order(colnames(coverm))];
covermt=t(coverm); 

# BRAY-CURTIS distance
bray.dist=vegdist(covermt, method="bray");

# run PERMANOVA
adonis2(bray.dist~   
         preym+
          mat+
         hyenaID2+
         year,
       data=meta, by="term",
       method = met[i],     
       permutations = 999);


################################################################################
#           3. PERMANOVA based on abundances of MAGs classified to Species
#           all other MAGs filtered from dataset
################################################################################

# load MAG taxonomy data
mag=read.csv("data/metagenomes/00_metabat_MAGs.csv", header=T, row.names=1); 

# retain MAGs classified to species
mag[mag==""] <- NA
mag2=mag[complete.cases(mag$Species),];
covermt2=covermt[,colnames(covermt) %in% rownames(mag2)];

# BRAY-CURTIS distance
bray.dist2=vegdist(covermt2, method="bray");
adonis2(bray.dist2~
          preym+
          mat+
         hyenaID2+
          year,
       data=meta, by="term",
       method = met[i],     
       permutations = 999);


################################################################################
#           4. PCoAs - all MAGs
################################################################################
# extract PCoA coordinates
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,meta,by="sampleID"); 
pcoa_met$hyenaID2=factor(pcoa_met$hyenaID2,
                         levels=c("D1","G1","D3","G3"));

# calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

# set color palette
hyena <- c("#253494","#7fcdbb", "tan3","navajowhite");

# plot the PCoA-color coded by hyenaID
pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=hyenaID2),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Hyena ID",
       title="PCoA of MAG Abundances")+
  scale_fill_manual(values=hyena)+
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

# save plot
ggsave(filename="06_MAG_PCoA.pdf",
       device="pdf",path="./figures/metagenomes",
       plot=pcoa1,
       width=5.6,
       height=4,
       units="in",
       dpi=500);


################################################################################
#           4. PCoAs - MAGs classified to Species
################################################################################

# extract PCoA coordinates
pcoa_dec=cmdscale(bray.dist2, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,meta,by="sampleID"); 
pcoa_met$hyenaID2=factor(pcoa_met$hyenaID2,
                         levels=c("D1","G1","D3","G3"));

# calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

# plot the PCoA-color coded by hyenaID
pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=hyenaID2),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Hyena ID",
       title="PCoA of Species-level MAG Abundances")+
  scale_fill_manual(values=hyena)+
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

# save plot
ggsave(filename="06_MAG_PCoA_species.pdf",
       device="pdf",path="./figures/metagenomes",
       plot=pcoa2,
       width=5.6,
       height=4,
       units="in",
       dpi=500);




