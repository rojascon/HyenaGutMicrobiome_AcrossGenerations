#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 3 Feb 2021
#                         Last updated: 20 June 2022
#
################################################################################

## CODE FOR: a) identifying the core bacterial genera found in >85% of samples
#   b) making heatmaps of their abundances across samples 
#   c) testing whether the abundances of core genera varied with 4 host predictors: 
#         matriline, age, prey density, and calendar year using LMMs
#   d) plotting the beta-coefficients of statistically significant predictors
#       as determined by the linear models above

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load filtered ASV table, ASV taxonomy, and 
#                       formatted sample metadata
################################################################################

load("data/01_ASV_table_filtered.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T); 
load("data/02_sample_metadata_formatted.Rdata");


################################################################################
#             2.  Identify the core bacterial genera found in 85% of samples
################################################################################
# add genus labels to ASV names
tax2=tax[,c(1,7)];
rownames(tax2)=tax2$ASV; tax2$ASV=NULL;
asv_tax=merge(asvf,tax2,by="row.names"); 
asv_tax=asv_tax[,-1];

# collapse ASV table by genera
gd=aggregate(.~Genus, asv_tax, sum); 
rownames(gd)=gd$Genus;
gd$Genus=NULL;

# convert abundance data to presence/absence data
gpad=(gd>0)*1; 
gpad=as.data.frame(gpad);

# identify the core genera found in 85% of samples
corenum=round(0.85*ncol(gpad));
gcor=gpad[rowSums(gpad)>corenum,];
print(paste("The number of core genera was:", nrow(gcor)));
print("The taxonomic labels of these core genera were:"); print(rownames(gcor));


################################################################################
#             3.  Make heatmaps showing the relative abundances 
#         of core genera with mean relative abundances of at least 0.5%
################################################################################

# obtain abundance data - GENERA
gd2=t(gd);
gd2=apply(gd2, 1, function(i) (i/sum(i))*100);
gd2=gd2[rownames(gd2) %in% rownames(gcor),];
gd2=gd2[rowMeans(gd2) > 0.5,]; gd2=data.frame(gd2);
gd2$tnames=rownames(gd2);

# melt data frame for plotting with ggplot2
gen_abun<-reshape2::melt(gd2, id.vars="tnames",value.name = "abun");
colnames(gen_abun)[2]="sampleID";
gen_abun=merge(gen_abun, metadf, by="sampleID");
gen_abun$abun2=gen_abun$abun;gen_abun$abun2[gen_abun$abun2>=50]=50;

# plot heatmap!
heat1=ggplot(data=gen_abun, aes(reorder(sampleID,ideal_order), 
                                y=tnames, fill=abun2)) + 
  geom_tile() +
  facet_grid(~hyenaID2,scales="free_x")+ 
  scale_fill_viridis_c(begin=0.18) +
  labs(y="",
       x="",
       fill="Relative Abundance (%)")+
  theme_bw()+
  theme(legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=13, face="bold"))+
  theme(strip.text.x = element_text(size = 13)); 
plot(heat1); 


################################################################################
#         4.  Test whether core genera abundances vary with 4 host predictors:
#                 matriline, age, prey density, and calendar year
################################################################################

# retain core genera with >1% mean relative abundance
gd2$tnames=NULL; 
gtest=gd2[rowMeans(gd2[-1]) > 1,]; 

# add sample metadata to genera abundances
rownames(metadf)=metadf$sampleID; metadf$sampleID=NULL;
gtest=data.frame(t(gtest)); 
gtest=merge(gtest, metadf, by="row.names"); 
colnames(gtest)[1]="sampleID";

# remove samples without prey availability info
gtest=gtest[complete.cases(gtest$preym),];

# Do bacterial genera vary with host predictors? 
# export model output to an empty data frame called gres
gres=c(); 
for(i in 2:18)
{
  print(paste("Linear mixed model for core genera:", colnames(gtest)[i]));
  gmod=lmer(gtest[,i]~
              gtest$mat+
              gtest$age_yrs+
              gtest$preym+
              gtest$year+
              (1|gtest$hyenaID2),
            data=gtest, na.action=na.omit);
  coefs <- data.frame(coef(summary(gmod)));
  coefs$df.Satt <- coef(summary(gmod))[, 3];
  coefs$p.Satt <- coef(summary(gmod))[, 5];
  tempdf=as.data.frame(coefs);
  tempdf$taxa=colnames(gtest)[i]; tempdf$terms=rownames(tempdf);
  rownames(tempdf)=NULL;
  gres=rbind(gres,tempdf);
};

# trim model output to only include Beta-coefficients and their pvalues
gmodtbl=gres[gres$terms!="(Intercept)",];
gmodtbl$padj=p.adjust(gmodtbl$p.Satt, n=nrow(gmodtbl),method="fdr");
#gmodtbl2=gmodtbl[gmodtbl$p.Satt < 0.05,];
gmodtbl2=gmodtbl[gmodtbl$padj < 0.05,];
#
#write.csv(gmodtlb, "model.output.coregenera.csv", row.names=F);


################################################################################
#           5.  Plot beta-coefficients of model terms with p.values a<0.05
################################################################################

# create color palette for each model variable
phy_col=c('#fc8d62','#7570b3');

# plot diverging dot plot
betagen<-ggplot(gmodtbl2, aes(x=reorder(taxa,Estimate),
                          y=Estimate)) +
  geom_point(aes(fill=terms),
             stat='identity',
             pch=21,
             size = 3)+
  scale_fill_manual(values=phy_col,
                    labels=c("PreyAbund","Year"))+
  coord_flip() +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "black", size=0.6)+
  labs(y="Estimate",
       x="", 
       fill="") +
  theme(legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=9, face="bold"),
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"),
        axis.text.x=element_text(size=10),
        strip.text = element_text(size =8, face="bold"))+
  guides(fill = guide_legend(nrow=1));plot(betagen);


################################################################################
#           6.  export heatmap and diverging plot
################################################################################

ggsave(filename="04_heatmap_coregenera.pdf",
       device="pdf",path="./figures",
       plot=heat1,
       width=9,
       height=5.5, 
       units="in",
       dpi=500);

ggsave(filename="04_coregenera_LMMbetacoeff.pdf",
       device="pdf",path="figures",
       plot=betagen,
       width=7,
       height=2.7,
       units="in",
       dpi=500);

