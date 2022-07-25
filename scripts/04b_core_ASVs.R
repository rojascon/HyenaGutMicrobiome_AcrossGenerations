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

## CODE FOR: a) identifying the core bacterial ASVs found in >85% of samples
#   b) making heatmaps of their abundances across samples 
#   c) testing whether the abundances of core ASVs varied with 4 host predictors: 
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
#             2.  Identify the core bacterial ASVs found in 85% of samples
################################################################################
# add asv taxonomy label to ASV names 
tax3=tax[,c(1,8)];
rownames(tax3)=tax3$ASV; tax3$ASV=NULL;
asv_tax=merge(asvf,tax3,by="row.names"); 
asv_tax=asv_tax[,-1];
rownames(asv_tax)=asv_tax$asv_genus;
asv_tax$asv_genus=NULL;

# convert abundance data to presence/absence data
apad=(asv_tax>0)*1; 
apad=as.data.frame(apad);

# identify the core ASVs found in 85% of samples
corenum=round(0.85*ncol(apad));
acor=apad[rowSums(apad)>corenum,];
print(paste("The number of core ASVs was:", nrow(acor)));
print("The taxonomic labels of these core ASVs were:"); print(rownames(acor));


################################################################################
#             3.  Make heatmaps showing the relative abundances 
#                         of core ASVs
################################################################################

# obtain abundance data - ASVs
ad2=t(asv_tax);
ad2=apply(ad2, 1, function(i) (i/sum(i))*100);
ad2=ad2[rownames(ad2) %in% rownames(acor),];
ad2=ad2[rowMeans(ad2) > 0.5,]; ad2=data.frame(ad2);
ad2$tnames=rownames(ad2);

# melt data frame for plotting with ggplot2
asv_abun<-reshape2::melt(ad2, id.vars="tnames",value.name = "abun");
colnames(asv_abun)[2]="sampleID";
asv_abun=merge(asv_abun, metadf, by="sampleID");
asv_abun$abun2=asv_abun$abun;asv_abun$abun2[asv_abun$abun2>=50]=50;

# plot heatmap!
heat1=ggplot(data=asv_abun, aes(reorder(sampleID,ideal_order), 
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
#         4.  Test whether core ASV abundances vary with 4 host predictors:
#                 matriline, age, prey density, and calendar year
################################################################################

# retain core ASVs with >1% mean relative abundances
ad2$tnames=NULL;
atest=ad2[rowMeans(ad2[-1]) > 1,]; 

# add sample metadata to ASV abundances
rownames(metadf)=metadf$sampleID; metadf$sampleID=NULL;
atest=data.frame(t(atest));
atest=merge(atest, metadf, by="row.names");
colnames(atest)[1]="sampleID";

# remove samples without prey availability info
atest=atest[complete.cases(atest$preym),];

# Do bacterial ASVs vary with host predictors? 
# export model output to an empty data frame called ares
ares=c(); 
for(i in 2:11)
{
  print(paste("Linear mixed model for core ASVs:", colnames(atest)[i]));
  amod=lmer(atest[,i]~
              atest$mat+
              atest$age_yrs+
              atest$preym+
              atest$year+
              (1|atest$hyenaID2),
            data=atest, na.action=na.omit);
  coefs <- data.frame(coef(summary(amod)));
  coefs$df.Satt <- coef(summary(amod))[, 3];
  coefs$p.Satt <- coef(summary(amod))[, 5];
  tempdf=as.data.frame(coefs);
  tempdf$taxa=colnames(atest)[i]; tempdf$terms=rownames(tempdf);
  rownames(tempdf)=NULL;
  ares=rbind(ares,tempdf);
};

# trim model output to only include Beta-coefficients and their pvalues
amodtbl=ares[ares$terms!="(Intercept)",];
amodtbl$padj=p.adjust(amodtbl$p.Satt, n=nrow(amodtbl),method="fdr");
#amodtbl2=amodtbl[amodtbl$p.Satt<0.05,];
amodtbl2=amodtbl[amodtbl$padj<0.05,];
#write.csv(amodtlb, "model.output.coreASVs.csv", row.names=F);


################################################################################
#           5.  Plot beta-coefficients of model terms with p.values a<0.05
################################################################################
# create color palette for each model variable
phy_col=c('#fc8d62','#7570b3');

# plot diverging dot plot
betaasv<-ggplot(amodtbl2, aes(x=reorder(taxa,Estimate),
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
             color = "black", size=0.7)+
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
  guides(fill = guide_legend(nrow=1));plot(betaasv);


################################################################################
#           6.  export heatmap and diverging plot
################################################################################

ggsave(filename="04_heatmap_coreASVs.pdf",
       device="pdf",path="./figures",
       plot=heat1,
       width=9,
       height=5.5, 
       units="in",
       dpi=500);

ggsave(filename="04_coreASVs_LMMbetacoeff.pdf",
       device="pdf",path="figures",
       plot=betaasv,
       width=7,
       height=2.3,
       units="in",
       dpi=500);

