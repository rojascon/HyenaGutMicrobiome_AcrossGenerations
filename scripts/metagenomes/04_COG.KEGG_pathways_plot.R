#################################################################################
#
#               Gut Microbiome Across 3 Generations of wild spotted hyenas
#                      
#              Rojas et al 2021.Gut Microbiome stability and function across
#                   three generations of wild spotted hyenas
#
#                               By: Connie Rojas
#                               Created: 20 Nov 2021
#                            Last updated: 21 July 2022
#
################################################################################

## CODE FOR: generating stacked barplots of KEGG and COG pathway abundances
              
source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load KEGG and COG data and sample metadata                
################################################################################
load("data/metagenomes/02_KEGG_abundances.Rdata");
load("data/metagenomes/02_COG_abundances.Rdata");
load("data/metagenomes/00_metagenomes_metadata.Rdata");


################################################################################
#             2. Barplot of KEGG abundances             
################################################################################

# remove unchar protein and putative transposases for plotting purposes
keggabun$avg=rowMeans(keggabun[,-1]);
top=keggabun[keggabun$avg<800,]; 
top=top[top$avg>150,]; top$avg=NULL;

# recode large and small subunit proteins
top$gf[16:23]="large subunit ribosomal protein";
top$gf[44:49]="small subunit ribosomal protein";
top=aggregate(.~gf,top, sum);
top=top[rowMeans(top[-1])>200,];

# melt data frame for ggplot
kbar<-reshape2::melt(top, id.vars="gf",value.name = "abun");
colnames(kbar)[1:2]=c("protein","sampleID")
kbar=merge(kbar, meta, by="sampleID");
kbar$hyenaID2=factor(kbar$hyenaID2,levels=c("D1","G1","D3","G3"))

# set color palette
c25 <- c(
  "dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00",
  "black", "gold1","skyblue2", "#FB9A99","palegreen2",
  "#CAB2D6", "#FDBF6F","gray70", "khaki2","maroon", "orchid1", 
  "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", 
  "yellow3","darkorange4", "brown");

# plot
bark=ggplot(data=kbar, 
                aes(reorder(sampleID,ideal_order),y=abun, fill=protein))+
  geom_bar(stat = "identity")+
  facet_grid(~hyenaID2,scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "TPM",
       fill="KEGG Protein")+
  scale_fill_manual(values=c25)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text.x = element_text(size = 13,face="bold"))+
  guides(fill = guide_legend(ncol=1));plot(bark)

# save plot
ggsave(filename="04_kegg_bar.pdf",
       device="pdf",path="figures/metagenomes",
       plot=bark,
       width=11,
       height=7.5,
       units="in",
       dpi=500);


################################################################################
#             3. Barplot of COG abundances             
################################################################################

# convert COG TPM abundances to proportions so all samples add up to 100
cogabun[,-1] <- lapply(cogabun[,-1], function(x) (x/sum(x))*100);

# remove less abundant pathways for plotting purposes
cogabun$avg=rowMeans(cogabun[,-1]);
top=cogabun[cogabun$avg>2,]; top$avg=NULL;
colnames(top)[1]="COG.Pathway";

# denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(top[2:ncol(top)])); 
top=rbind(top, newrow); 
top$COG.Pathway=as.character(top$COG.Pathway);
top[nrow(top),1]="Other";

# melt data frame for ggplot
cbar<-reshape2::melt(top, id.vars="COG.Pathway",
                     value.name = "TPM");
colnames(cbar)[2]="sampleID";
cbar=merge(cbar, meta, by="sampleID");
cbar$hyenaID2=factor(cbar$hyenaID2,levels=c("D1","G1","D3","G3"));

# set color palette
cogcol=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
       "#117777", "#44AAAA", "#737373","#117744","#88CCAA", "#777711", 
       "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", 
       "#AA4455", "#DD7788","black");
cogcol2=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#737373", "#44AA77", "#88CCAA", "#777711", "#AAAA44", 
          "#DDDD77","#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "beige",
          "grey","black");

# create plot
barc=ggplot(data=cbar, 
               aes(x=reorder(sampleID,ideal_order),y=TPM, fill=COG.Pathway))+ 
  geom_bar(stat = "identity")+
  facet_grid(~hyenaID2,scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="COG Category")+ #COG Pathway COG Category
  scale_fill_manual(values=cogcol)+ #fam_col; bcol
  theme(legend.position="right", 
        legend.text = element_text(size=10),
        legend.title = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13,face="bold"))+
  guides(fill=guide_legend(ncol=1));plot(barc);

# save plot
ggsave(filename="04_cog_bar.pdf",
       device="pdf",path="./figures/metagenomes",
       plot=barc,
       width=9.5,  
       height=6,
       units="in",
       dpi=500);