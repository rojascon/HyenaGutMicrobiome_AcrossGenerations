#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 21 Jan 2021
#                         Last updated: 15 March 2021
#
################################################################################

## CODE FOR: generating stacked barplots showing the relative abundances of 
#  bacterial phyla, orders, and genera

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load filtered ASV table, ASV taxonomy, and 
#                       formatted sample metadata
################################################################################

load("data/01_ASV_table_filtered.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T); #row.names=1
load("data/02_sample_metadata_formatted.Rdata");

# attach taxonomy to the ASV table 
rownames(tax)=tax$ASV; tax$ASV=NULL;
asv_tax=merge(asvf,tax,by="row.names"); 
colnames(asv_tax)[1]="ASV";

# make a vector of your sample names for later
samples=metadf$sampleID;


################################################################################
#             2. Create Phylum level barplots                 
################################################################################

# select bacterial taxonomic rank
phylum=asv_tax[,which(names(asv_tax) 
                      %in% c(samples, "Phylum"))];
colnames(phylum)[ncol(phylum)]="taxa";

# calculate taxa relative abundances 
phylum=aggregate(.~taxa, phylum, sum);  
phylum[,-1] <- lapply(phylum[,-1], function(x) (x/sum(x))*100);
#print(colSums(phylum[-1]));

# keep phyla >1% relative abundance across samples
phylum$AVG=rowMeans(phylum[,-1]);
phylum=phylum[phylum$AVG>1,];
phylum$AVG=NULL;

# denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(phylum[2:ncol(phylum)])); 
phylum=rbind(phylum, newrow); 
phylum$taxa=as.character(phylum$taxa);
phylum[nrow(phylum),1]="Other";

# melt data frame for ggplot
pbar<-reshape2::melt(phylum, id.vars="taxa",value.name = "abun");
colnames(pbar)[2]="sampleID";
pbar=merge(pbar, metadf, by="sampleID");

# set color-palette
phy_col=c('#66c2a5','#fc8d62','#7570b3','#e78ac3','#a6d854','#ffd92f');

# create plot; 
# organize samples by date (earliest sample to most recent)
barphy=ggplot(data=pbar, 
              aes(reorder(sampleID,sample_date),y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Phylum")+
  scale_fill_manual(values=phy_col)+
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
        axis.title.x = element_text(size=14, face="bold"));

plot(barphy);

# save image 
ggsave(filename="03_barplot_phylum.pdf",
       device="pdf",path="./figures",
       plot=barphy,
       width=8,
       height=5,
       units="in",
       dpi=500);


################################################################################
#             3. Create Order level barplots                 
################################################################################

# select bacterial taxonomic rank 
order=asv_tax[,which(names(asv_tax) 
                   %in% c(samples, "Order"))];
colnames(order)[ncol(order)]="taxa";

# calculate taxa relative abundances 
order=aggregate(.~taxa, order, sum);  
order[,-1] <- lapply(order[,-1], function(x) (x/sum(x))*100);
#print(colSums(order[-1]));

# keep orders >0.6% relative abundance across samples
order$AVG=rowMeans(order[,-1]);
order=order[order$AVG>0.6,];
order$AVG=NULL;

# denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(order[2:ncol(order)])); 
order=rbind(order, newrow); 
order$taxa=as.character(order$taxa);
order[nrow(order),1]="Other";

# melt data frame for ggplot
cbar<-reshape2::melt(order, id.vars="taxa",value.name = "abun");
colnames(cbar)[2]="sampleID";
cbar=merge(cbar, metadf, by="sampleID");

# set color-palette
order_col=c("#8c6bb1","#6baed6","#238b45","#fc9272",
          "#525252","#e7298a","#08589e","#bf812d","#a50f15");

# create plot
barorder=ggplot(data=cbar, 
              aes(reorder(sampleID,sample_date),y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Order")+
  scale_fill_manual(values=order_col)+
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
        axis.title.x = element_text(size=14, face="bold"));

plot(barorder);

# save image 
ggsave(filename="03_barplot_order.pdf",
       device="pdf",path="./figures",
       plot=barorder,
       width=8,
       height=5,
       units="in",
       dpi=500);


################################################################################
#             4. Create Genus level barplots                 
################################################################################

# select bacterial taxonomic rank 
gen=asv_tax[,which(names(asv_tax) 
                   %in% c(samples, "Genus"))];
colnames(gen)[ncol(gen)]="taxa";

# calculate taxa relative abundances 
gen=aggregate(.~taxa, gen, sum);  
gen[,-1] <- lapply(gen[,-1], function(x) (x/sum(x))*100);
#print(colSums(gen[-1]));

# keep genera >3% relative abundance across samples
gen$AVG=rowMeans(gen[,-1]);
gen=gen[gen$AVG>3,];
gen$AVG=NULL;

# denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(gen[2:ncol(gen)])); 
gen=rbind(gen, newrow); 
gen$taxa=as.character(gen$taxa);
gen[nrow(gen),1]="Other";

# melt data frame for ggplot
gbar<-reshape2::melt(gen, id.vars="taxa",value.name = "abun");
colnames(gbar)[2]="sampleID";
gbar=merge(gbar, metadf, by="sampleID");

# set color-palette
gen_col=c("#771155", "#f768a1", "#CAB2D6", "#114477", "#77AADD", 
          "#774411", "#AA7744","#DDAA77", 
          "#016c59","#44AA77",
          "gray70","grey35"); 

# create plot
bargen=ggplot(data=gbar, 
                aes(reorder(sampleID,sample_date),y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Genus")+
  scale_fill_manual(values=gen_col)+
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
        axis.title.x = element_text(size=14, face="bold"));

plot(bargen);

# save image 
ggsave(filename="03_barplot_genus.pdf",
       device="pdf",path="./figures",
       plot=bargen,
       width=9,
       height=5.5,
       units="in",
       dpi=500);



