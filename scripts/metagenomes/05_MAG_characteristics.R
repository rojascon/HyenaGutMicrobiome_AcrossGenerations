#################################################################################
#
#               Gut Microbiome Across 3 Generations of wild spotted hyenas
#                      
#              Rojas et al 2021.Gut Microbiome stability and function across
#                   three generations of wild spotted hyenas
#
#                               By: Connie Rojas
#                               Created: 19 Feb 2022
#                            Last updated: 21 July 2022
#
################################################################################

#           "METAGENOME ASSEMBLED GENOMES - MAGS"
##CODE FOR: making plots of 
#           A) MAG contamination x MAG completeness
#           B) MAG Sizes (Mb)
#           C) Annotation Rate (by GTDB-Tk)
#           C) Taxonomic classification

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load file with list of 149 MAGs                 
################################################################################

mag=read.csv("data/metagenomes/00_metabat_MAGs.csv", header=T, row.names=1); 


################################################################################
#             2. scatter of MAG completeness x contamination
################################################################################

magp1=ggplot(mag,aes(x=Completeness, y=Contamination))+
  geom_point(size=2, color="#969696")+ 
  geom_smooth(method=lm, color="#88419d", se=TRUE)+
  labs(y="MAG Contamination",
       x="MAG Completeness",
       title="A)")+
  theme_classic() + 
  theme(axis.title = element_text(size = 12, face="bold"), 
        axis.text = element_text(size = 11),
        plot.title=element_text(size = 12, face="bold"),
        legend.position="none"); plot(magp1);


################################################################################
#             3. histogram of MAG Sizes
################################################################################

magp2=ggplot() + aes(mag$Size_Mb)+
  geom_histogram(bins=8,             
                 color="black",
                 fill="#decbe4") +
  scale_x_continuous(breaks=seq(0.56,3.93,by=0.8))+
  theme_classic() +
  labs(x = "MAG Size (Kb)",
       y = "Number of MAGs",
       title="B)")+
  theme(axis.title = element_text(size = 12,face="bold"), 
        axis.text = element_text(size = 11),
        plot.title=element_text(size = 12, face="bold"));plot(magp2);


################################################################################
#         4. Barplot of annotation rate (% assigned to a phylum, a class, etc)
################################################################################

# make data frame of # of MAGs assigned to Family, # of MAGs assigned to species, etc
categs=c("Phylum","Class","Order","Family","Genus","Species");
nums=c(149,149,149,148,101,31);
pers=(nums/149) * 100;
andf=data.frame(categs,nums,pers);
andf$categs=factor(andf$categs, levels=c("Phylum","Class","Order",
                                   "Family","Genus","Species"));

# set color palette
andf$categs=forcats::fct_rev(factor(andf$categs));
new_col=c("#bebada","#fdcdac");

# plot!
magp3=ggplot(data=andf, 
            aes(x=pers,y=categs))+
  geom_bar(stat = "identity",position="dodge",fill="#fdcdac",colour="black")+
  theme_bw()+ 
  labs(x = "Annotation Rate (%)",
       y = "",
       title="C)")+
  theme(legend.position="none", 
        plot.title=element_text(size = 12, face="bold"),
        panel.background = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        axis.title.x = element_text(size=11, face="bold"),
        axis.text.x = element_text(size=12));plot(magp3);


################################################################################
#         5. Pie of MAG Family abundances
################################################################################

# showcase MAG families with counts of 2 or more
famdf=as.data.frame(table(mag$Family),stringsAsFactors = FALSE);
num.other=149-sum(famdf$Freq[famdf$Freq>2]);

# add a new row to capture "Other" category
famdf$Var1[1]="Other";
famdf$Freq[famdf$Var1=="Other"]=num.other;
famdf=famdf[famdf$Freq>2,];

# convert counts to proportions
famdf$Freq2=(famdf$Freq/sum(famdf$Freq))*100;
famdf$sample="S1";

# select color-palette
piecol=c("maroon", "darkturquoise", "#225ea8", "gold1", "#FF7F00","yellow3",
         "brown", "gray70", "green4", "#FB9A99", "darkorange4", "deeppink1",
         "#6A3D9A", "steelblue4", "skyblue2", "yellow4", "palegreen2", 
         "black", "khaki2", "#CAB2D6", "#bf812d", "#FDBF6F");

# plot piechart in ggplot2
magp4=ggplot(famdf, aes(x=sample, y=Freq2, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=piecol)+
  theme_void()+
  labs(fill="",
       x="",
       title="C) MAG Families")+
  theme(axis.text.x=element_blank(),
        text=element_text(size=13, face="bold"),
        legend.position = "bottom",
        legend.text = element_text(size=10),
        axis.title.y=element_text(size=12, face="bold"))+
  guides(fill = guide_legend(ncol=3));plot(magp4);


################################################################################
#         6. Save figures
################################################################################
sFigs=arrangeGrob(magp1,magp2,magp3, nrow=2); 

ggsave(filename="05_MAG_characteristics.pdf",
       device="pdf",path="./figures/metagenomes",
       plot=sFigs,
       width=7,
       height=5.5,
       units="in",
       dpi=500);

ggsave(filename="05_MAG_composition.pdf",
       device="pdf",path="./figures/metagenomes",
       plot=magp4,
       width=7,
       height=8,
       units="in",
       dpi=500);