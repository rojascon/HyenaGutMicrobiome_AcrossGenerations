#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 21 Jan 2021
#                         Last updated: 19 Jul 2022
#
################################################################################

## CODE FOR: A) generating rarefaction curves of ASV richness after subsampling
#            B) histogram showing number of reads / sample

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Generate rarefaction curve data using mothur software                 
################################################################################

# Step 0: subsample ASV table to 2900 reads / sample (see script 06a)

# Step1: open subsampled data file "data/06_output_mothur.txt" in Excel 
# and save as .txt. For some reason, mothur insists on this step.

# Step2: open mothur and run command:
#       rarefaction.single(shared=05_output_mothur.txt, calc=sobs, freq=2)
# code will take ~20 mins to run

# Step3: rename output file as "07_mothur_rarefaction.txt" and place in the 
#  data directory


################################################################################
#             2. Clean and arrange rarefaction data for ggplot2                 
################################################################################

# read in rarefaction data output by mothur
mrare=read.table("data/07_mothur_rarefaction.txt", 
                              sep="\t", header=TRUE);

# load sample metadata
load("data/02_sample_metadata_formatted.Rdata");

# clean rarefaction data 
raref=mrare[, grepl("X0.03.*",names(mrare))]; 
colnames(raref)=gsub("X0.03.", "", colnames(raref)); 

raref=cbind(mrare$numsampled,raref);
colnames(raref)[1]="numsampled";

# melt data for ggplot2 (timepoint|sample|numOTUs)
rf<- reshape2::melt(raref, id.vars="numsampled",
                    value.name = "ASV_Richness");
colnames(rf)=c("TimeSampled", "sampleID", "ASV_Richness");
rare_meta=merge(rf, metadf, by="sampleID");

# add a Y=X line
b=rare_meta[rare_meta$sampleID=="P1",];
b$sampleID="P999"; b$ASV_Richness=b$TimeSampled; 
b$hyenaID2="Y=X";
rare_meta$hyenaID2=as.character(rare_meta$hyenaID2);
rare_meta=rbind(rare_meta,b);

# make hyenaID a factor for plotting purposes
rare_meta$hyenaID2=factor(rare_meta$hyenaID2, levels=c("M1","D1","G1",
                                                     "M2","D2","G2",
                                                     "M3","D3","G3",
                                                     "M4","D4","G4","Y=X"));


################################################################################
#             3. Make plot of sample rarefaction curves                 
################################################################################

# set color-palette (1 color for each hyena individual)
my_col=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
         "#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f","black");

# set tick marks for x-axis 
reads=seq(0, 3000, by=500);

# plot! 
rfc=ggplot(data=rare_meta) + 
  geom_line(mapping=aes(x=TimeSampled, 
                        y=ASV_Richness, 
                        col=hyenaID2, 
                        group=sampleID))+ 
  scale_x_continuous(breaks=reads)+
  lims(y=c(0,1500))+
  scale_colour_manual(values=my_col)+
  labs(x = "Number of Reads",
       y = "Number of ASVs",
       colour= "")+
  theme_bw() +
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme(panel.background= element_rect(colour = "black", size=2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right", 
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.x = element_text(size=14, face="bold")); plot(rfc);


################################################################################
#             4. Make histogram showing number of reads / sample               
################################################################################

# load ASV table to obtain total number of reads for each sample
load("data/01_ASV_table_filtered.Rdata");
tots=colSums(asvf);

# plot!
readseq = ggplot() + aes(tots)+
  geom_histogram(bins=10,             
                 color="black",
                 fill="#decbe4") +
  scale_x_continuous(breaks=seq(1116,30482,by=5000))+
  theme_classic() +
  labs(x = "Number of Reads per Sample",
       y = "Number of Samples")+
  theme(axis.title = element_text(size = 12, face="bold"), 
        axis.text = element_text(size = 11));plot(readseq);


################################################################################
#             4. save plots             
################################################################################

ggsave(filename="07_rarefaction_curves.pdf",
       device="pdf",path="./figures",
       plot=rfc,
       width=9,
       height=6,
       units="in",
       dpi=500);

ggsave(filename="07_NumReads_plot.pdf",
       device="pdf",path="./figures",
       plot=readseq,
       width=5,
       height=4,
       units="in",
       dpi=500);