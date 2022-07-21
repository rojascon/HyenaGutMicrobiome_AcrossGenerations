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

## CODE FOR: removing contaminant ASVs from dataset using sequences obtained from 
#  extraction kit controls

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load ASV Table, ASV Taxonomy, and sample metadata                 
################################################################################

# load ASV table output by DADA2 R tutorial 
#https://benjjneb.github.io/dada2/tutorial.html
load("data/00_seqtab_dada2.Rdata");

# load sample metadata 
meta=read.csv("data/00_sample_metadata.csv", stringsAsFactors = F);

# load ASV taxonomy file
#ASVs classified as Eukarya, Mitochondria, Chloroplast, and Unknown were removed
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T); 

# remove the 2 biological samples that did not amplify well [had <100 reads]
bad=c("P100","P168");
asvtbl=seqtab.nochim[!rownames(seqtab.nochim) %in% bad,];


################################################################################
#             2.  Identify potential contaminant bacterial ASVs using
#                                   R decontam 
################################################################################

# combine the 3 tables into a phyloseq object
rownames(tax)=tax$ASV; tax$ASV=NULL; tax=as.matrix(tax);
rownames(meta)=meta$sampleID; meta$sampleID=NULL;

ps=phyloseq(otu_table(asvtbl, 
                      taxa_are_rows =F),
            tax_table(tax), 
            sample_data(meta));

# identify potential contaminant ASVs using R decontam package
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# using the "prevalence" method & a 0.5 decontam score threshold
sample_data(ps)$is.neg=sample_data(ps)$type == "control";
scores <- isContaminant(ps, method="prevalence", neg="is.neg",
                        threshold = 0.5, detailed = TRUE);

# view # of contaminants identified by decontam
table(scores$contaminant);

# generate a histogram of ASV decontam scores to determine whether the
# 0.5 threshold was appropriate
scores$p.freq=NULL;
scores$prev_cat <- cut(scores$prev, c(0, 2, 5, 10, 9999), 
                       labels = c("2", "3-5", "6-10", "11+"));
prev_col <- c("2" = "#edf8e9", "3-5" = "#bae4b3", 
              "6-10" = "#74c476", "11+" = "#238b45");

prevhisto <-ggplot(scores,aes(x = p.prev, fill = prev_cat))+
  geom_histogram() + 
  labs(x='ASV decontam Score', 
       y='Number ASVs',
       fill="ASV Prevalence (# samples)") +
  theme_bw()+
  scale_fill_manual(values=prev_col) +
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.05)) +
  scale_y_continuous(breaks=c(0, 25,50, 100, 150, 200,235))+
  theme(legend.position = "bottom");

plot(prevhisto);
# the majority of ASVs have decontam scores > 0.5, meaning they are not contaminants
# 0.5 seems like an appropriate threshold , not too restrictive, not too flexible


################################################################################
#             3.  Decide on which contaminant ASVs will be 
#                         removed from dataset 
#################################################################################

# Criteria: 
# ASVs must have been identified as contaminants by decontam, 
# are present in at least 50% of control samples,
# have relative abundances >1% in control samples

# save ASVs with decontam scores <0.5 (these are the contaminant ASVs)
contam=scores[scores$contaminant=="TRUE",];

# filter ASV table to only include control samples
asvt=as.data.frame(t(asvtbl));
iscontrols=grep(pattern = "PB", x=colnames(asvt));
controls=asvt[, iscontrols]; 
controls=as.data.frame(controls);
save(controls, file="data/01_controls_ASV_table.Rdata");  ##will use in another script later

# convert that ASV table to presence absence
controlspa=(controls>0)*1; 
controlspa=as.data.frame(controlspa);

# retain ASVs present in at least 50% of control samples
controlspa$cprev=rowSums(controlspa);
controlspa=controlspa[controlspa$cprev>7,];

# retain ASVs that have >1% relative abundance across control samples
controls=apply(controls, 2, function(i) (i/sum(i))*100);

controls=controls[row.names(controls) 
                  %in% row.names(controlspa),];
controls=as.data.frame(controls);
controls$Abun=rowMeans(controls); 
controls=controls[controls$Abun>1,];

# retain ASVs with decontam scores <0.5
controls=controls[row.names(controls) 
                  %in% row.names(contam),];

# view the taxonomic classification of these ASVs
View(tax[row.names(tax) 
         %in% row.names(controls),]);


################################################################################
#             4.  Remove the selected contaminant ASVs from ASV table
#################################################################################

# remove contaminant ASVs
asvf=as.data.frame(t(asvtbl))
asvf=asvf[!row.names(asvf) %in% 
            row.names(controls),];

# also remove ASVs classified as chloroplast, mitochondria, and unknown
asvf=asvf[row.names(asvf) %in% row.names(tax),];

# remove control samples from ASV table now that they are no longer needed
asvf=asvf[,-iscontrols]; 
asvf=as.data.frame(asvf);

# remove ASVs with 0 reads
asvf$sm=rowSums(asvf);
asvf=asvf[asvf$sm>0,];
asvf$sm=NULL;

################################################################################
#             5.  save filtered ASV abundance table
#################################################################################

save(asvf, file="data/01_ASV_table_filtered.Rdata");

