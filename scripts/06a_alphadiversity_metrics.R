#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 10 Feb 2021
#                         Last updated: 17 Feb 2021
#
################################################################################

## CODE FOR: 
#  A) subsampling samples to 2,900 sequences each using mothur software
#  B) calculating Chao1 richness and Shannon diversity using phyloseq
#  C) calculating Faith's PD using picante

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Subsample samples to 2900 sequences using mothur   
#         this number was chosen because it was the 2nd lowest number of 
#                   sequences found in a biological sample 
################################################################################

# load ASV table and ASV taxonomy
load("data/01_ASV_table_filtered.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T);

# make a df of ASV sequences and ASV labels
asvdf=as.data.frame(t(asvf));
names=data.frame("seqs" = colnames(asvdf), 
                 "names" = sprintf("ASV%s",
                                   seq(1:ncol(asvdf))));

# format ASV table for mothur
# step1: replace colnames of ASV table with ASV labels from above
# step2:add the label, Group, and numOTUs columns to ASV table
colnames(asvdf)=names$names[match(names(asvdf),names$seqs)];
tempdf=data.frame("label" = rep(0.03,nrow(asvdf)), 
                  "Group" = row.names(asvdf), 
                  "numOtus" = rep(ncol(asvdf),nrow(asvdf)));
mothur=cbind(tempdf, asvdf);
write.table(mothur, file="data/06_input_mothur.txt", row.names=F, sep="\t");

# run mothur
# step1: for some reason, must open the above .txt file in Excel, and 
#           save as Tab delimited file (.txt)
# step2: download and open mothur software (Schloss et.al 2009)
# step3: run the command: 
#           sub.sample(shared=06_input_mothur.txt, size=2900, persample=true)
# step4: rename output file as "05_output_mothur.txt") and place in the data directory


################################################################################
#             2. Clean mothur file output
################################################################################

# read in output from mothur
mothur2=read.table("data/06_output_mothur.txt", sep="\t",header=T);
rownames(mothur2)=mothur2$Group;
mothur2=mothur2[,-c(1:3)];

# replace ASV labels with ASV sequence names
colnames(mothur2)=names$seqs[match
                             (names(mothur2),names$names)];
mothur2=as.matrix(mothur2);


################################################################################
#             3. alpha-diversity metrics with phyloseq
################################################################################

# make a phyloseq object and calculate 4 metrics of alpha-diversity
ps<- phyloseq(otu_table(mothur2, taxa_are_rows=FALSE));
alpha=estimate_richness(ps,split = TRUE, 
                        measures = c("Chao1","Shannon"));

# calculate Good's coverage and append to alpha df
gcov=goods(mothur2);
alpha=merge(alpha, gcov, by="row.names", all=TRUE); 
colnames(alpha)[1]="sampleID";


################################################################################
#             4. Calculate Faith's Phylogenetic diversity using picante
################################################################################

#load phylogenetic tree of ASVs sequences
#see script "05a_betadiversity_distances.R" for how to generate
load("data/05_ASV_phylotree.Rdata"); 

#calculate Faith's phylogenetic diversity and append to alpha df
faith=pd(mothur2, phy_tree(fitGTR$tree), include.root=F);
faith$sampleID=row.names(faith);
alpha=inner_join(alpha, faith, by="sampleID");


################################################################################
#             5. save df of alpha diversity values 
################################################################################

save(alpha, file="data/06_alpha_values.Rdata");


