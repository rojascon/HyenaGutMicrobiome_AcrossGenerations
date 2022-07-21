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

##CODE FOR: calculating 4 distance matrices for beta-diversity analyses:
#  Bray-Curtis (bray), Jaccard (jac), 
#  Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         formatted metadata table                 
################################################################################

load("data/01_ASV_table_filtered.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T); 
load("data/02_sample_metadata_formatted.Rdata");


################################################################################
#             2. Create phylogenetic tree of ASV sequences;
#     this is necessary for calculating Unifrac distances 
################################################################################

#the fitGTR steps are computationally expensive and may need high performance 
#computing.

# asvf=t(asvf);
# seqs <- getSequences(as.matrix(asvf)); 
# names(seqs) <- seqs;
# alignment <- AlignSeqs(DNAStringSet(seqs), 
#                         anchor=NA,verbose=FALSE);
# phangAlign <- phyDat(as(alignment, "matrix"), type="DNA"); 
# dm <- dist.ml(phangAlign); 
# treeNJ <- NJ(dm); 
# fit = pml(treeNJ, data=phangAlign);
# fitGTR <- update(fit, k=4, inv=0.2); 
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, 
#                     optGamma=TRUE, rearrangement = "stochastic", 
#                     control = pml.control(trace = 0));

#save the output as an .Rdata file
#save(fitGTR, file="data/05_ASV_phylotree.Rdata");


################################################################################
#             3. Remove singleton and doubleton ASVs
################################################################################

#these are ASVs that only appear once or twice in entire dataset
tempdf=(asvf>0)*1; 
tempdf=as.data.frame(tempdf);
tempdf=tempdf[rowSums(tempdf)>2,];

asvf=asvf[rownames(asvf) %in% rownames(tempdf),];

##remove samples that do not have prey availability data
rem=c("P173","P193","P253","P302");
asvf=as.data.frame(asvf[,
                        -which(names(asvf) %in% rem)]);


################################################################################
#             4. Generate distance matrices (bray, jaccard, unifrac)
################################################################################

#transpose ASV table so ASVs are columns
asvf=t(asvf);

# BRAY-CURTIS distance
bray<-apply(asvf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

# JACCARD distance
jac=(asvf>0)*1;
print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

# create phyloseq object for Unifrac distances
load("data/05_ASV_phylotree.Rdata");
wuni<-phyloseq(otu_table(asvf, taxa_are_rows=FALSE), 
               phy_tree(fitGTR$tree));

set.seed(711);
phy_tree(wuni) <- root(phy_tree(wuni), 
                       sample(taxa_names(wuni), 1), resolve.root = TRUE)
is.rooted(phy_tree(wuni));

# WEIGHTED UNIFRAC distance
wuni.dist=UniFrac(wuni, 
                  weighted=TRUE, 
                  normalized=TRUE);

# UNWEIGHTED UNIFRAC distance
unwuni.dist=UniFrac(wuni, 
                    weighted=FALSE, 
                    normalized=TRUE);


################################################################################
#             5. save the matrices as an R object
################################################################################

save("bray.dist", "jac.dist",
     "unwuni.dist","wuni.dist", 
     file = "data/05_dissimilarity_distances.Rdata");


