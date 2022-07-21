#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 21 Jan 2021
#                         Last updated:  18 Feb 2021
#
################################################################################

## CODE FOR: filtering sample metadata of the larger dataset to only include the 
#             32 samples that have metagenomic data

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load sample metadata and only retain samples of interest                
################################################################################

# load metadata
load("data/02_sample_metadata_formatted.Rdata");

# retain samples of interest
mysamples=c("P102","P107","P14","P143","P149","P167","P178","P179","P186",
                 "P190","P191","P196","P212","P225","P261","P265","P266","P276",
                 "P291","P294","P295","P301","P304","P41","P43","P46","P47",
                 "P74","P75","P82","P84","P9");
meta=metadf[metadf$sampleID %in% mysamples,];

meta$hyenaID2<-factor(meta$hyenaID2, levels=c("D1","G1","D3","G3"));
meta=meta[order(meta$sampleID),];

# save metadata file
save(meta, file="data/metagenomes/00_metagenomes_metadata.Rdata");
