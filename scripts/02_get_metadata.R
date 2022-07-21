#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 21 Jan 2021
#                         Last updated: 30 April 2021
#
################################################################################

## CODE FOR: formatting sample metadata (remove control samples, add prey
#  abundance data, and create factors)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load sample metadata and only retain samples of interest                
################################################################################

# load metadata
meta=read.csv("data/00_sample_metadata.csv", stringsAsFactors = F);

# remove control samples 
meta=meta[meta$type=="sample",];

# remove the 2 biological samples that did not amplify/sequence well (<400 reads)
bad=c("P100","P168");
meta=meta[!meta$sampleID %in% bad,];
meta$type=NULL;
metadf=meta;


################################################################################
#             2.  Format metadata factors                
################################################################################

metadf$birthdate=as.Date(metadf$birthdate, format="%d-%b-%y");
metadf$sample_date=as.Date(metadf$sample_date, format="%d-%b-%y");

metadf$hyenaID<-factor(metadf$hyenaID, levels=c("bsh","mrph","adon","mp","ffl",
                                                "tilt","coch","nav","roos","ua",
                                                "bail","hex"));
metadf$mat=factor(metadf$mat);

metadf$hyenaID2<-factor(metadf$hyenaID2, levels=c("M1","D1","G1",
                            "M2","D2","G2",
                            "M3","D3","G3",
                            "M4","D4","G4"));
metadf=metadf[order(metadf$sampleID),];


################################################################################
#             3.  Add mean monthly prey abundance data
################################################################################

# load prey density data from Holekamp lab data repository (private)
library(hyenadata); load_all_tables();
prey=tblPreyCensus; prey=prey[prey$region=="Narok",];
prey=prey[prey$transect %in% 
            c("WLOW","WHIGH","W3"),];
prey=prey[,c(3,6,9:25,27,38:39)];
prey[,3:22] <- sapply(prey[,3:22],as.numeric);

# sum across herbivore species to get a total # of herbivores for that day
prey$herb=rowSums(prey[ , c(3:20)], na.rm=TRUE);

# calculate mean monthly prey abundance
# e.g. mean number of herbivore prey for 30 days preceding sample date
metadf$preym=c();

for(i in 1:nrow(metadf))
{
  myv=metadf$sample_date[i]-30;
  A=mean(prey$herb[(prey$date>=myv) & 
                     (prey$date<=metadf$sample_date[i])], 
         na.rm=T);
  metadf$preym[i]=A;
};

# 4 samples did not have prey data --> NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
metadf$preym[is.nan(metadf$preym)] <- NA;


################################################################################
#             4. save formatted metadata file
################################################################################

save(metadf, file="data/02_sample_metadata_formatted.Rdata");

