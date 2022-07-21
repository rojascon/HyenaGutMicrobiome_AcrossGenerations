#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 10 Feb 2021
#                         Last updated:  27 May 2021
#
################################################################################

## CODE FOR: a) conducting PERMANOVA tests to determine which host factors predict
#                 gut microbiome similarity
#  b) plotting the model output of PERMANOVA (e.g. R2 values) as a barplot

# DISTANCES USED:
# Bray-Curtis (bray), Jaccard (jac), and Weighted Unifrac (wuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load dissimilarity matrices generated in script 05a
#                         & formatted metadata table                 
################################################################################

load("data/02_sample_metadata_formatted.Rdata");
load("data/05_dissimilarity_distances.Rdata");

# remove samples without prey availability from metadata
meta2=metadf[complete.cases
             (metadf$preym),];

# rename matrilines (e.g. 1 to "M1"
meta2$mat=as.character(meta2$mat);
meta2$mat2=paste("M",meta2$mat,sep="");
meta2$mat2<-factor(meta2$mat2, levels=c("M1","M2","M3","M4"));


################################################################################
#                     2. run marginal PERMANOVAs
#                       microbiota similarity ~ 
#     hyena ID + age + maternal relatedness + monthly prey + calendar year
################################################################################

# set variables for loop
mydist=list(jac.dist,bray.dist, wuni.dist);
names=c("Jaccard","Bray-Curtis","W Unifrac");
met=c("jaccard","bray","bray") ;

# run for loop to conduct the 3 PERMANOVA tests
for(i in 2:3)
{
  print(paste("PERMANOVA test, N=301, using:", names[i]));
  print(adonis(mydist[[i]]~ 
                 mat2+
                 hyenaID2+
                 age_yrs+
                 preym+
                 year,
               data=meta2,
               method = met[i],
               by="margin",
               permutations = 999));
};


#################################################################################
#             3. Barplots of % variance explained by host predictors
##                      (e.g. PERMANOVA R2 values)
################################################################################

# make dataframe of predictors and their R2 values
df=data.frame("Predictors"=c("Hyena ID","Age", "Matriline",
                             "Year","PreyAbun"), 
               "variance"=c(11.31,2.55,4.34,0.00001,0.000001));

# plot
perm.ecol=ggplot(df,aes(x=reorder(Predictors, -variance), y=variance))+
  geom_bar(stat="identity", 
           color="black",
           fill="#9ecae1")+ ##fcbba1
  labs(y="% variance explained",
       x="Host Predictors")+
  theme_classic() + 
  theme(axis.title = element_text(size = 12, face="bold"), 
        axis.text = element_text(size = 11),
        axis.text.x=element_text(angle = 55, vjust = 0.60),
        legend.position="none");
plot(perm.ecol);


# save plot
ggsave(filename="05_PERMANOVA_plot.pdf",
       device="pdf",path="./figures",
       plot=perm.ecol,
       width=3,
       height=3.6,
       units="in",
       dpi=500);
