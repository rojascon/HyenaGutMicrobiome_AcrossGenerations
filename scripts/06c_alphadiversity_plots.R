#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         By: Connie Rojas
#                         Created: 17 Feb 2022
#                         Last updated: 30 May 2022
#
################################################################################

## CODE FOR: 
  # A) plots of microbiome alpha-diversity ~ hyena identity
  # B) plots of microbiome alpha-diversity ~ age (yrs)
  # C) plots of microbiome alpha-diversity ~ prey density
# all based on statistical output from script 06b

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load sample metadata, and table of alpha diversity metrics
################################################################################

load("data/02_sample_metadata_formatted.Rdata");
load("data/06_alpha_values.Rdata");

# add metadata to alpha-diversity data
# remove 3 samples without prey data
am=inner_join(metadf[complete.cases(metadf$preym),], 
              alpha[,c(1,2,4,8)], by="sampleID");


################################################################################
#             2. mean microbiome diversity ~ host identity
################################################################################

# calculate mean chao richness for every hyena individual
df=am[,c(5,12)]; 
b=aggregate(Chao1~hyenaID2, df, mean); colnames(b)[2]="mean_alpha";
c=aggregate(Chao1~hyenaID2, df, sd); colnames(c)[2]="sd_alpha";
d=merge(b,c,by="hyenaID2");

# calculate standard error of the mean
d$size=table(am$hyenaID2);
d$err=d$sd_alpha/sqrt(d$size);
d$upper=d$mean_alpha+ d$err;
d$lower=d$mean_alpha-d$err;

# plot mean alpha diversity in ggplot
alpha1=
  ggplot(d,aes(x=hyenaID2, y=mean_alpha,ymin=100))+ 
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  theme_classic()+
  labs(y="Chao 1 Richness",
       x="Hyena identity",
       title="A)")+
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=15, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"));plot(alpha1);


################################################################################
#             2. mean microbiome diversity ~ hyena age (yrs)
################################################################################

# round host age (11.5 yrs > 11 yrs)
am$agro=round(am$age_yrs);
am$agro[am$agro>14]=15;

# calculate average chao richness for every age
df=am[,c(12,15)]; 
b=aggregate(Chao1~agro, df, mean); colnames(b)[2]="mean_alpha";
c=aggregate(Chao1~agro, df, sd); colnames(c)[2]="sd_alpha";
d=merge(b,c,by="agro");

# calculate standard error of the mean
d$size=table(am$agro);
d$err=d$sd_alpha/sqrt(d$size);
d$upper=d$mean_alpha+ d$err;
d$lower=d$mean_alpha-d$err;

# create plot of the means using ggplot
alpha2=
  ggplot(d,aes(x=agro, y=mean_alpha,ymin=100))+ 
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  geom_smooth(method=lm, color="#737373", fill="#bdbdbd",se=TRUE)+
  theme_classic()+
  labs(y="Phylogenetic Diversity",
       x="Age (yrs)",
       title="B)")+
  scale_x_continuous(breaks=seq(2,15, by=4))+
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=15, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"));plot(alpha2);


################################################################################
#             3. mean microbiome diversity ~ prey density
################################################################################

# categorize prey density values into 7 bins
breaks <- seq(0,1094, by = 100); breaks=c(breaks, 1094);
breaks=c(0,50,100,200,300,400,500,1094);

am$preyb=cut(round(am$preym),
           breaks=breaks,
           include.lowest=TRUE,
           labels = c("<50", "50-100", "100-200", 
                      "200-300","300-400","400-500",
                      ">500"));

# calculate average chao richness for every prey density bin
df=am[,c(12,16)]; 
b=aggregate(Chao1~preyb, df, mean); colnames(b)[2]="mean_alpha";
c=aggregate(Chao1~preyb, df, sd); colnames(c)[2]="sd_alpha";
d=merge(b,c,by="preyb");

# calculate standard error of the mean
d$size=table(am$preyb);
d$err=d$sd_alpha/sqrt(d$size);
d$upper=d$mean_alpha+ d$err;
d$lower=d$mean_alpha-d$err;

# assign bins a numerical value so that you can add a trendline to plot
d$preyb2[d$preyb=="<50"]=1;d$preyb2[d$preyb==">500"]=7; 
d$preyb2[d$preyb=="100-200"]=3; d$preyb2[d$preyb=="200-300"]=4;
d$preyb2[d$preyb=="300-400"]=5;d$preyb2[d$preyb=="400-500"]=6;
d$preyb2[d$preyb=="50-100"]=2;

# create plot of the means using ggplot
alpha3=
  ggplot(d,aes(x=preyb2, y=mean_alpha, ymin=100))+
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  geom_smooth(method=lm, color="#737373", fill="#bdbdbd", se=TRUE)+
  scale_x_continuous(breaks=d$preyb2,labels=d$preyb)+
  theme_classic()+
  labs(y="Chao 1 Richness",
       x="Mean Monthly Prey Abundance",
       title="C)")+
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=15, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11, angle = 45, vjust = 0.66),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"));plot(alpha3);


################################################################################
#             4. save plots into a single figure
################################################################################
aps=arrangeGrob(alpha1,alpha2,alpha3, nrow=1); 

ggsave(filename="06_alphadiv_plots.pdf",
       device="pdf",path="./figures",
       plot=aps,
       width=9.8,
       height=3.2,
       units="in",
       dpi=500);
