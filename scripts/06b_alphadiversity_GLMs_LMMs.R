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
#  A) running GLMs to see which host factors predict gut microbiota 
#     alpha-diversity

#  B) running LMMs to see which host factors predict gut microbiota 
#     alpha-diversity after controling for host identity and sample year

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

# calculate means and standard deviations
aggregate(am$PD, list(am$hyenaID2), FUN=mean);
aggregate(am$PD, list(am$hyenaID2), FUN=sd);


################################################################################
#             2. Run linear mixed models
#    alpha-diversity ~ hyena ID + matriline + age + prey density + year
################################################################################

full.model.chao <- glm(log(Chao1) ~hyenaID2+ mat+ age_yrs+ preym + year, 
                       data = am, na.action = "na.omit"); 

full.model.shannon <- glm(log(Shannon) ~hyenaID2+ mat+ age_yrs+ preym + year, 
                          data = am, na.action = "na.omit");

full.model.pd <- glm(log(PD) ~hyenaID2+ mat+ age_yrs+ preym + year, data = am,
                      na.action = "na.omit");

##assess significance of host predictors
Anova(full.model.chao);
Anova(full.model.shannon);
Anova(full.model.pd);

#make sure model residuals look normal
qqnorm(resid(full.model.chao));
qqnorm(resid(full.model.shannon));
qqnorm(resid(full.model.pd));


################################################################################
#             3. Run linear mixed models
#    alpha-diversity ~  matriline + age + prey density
#       hyena identity and year as a random factor
################################################################################
reduced.model.chao <- lmer(log(Chao1) ~ mat+ age_yrs + preym +(1|year)+(1|hyenaID), 
                           data = am, na.action = "na.fail");

reduced.model.shannon <- lmer(log(Shannon) ~ mat+ age_yrs + preym +(1|year)+
                                (1|hyenaID), data = am, na.action = "na.fail"); 

reduced.model.pd <- lmer(log(PD) ~ mat+ age_yrs + preym +(1|year)+(1|hyenaID), 
                         data = am, na.action = "na.fail");

##assess significance of host predictors
Anova(reduced.model.chao);
Anova(reduced.model.shannon);
Anova(reduced.model.pd);

#make sure model residuals look normal
qqnorm(resid(reduced.model.chao));
qqnorm(resid(reduced.model.shannon));
qqnorm(resid(reduced.model.pd));
