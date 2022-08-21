#RCode for Castillioni et al. (2022). Effects of drought on grassland phenology depend on functional types. New Phytologist.

library(dplyr)
library(ggplot2)
library(lmerTest)
library(lme4)
library(car)
require(devtools)
library(ggpubr)
library(MASS)
library(performance)
library(glmmTMB)
library(tidyverse)
library(lubridate)
library(bbmle)
library(emmeans)
library(RAM)
library(report)
library(sjPlot)
library(sjmisc)
theme_set(theme_sjplot())

#dataset#
phenology_data <- read.csv("Castillioni_etal_2022_phenology_data.csv", header= TRUE)
loral_abundance_data <- read.csv("Castillioni_etal_2022_floral_abundance.csv", header = TRUE)
fruit_abundance_data <- read.csv("Castillioni_etal_2022_fruit_abundance.csv", header = TRUE)
dat_seed_data <- read.csv("Castillioni_etal_2022_seed_viability.csv", header = TRUE)

#-------------- G R O W I N G   D E G R E E   D A Y   --------------#
     gdd_model <- lmer(log(gdd)~precipitation + (1|block:plot), data=phenology_data)
     qqnorm(resid(gdd_model))
     qqline(resid(gdd_model))
     anova(gdd_model)
     summary(gdd_model)
     r2(gdd_model)
     
#models:block or plot is removed if variance and stdv is zero (check model summary)     
#-------------- C O M M U N I T Y  -  L E V E L --------------#
     #--PEAK FLOWER--#     
     peakflower_community_glmmTMB <- glmmTMB(peak_flower ~ precipitation + (1|block/plot) + (1|species), family= poisson(link=log), data = phenology_data)
     check_overdispersion(peakflower_community_glmmTMB)
     Anova(peakflower_community_glmmTMB, type="III")
     
     #--DURATION FLOWER--#
     duration_flower_community_glmmTMB <- glmmTMB( duration_flower ~  precipitation  + (1|block/plot) + (1|species), family=nbinom2, data = phenology_data)
     check_overdispersion(duration_flower_community_glmmTMB)
     Anova(duration_flower_community_glmmTMB, type="III")
     
     #--DURATION FRUIT--#
     duration_fruit_community_glmmTMB <- glmmTMB(duration_fruit ~  precipitation + (1|block/plot) + (1|species), family=nbinom2, data = phenology_data)
     check_overdispersion(duration_fruit_community_glmmTMB)
     Anova(duration_fruit_community_glmmTMB, type="III")
     
     #--FLOWER ABUNDANCE--#
     floral_abundance_community_glmmTMB <- glmmTMB( flower ~  precipitation +  (1|species) + (1|block/plot), family=nbinom2,  data= floral_abundance_data) 
     check_overdispersion(floral_abundance_community_glmmTMB)
     Anova(floral_abundance_community_glmmTMB, type="III")
     
     #--FRUIT ABUNDANCE--#
     fruit_abundance_community_glmmTMB <- glmmTMB( fruit ~  precipitation +  (1|species) + (1|block/plot), family=nbinom2, data= fruit_abundance_data) 
     check_overdispersion(fruit_abundance_community_glmmTMB)
     Anova(fruit_abundance_community_glmmTMB, type="III")
     
     #--SEED VIABILITY--#
     proportion_viable_community_glmmTMB <- glmmTMB(proportion_viable ~ precipitation  + (1|block/plot) + (1|species), family=betabinomial,  weights = total_seeds_sum, data = dat_seed_data)
     check_overdispersion(proportion_viable_community_glmmTMB)
     Anova(proportion_viable_community_glmmTMB, type="III")
     
     #-------------- B L O O M   T I M E --------------#
    #--PEAK FLOWER--#     `
     peakflower_bloomtime_glmmTMB <- glmmTMB(peak_flower ~  scale(precipitation) * scale(first_flower_mean) +  (1|species) + (1|block/plot), family= poisson(link=log),  data = phenology_data)
     check_overdispersion(peakflower_bloomtime_glmmTMB)
     Anova(peakflower_bloomtime_glmmTMB, type="III")
     summary(peakflower_bloomtime_glmmTMB)
    
    #Shift in phenology: 
       #early
     summary(glmmTMB(peak_flower~precipitation+(1|species)+(1|block/plot),family= poisson(link=log),data = filter(phenology_data, early_late == "early")))
     exp(5.2094227 + (0.0002682*-60)) - exp(5.2094227 + (0.0002682*0)) #[1] -2.921083
        #late
     summary(glmmTMB(peak_flower~precipitation+(1|species)+(1|block/plot),family= poisson(link=log),data = filter(phenology_data, early_late == "late")))
     exp(5.628e+00 + (-1.778e-04*-60)) - exp(5.628e+00 + (-1.778e-04*0)) #[1] 2.982709
  
     #--DURATION FLOWER--#
     duration_flower_bloomtime_glmmTMB <- glmmTMB( duration_flower ~  scale(precipitation) * scale(first_flower_mean) +  (1|species) + (1|block/plot),  family=nbinom2, data = phenology_data)
     check_overdispersion(duration_flower_bloomtime_glmmTMB)
     Anova(duration_flower_bloomtime_glmmTMB, type="III")
     
     #--DURATION FRUIT--#
     duration_fruit_bloomtime_glmmTMB <- glmmTMB(duration_fruit ~  scale(precipitation) * scale(first_flower_mean)  +  (1|species) + (1|block/plot),  family=nbinom2, data = phenology_data)
     check_overdispersion(duration_fruit_bloomtime_glmmTMB)
     Anova(duration_fruit_bloomtime_glmmTMB, type="III")
     
     #--FLOWER ABUNDANCE--#
     floral_abundance_bloomtime_glmmTMB <- glmmTMB( flower ~  scale(precipitation) * scale(first_flower_mean) +  (1|species) + (1|block/plot), family=nbinom2,   data= floral_abundance_data) #keep the interaction model
     check_overdispersion(floral_abundance_bloomtime_glmmTMB)
     Anova(floral_abundance_bloomtime_glmmTMB, type="III")
     
     #--FRUIT ABUNDANCE--#
     fruit_abundance_bloomtime_glmmTMB <- glmmTMB( fruit ~  scale(precipitation) * scale(first_flower_mean) +  (1|species) + (1|block/plot), family=nbinom2,  data= fruit_abundance_data) #keep the interaction model
     check_overdispersion(fruit_abundance_bloomtime_glmmTMB)
     Anova(fruit_abundance_bloomtime_glmmTMB, type="III")
     
     #--SEED VIABILITY--#
     proportion_viable_bloomtime_glmmTMB <- glmmTMB(proportion_viable ~  scale(precipitation) * scale(first_flower_mean) +  (1|species) + (1|block/plot),  family=betabinomial, weights = total_seeds_sum, data = dat_seed_data)
     check_overdispersion(proportion_viable_bloomtime_glmmTMB)
     Anova(proportion_viable_bloomtime_glmmTMB, type="III")
     
     #-------------- P H O T O S Y N T H E T I C   P A T H W A Y ---------------#
    #--PEAK FLOWER--#     
         peakflower_func_group_glmmTMB <- glmmTMB(peak_flower~precipitation*func_group+(1|species)+(1|block/plot),family= poisson(link=log),data = phenology_data)
         check_overdispersion(peakflower_func_group_glmmTMB)
         Anova(peakflower_func_group_glmmTMB, type="III")
         summary(peakflower_func_group_glmmTMB)
         
     #Shit in phenology
         #C3
         peakflower_c3_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|species)+(1|block/plot),family= poisson(link=log),data = filter(phenology_data, func_group=="C3"))
         summary(peakflower_c3_glmmTMB)
         exp(5.3243842 + (0.0001296*-60)) - exp(5.3243842 + (0.0001296*0)) #[1] -1.590082
         #C4
         peakflower_c4_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|species)+(1|block/plot),family= poisson(link=log),data = filter(phenology_data, func_group=="C4"))
         summary(peakflower_c4_glmmTMB)
         exp(5.6282545 + (-0.0001971*-60)) - exp(5.6282545 + (-0.0001971*0)) #[1] 3.30924
         
     #--DURATION FLOWER--#
         duration_flower_func_group_glmmTMB <- glmmTMB( duration_flower~precipitation*func_group+ (1|species)+(1|block/plot),family=nbinom2,data = phenology_data)
         check_overdispersion(duration_flower_func_group_glmmTMB)
         Anova(duration_flower_func_group_glmmTMB, type="III")
         
     #--DURATION FRUIT--#
         duration_fruit_func_group_glmmTMB <- glmmTMB(duration_fruit~precipitation*func_group+(1|species)+(1|plot),family= nbinom2,data = phenology_data)
         check_overdispersion(duration_fruit_func_group_glmmTMB)
         Anova(duration_fruit_func_group_glmmTMB, type="III")
         summary(duration_fruit_func_group_glmmTMB)
         
      #shift in phenology
         #C3
         duration_fruit_C3_glmmTMB <- glmmTMB(duration_fruit~precipitation+(1|species)+(1|plot),family= nbinom2,data = filter(phenology_data, func_group=="C3"))
         summary(duration_fruit_C3_glmmTMB)
         exp(3.7971821 + (-0.0002545*-60)) - exp(3.7971821 + (-0.0002545*0)) #[1] 0.6858898
         #C4
         duration_fruit_C4_glmmTMB <- glmmTMB(duration_fruit~precipitation+(1|species)+(1|plot),family= nbinom2,data = filter(phenology_data, func_group=="C4"))
         summary(duration_fruit_C4_glmmTMB)
         exp(3.287563 + (0.003167*-60)) - exp(3.287563 + (0.003167*0)) #[1] -4.634049
         
    #--FLOWER ABUNDANCE--#
         floral_abundance_func_group_glmmTMB <- glmmTMB(flower~precipitation*func_group+(1|species)+(1|block/plot),family=nbinom2,data= floral_abundance_data)
         check_overdispersion(floral_abundance_func_group_glmmTMB)
         Anova(floral_abundance_func_group_glmmTMB, type="III")
         summary(floral_abundance_func_group_glmmTMB)
      
      #Shift in phenology
         #C3
         floral_abundance_C3_glmmTMB <- glmmTMB(flower~precipitation+(1|species)+(1|block/plot),family=nbinom2,data= filter(floral_abundance_data, func_group=="C3"))
         summary(floral_abundance_C3_glmmTMB)
         exp(1.952631 + (-0.001434*-60)) - exp(1.952631 + (-0.001434*0)) #[1] 0.6331908
         #C4
         floral_abundance_C4_glmmTMB <- glmmTMB(flower~precipitation+(1|species)+(1|block/plot),family=nbinom2,data= filter(floral_abundance_data, func_group=="C4"))
         summary(floral_abundance_C4_glmmTMB)
         exp(2.219434 + (0.005696*-60)) - exp(2.219434 + (0.005696*0)) #[1] -2.663842
         
    #--FRUIT ABUNDANCE--#
         fruit_abundance_func_group_glmmTMB <- glmmTMB(fruit~precipitation*func_group+(1|block/plot),family=nbinom2,data= fruit_abundance_data) 
         check_overdispersion(fruit_abundance_func_group_glmmTMB)
         Anova(fruit_abundance_func_group_glmmTMB, type="III")
         
    #--SEED VIABILITY--#
         proportion_viable_func_group_glmmTMB <- glmmTMB(proportion_viable~ precipitation*func_group+ (1|species)+(1|block/plot),weights=total_seeds_sum,family= betabinomial, data = dat_seed_data)
         check_overdispersion(proportion_viable_func_group_glmmTMB)
         Anova(proportion_viable_func_group_glmmTMB, type="III")
     
#-------------- F O R B  -  G R A S S --------------#
      #--PEAK FLOWER--#     
          peakflower_forb_grass_glmmTMB <- glmmTMB(peak_flower~precipitation*forb_grass+(1|species)+(1|block/plot),family= poisson(link=log),data = phenology_data)
           check_overdispersion(peakflower_forb_grass_glmmTMB)
           Anova(peakflower_forb_grass_glmmTMB, type="III")
           
      #--DURATION FLOWER--#
      duration_flower_forb_grass_glmmTMB <- glmmTMB( (duration_flower) ~  precipitation * forb_grass +  (1|species) + (1|block/plot),  family=nbinom2, data = phenology_data)
           check_overdispersion(duration_flower_forb_grass_glmmTMB)
           Anova(duration_flower_forb_grass_glmmTMB, type="III")
           
      #--DURATION FRUIT--#
      duration_fruit_forb_grass_glmmTMB <- glmmTMB((duration_fruit) ~  precipitation * forb_grass +  (1|species) + (1|block/plot),  family=nbinom2, data = phenology_data)
           check_overdispersion(duration_fruit_forb_grass_glmmTMB)
           Anova(duration_fruit_forb_grass_glmmTMB, type="III")
     
#--FLOWER ABUNDANCE--#
    floral_abundance_forb_grass_glmmTMB <- glmmTMB(flower~precipitation*forb_grass+(1|species)+(1|block/plot),family=nbinom2,data= floral_abundance_data) 
     check_overdispersion(floral_abundance_forb_grass_glmmTMB)
     Anova(floral_abundance_forb_grass_glmmTMB, type="III")
     summary(floral_abundance_forb_grass_glmmTMB)
     
     #forb
     floral_abundance_forb_glmmTMB <- glmmTMB(flower~precipitation+(1|species)+(1|block/plot),family=nbinom2,data= filter(floral_abundance_data, forb_grass == "forb")) 
     summary(floral_abundance_forb_glmmTMB)
     exp(2.101406 + (-0.002381*-60)) - exp(2.101406 + (-0.002381*0)) #[1] 1.255829
     #grass
     floral_abundance_grass_glmmTMB <- glmmTMB(flower~precipitation+(1|species)+(1|plot),family=nbinom2,data= filter(floral_abundance_data, forb_grass == "grass")) 
     summary(floral_abundance_grass_glmmTMB)
     exp(1.961707 + (0.005510*-60)) - exp(1.961707 + (0.005510*0)) #[1] -2.001928
     
  #--FRUIT ABUNDANCE--#
     fruit_abundance_forb_grass_glmmTMB <- glmmTMB(fruit~precipitation*forb_grass+(1|species)+(1|block/plot),family=nbinom2,data=fruit_abundance_data) 
     check_overdispersion(fruit_abundance_forb_grass_glmmTMB)
     Anova(fruit_abundance_forb_grass_glmmTMB, type="III")
     summary(fruit_abundance_forb_grass_glmmTMB)
     
     #forb
     fruit_abundance_forb_glmmTMB <- glmmTMB(fruit~precipitation+(1|species)+(1|block/plot),family=nbinom2,data=filter(fruit_abundance_data, forb_grass =="forb")) 
     summary(fruit_abundance_forb_glmmTMB)
     exp(2.703120 + (-0.001348*-60)) - exp(2.703120 + (-0.001348*0)) #[1] 1.257397
     #grass
     fruit_abundance_grass_glmmTMB <- glmmTMB(fruit~precipitation+(1|species)+(1|block/plot),family=nbinom2,data=filter(fruit_abundance_data, forb_grass =="grass")) 
     summary(fruit_abundance_grass_glmmTMB)
     exp(1.776232 + (0.005598*-60)) - exp(1.776232 + (0.005598*0)) #[1] -1.685373
     
    #--SEED VIABILITY--#
    proportion_viable_forb_grass_glmmTMB <- glmmTMB(proportion_viable~precipitation*forb_grass+(1|species)+(1|block/plot),family=betabinomial, weights = total_seeds_sum, data = dat_seed_data)
     check_overdispersion(proportion_viable_forb_grass_glmmTMB)
     Anova(proportion_viable_forb_grass_glmmTMB, type="III")
     
#####- S P E C I E S   S P E C I F I C  -  L E V E L --------------####
     #species: Levels: AP BI CF CM CS DO ES SC SE SN SS

      #--PEAK FLOWER--#     
    peakflower_glmmTMB <- glmmTMB(peak_flower ~  precipitation*species + (1|block/plot) ,  family= poisson(link=log), data= phenology_data) 
    check_overdispersion(peakflower_glmmTMB)
    summary(peakflower_glmmTMB)
    Anova(peakflower_glmmTMB, type="III")

          #phenological shift
          peakflower_AP_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|plot),family=poisson(link=log),data = subset(phenology_data, species=="AP")) #subset focal species
          summary(peakflower_AP_glmmTMB)
          #phenology shift calculation:
          exp((5.5521669)+(-0.0004530*-60))-exp((5.5521669)+(-0.0004530*0)) #AP
          #[1] 7.102976
          
          peakflower_BI_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|block/plot),family=poisson(link=log),data = subset(phenology_data, species=="BI"))
          summary(peakflower_BI_glmmTMB)
          #phenology shift calculation:
          exp((5.6483961)+(-0.0002282*-60))-  exp((5.6483961)+(-0.0002282*0)) #BI
          #[1] 3.913008
          
          peakflower_CF_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|block/plot),family=poisson(link=log),data = subset(phenology_data, species=="CF"))
          summary(peakflower_CF_glmmTMB)
          #phenology shift calculation:
          exp((5.3703534)+(-0.0003753*-60)) - exp((5.3703534)+(-0.0003753*0)) #CF
          #[1] 4.894897
          
          peakflower_CM_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|block/plot),family=poisson(link=log),data = subset(phenology_data, species=="CM"))
          summary(peakflower_CM_glmmTMB)
          #phenology shift calculation:
          exp((5.5068077)+(0.0009393*-60)) - exp((5.5068077)+(0.0009393*0))#CM
          #[1] -13.50054
          
          peakflower_CS_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|block/plot),family=poisson(link=log),data = subset(phenology_data, species=="CS"))
          summary(peakflower_CS_glmmTMB)
          #phenology shift calculation:
          exp((5.126e+00)+(-6.822e-05*-60)) - exp((5.126e+00)+(-6.822e-05*0))#CS
          #[1] 0.6904713
          
          peakflower_DO_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|plot),family=poisson(link=log),data = subset(phenology_data, species=="DO"))
          summary(peakflower_DO_glmmTMB)
          #phenology shift calculation:
          exp((5.0042558)+(0.0002202*-60)) - exp((5.0042558)+(0.0002202*0)) #DO
          #[1] -1.956246
          
          peakflower_ES_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|plot),family=poisson(link=log),data = subset(phenology_data, species=="ES"))
          summary(peakflower_ES_glmmTMB)
          #phenology shift calculation:
          exp((5.0329275)+(-0.0001301*-60)) - exp((5.0329275)+(-0.0001301*0)) #ES
          #[1] 1.20198
          
          peakflower_SE_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|plot),family=poisson(link=log),data = subset(phenology_data, species=="SE"))
          summary(peakflower_SE_glmmTMB)
          #phenology shift calculation:
          exp((5.672e+00)+(-5.302e-05*-60)) - exp((5.672e+00)+(-5.302e-05*0)) #SE
          #[1] 0.9259771
          
          peakflower_SN_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|plot),family=poisson(link=log),data = subset(phenology_data, species=="SN"))
          summary(peakflower_SN_glmmTMB)
          #phenology shift calculation:
          exp((5.6120431)+(-0.0001293*-60)) - exp((5.6120431)+(-0.0001293*0))#SN
          #[1] 2.131645
          
          peakflower_SS_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|plot),family=poisson(link=log),data = subset(phenology_data, species=="SS"))
          summary(peakflower_SS_glmmTMB)
          #phenology shift calculation:
          exp((5.6223570)+(-0.0002854*-60)) - exp((5.6223570)+(-0.0002854*0))#SS
          #[1] 4.776256

    #--PEAK FRUIT--#
    peak_fruit_glmmTMB <- glmmTMB(peak_fruit ~  precipitation * species  + (1|plot) , family= poisson(link=log), data= phenology_data) 
    check_overdispersion(peak_fruit_glmmTMB)
    Anova(peak_fruit_glmmTMB, type="III")
    summary(peak_fruit_glmmTMB)
    
    #--DURATION FLOWER--#
    duration_flower_glmmTMB <- glmmTMB( duration_flower ~  precipitation * species , family= nbinom2, data= phenology_data) 
    check_overdispersion(duration_flower_glmmTMB)
    Anova(duration_flower_glmmTMB, type="III")
    
    #--FLOWER ABUNDANCE--#
    floral_abundance_glmmTMB <- glmmTMB( (flower) ~  precipitation * species +  (1|block/plot) , family=nbinom2, data= floral_abundance_data) 
    check_overdispersion(floral_abundance_glmmTMB)
    Anova(floral_abundance_glmmTMB, type="III")
    
    #--FRUIT ABUNDANCE--#
    fruit_abundance_glmmTMB <- glmmTMB( fruit ~  precipitation * species +  (1|block/plot) , family=nbinom2,  data= fruit_abundance_data) 
    check_overdispersion(fruit_abundance_glmmTMB)
    Anova(fruit_abundance_glmmTMB, type="III")
    
    #--DURATION FRUIT--#
    duration_fruit_glmmTMB <- glmmTMB(duration_fruit ~  precipitation * species +  (1|block/plot)  , family=nbinom2, data= phenology_data) 
    check_overdispersion(duration_fruit_glmmTMB)
    summary(duration_fruit_glmmTMB)
    Anova(duration_fruit_glmmTMB, type="III")
    #summary(duration_fruit_glmmTMB)$coefficients
    plotResiduals(duration_fruit_glmmTMB)
    
    #--SEED VIABILITY--#
    proportion_viable_glmmTMB <- glmmTMB(proportion_viable ~  precipitation * species + (1|block/plot) ,family=betabinomial, weights = total_seeds_sum, data= dat_seed_data) 
    Anova(proportion_viable_glmmTMB, type="III")

#-------------- L I F E   H I S T O R Y  -  L E V E L --------------#
     #--PEAK FLOWER--#     
      peakflower_life_history_glmmTMB <- glmmTMB(peak_flower ~  precipitation * life_history +  (1|species) + (1|block/plot),  family= poisson(link=log), data = phenology_data)
      check_overdispersion(peakflower_life_history_glmmTMB)
      summary(peakflower_life_history_glmmTMB)
      Anova(peakflower_life_history_glmmTMB, type="III")
 
        #phenological shift
       peakflower_annual_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|species)+(1|block/plot),family=poisson(link=log),data=subset(phenology_data,life_history=="annual"))
       summary(peakflower_annual_glmmTMB)
       #phenology shift calculation:
       exp(5.3090638 + 0.0004481*-60) - exp(5.3090638 + 0.0004481*0)
       #[1] -5.362881
       
       peakflower_perennial_glmmTMB <- glmmTMB(peak_flower~precipitation+(1|species)+(1|plot),family=poisson(link=log),data=subset(phenology_data,life_history=="perennial"))
       summary(peakflower_perennial_glmmTMB)
       #phenology shift calculation:
       exp(5.489e+00 -1.284e-04*-60) - exp(5.489e+00 -1.284e-04*0)
       #[1] 1.871685

      #--FLOWER ABUNDANCE--#
      floral_abundance_life_history_glmmTMB <- glmmTMB( (flower) ~  scale(precipitation) * life_history +  (1|species) + (1|block/plot), family=nbinom2,   data= floral_abundance_data) #keep the interaction model
      check_overdispersion(floral_abundance_life_history_glmmTMB)
      as.glht(emmeans(floral_abundance_life_history_glmmTMB, poly ~ life_history))
      Anova(floral_abundance_life_history_glmmTMB, type="III")
      summary(floral_abundance_life_history_glmmTMB)$coefficients
      
      #--FRUIT ABUNDANCE--#
      fruit_abundance_life_history_glmmTMB <- glmmTMB( (fruit) ~  scale(precipitation) * life_history +  (1|species) + (1|block/plot), family=nbinom2,  data= fruit_abundance_data) #keep the interaction model
      check_overdispersion(fruit_abundance_life_history_glmmTMB)
      Anova(fruit_abundance_life_history_glmmTMB, type="III")
      
      #--DURATION FLOWER--#
      duration_flower_life_history_glmmTMB <- glmmTMB( (duration_flower) ~  scale(precipitation) * life_history +  (1|species) + (1|block/plot),  family=nbinom2, data = phenology_data)
      check_overdispersion(duration_flower_life_history_glmmTMB)
      #as.glht(emmeans(duration_flower_life_history_glmmTMB, poly ~ life_history))
      Anova(duration_flower_life_history_glmmTMB, type="III")
      #summary(duration_flower_life_history_glmmTMB)$coefficients
      
      #--DURATION FRUIT--#
      duration_fruit_life_history_glmmTMB <- glmmTMB((duration_fruit) ~  scale(precipitation) * life_history +  (1|species) + (1|block/plot),  family=nbinom2, data = phenology_data)
      check_overdispersion(duration_fruit_life_history_glmmTMB)
      #as.glht(emmeans(duration_fruit_life_history_glmmTMB, poly ~ life_history))
      Anova(duration_fruit_life_history_glmmTMB, type="III")
      #summary(duration_fruit_life_history_glmmTMB)$coefficients
      
      #--SEED VIABILITY--#
      proportion_viable_life_history_glmmTMB <- glmmTMB(proportion_viable ~  precipitation * duration +  (1|species) + (1|block/plot),  family= betabinomial, weights = total_seeds_sum, data = dat_seed_data)
      check_overdispersion(proportion_viable_life_history_glmmTMB)
      Anova(proportion_viable_life_history_glmmTMB, type="III")
  #-----------------------------------------end-----------------------------------------#


 