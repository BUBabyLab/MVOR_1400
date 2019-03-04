##  Divergence analysis from EyeTrackingR
#  PG - Oct, 2018

library(reshape2)
library(gtools) 
library(lme4)
library(plyr)
library(lattice)
library(car)
library(HH)
library(agricolae) 
library(multcomp) 
library(MASS)
library(heplots) 
library (coin) 
library(ggplot2)
library(Hmisc)
library(lmerTest)
library(LMERConvenienceFunctions)
library(languageR)
library(doBy)
library(RColorBrewer)
library("eyetrackingR")
library(pbapply)

#source ("eyeRTrackLossCreateColumn.r")
source("MVOR-Scripts-fun.r") # PG-defined functions for script

MV2A <- read.csv(file = "MV2A_EyetrackingR_Adults_1400s_FOR_DEV.csv", header=TRUE)
unique (MV2A$Subject)
unique (MV2A$Trial)
head(MV2A, 5)


####create file for eyetrackingR
MV2A_ETR_ADULT <- make_eyetrackingr_data(MV2A, 
                                         participant_column = "Subject",
                                         trial_column= "Trial",
                                         time_column = "SampleTime",
                                         trackloss_column = "SampleLost",
                                         aoi_columns = c('Target','Complement', 'DiffExemplar', 'Distractor'),
                                         treat_non_aoi_looks_as_missing = TRUE)

###########################################################################################
#######Look at track loss, and remove participants and trials with too much track loss#####
###########################################################################################
#Create a file removing any time points after 1600ms 
MV2A_ETR_1400 <- subset_by_window(MV2A_ETR_ADULT,
                                  window_start_time = 0, 
                                  window_end_time = 1600, 
                                  rezero = TRUE, 
                                  remove = TRUE)

#How many subjects
unique(MV2A_ETR_1400[, c('Subject')])

###Look at track loss data (data where eyetracker lost the eye)
trackloss <- trackloss_analysis(data = MV2A_ETR_1400)
trackloss

###Remove subjects and trials with over a proportion certain proportion of track loss
MV2A_Clean <- clean_by_trackloss(MV2A_ETR_1400, participant_prop_thresh = 0.75, 
                                          trial_prop_thresh = .35)

unique(MV2A_Clean$Subject)
########################################################
  ###     Descriptives and figure time!      ###########
########################################################
#
# Create proportion of looks across two AOIs for each of the features, parts and 
#   target-only ("fillers") conditions, in time bins of 25 msec
#

#
# For PARTS
#
TimeBins <- make_time_sequence_data(MV2A_Clean, 
                                    time_bin_size = 25, 
                                    predictor_columns = c("Condition", "TrialType"),
                                    aois = c("Target", "Complement", "DiffExemplar", "Distractor"),
                                    summarize_by = "Subject"
                                    )
colnames(TimeBins)
TimeBins$AOI<-as.factor(TimeBins$AOI)

####Remove NaN's from proportion
KTB_NoNA <- TimeBins[TimeBins$Prop != 'NaN', ]
KTB_NoNA

# How many NAs?
KTB_NAs <- TimeBins[TimeBins$Prop == 'NaN', ]
KTB_NAs

###Create df with only Experimental Trials
KTB_Exp <- KTB_NoNA[KTB_NoNA$TrialType == 'Experimental', ]

# Within that df, restric to Parts only
KTB_Exp_Parts <- KTB_Exp[KTB_Exp$Condition == 'Parts', ]

# Within that df, restric to Features only
KTB_Exp_Feat <- KTB_Exp[KTB_Exp$Condition == 'Features', ]

###Create file with only Filler Trials
KTB_Filler <- KTB_NoNA[KTB_NoNA$TrialType == 'Filler', ]

######################
##  PARTS:
#####################

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_TvsC <- KTB_Exp_Parts[(KTB_Exp_Parts$AOI == 'Target' | KTB_Exp_Parts$AOI == 'Complement') , ]
myClusterTTest (LookProp_Exp_Parts_TvsC)

############################################
######## Target vs, DiffExempar analysis####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_TvsDE <- KTB_Exp_Parts[(KTB_Exp_Parts$AOI == 'Target' | KTB_Exp_Parts$AOI == 'DiffExemplar') , ]
myCluster_LMER_Compare (LookProp_Exp_Parts_TvsDE)

LookProp_Exp_Parts_TvsDE <- KTB_Exp_Parts[(KTB_Exp_Parts$AOI == 'Target' | KTB_Exp_Parts$AOI == 'DiffExemplar') , ]
#colnames(LookProp_Exp_Parts_TvsDE)
#unique(LookProp_Exp_Parts_TvsDE$AOI)
myClusterTTest (LookProp_Exp_Parts_TvsDE)

############################################
######## Target vs. Distractor analysis#####
############################################

LookProp_Exp_Parts_TvsDist <- KTB_Exp_Parts[(KTB_Exp_Parts$AOI == 'Target' | KTB_Exp_Parts$AOI == 'Distractor') , ]
myClusterTTest (LookProp_Exp_Parts_TvsDist)
myClusterCompare (LookProp_Exp_Parts_TvsDist)

#########################
###  FEATURES
########################

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Feat_TvsC <- KTB_Exp_Feat[(KTB_Exp_Feat$AOI == 'Target' | KTB_Exp_Feat$AOI == 'Complement') , ]
myClusterTTest (LookProp_Exp_Feat_TvsC)

LookProp_Exp_Feat_TvsC$TargetRef<- relevel(LookProp_Exp_Feat_TvsC$AOI, ref = "Target")

############################################
######## Target vs, DiffExempar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Feat_TvsDE <- KTB_Exp_Feat[(KTB_Exp_Feat$AOI == 'Target' | KTB_Exp_Feat$AOI == 'DiffExemplar') , ]
myClusterCompare (LookProp_Exp_Feat_TvsDE)

#####################
FixProp <- ddply(KTB_Exp, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
#FixProp

###################################
###Parts x AOI Figure##############
###################################

FixProp_Parts <- FixProp[FixProp$Condition == 'Parts', ]

###Reorder bars###
FixProp_Parts <- FixProp_Parts[c(97:128, 1:32, 33:63, 65:96), ]
FixProp_Parts$AOI <- factor(FixProp_Parts$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))

ggplot(data=FixProp_Parts, aes(x=TimeBin+24, y=MeanProp, colour=AOI, ymin = MeanProp - SEProp, ymax = MeanProp + SEProp)) + 
  geom_smooth(stat="identity", position = "dodge", size=2.0)+
  xlab("Time-bin Number (25 ms bins)") + ylab("Proportion of Looks") + ggtitle("Adults: Parts Exp. Trials") +
  coord_cartesian(xlim = c(20, 60), ylim = c(0.0, 0.90)) +        
  scale_fill_brewer(palette="Dark2", name = "Orientation") +
  theme(plot.title = element_text(colour="Black", size=28, face="bold"), 
        legend.title = element_text(colour="Black", size=20, face="bold"), 
        legend.text = element_text(colour="Black", size=18, face="bold"),
        axis.title.x = element_text(colour="Black", size=20, face="bold"), 
        axis.title.y = element_text(colour="Black", size=20, face="bold"),
        axis.text.x = element_text(colour="Gray48", size=18, face="bold"),
        axis.text.y = element_text(colour="Gray48", size=18, face="bold"))
