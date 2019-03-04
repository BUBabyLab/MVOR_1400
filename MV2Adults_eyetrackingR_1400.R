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
library(eyetrackingR)
library(pbapply)

#  Make one of these for the baby lab
source("MVOR-Scripts-fun.r") # PG-defined functions for script

######################################################
######################################################
######################################################

######            Data Analysis             !########

######################################################
######################################################
######################################################

# First get files and combine
MV2_SMP_FILLER <- read.csv(file = "MV2_SMP_FILLER_TSTAMP_ADULTS.csv", header = T)
MV2_SMP_EXP <- read.csv(file = "MVOR2_1400ms_Adults_TSTAMP.csv", header = T)
MV2A_SMP = rbind(MV2_SMP_EXP, MV2_SMP_FILLER)

###Creating a column that specifies whether the eye tracker lost the eye for a given sample
MV2A_SMP$SampleLost <- 999

#Attempt at lapply use in this setup
#MV2A_SMP$SampleLost <- lapply( MV2A_SMP, function (x) { (MV2A_SMP$L.Raw.X..px. == 0) & (MV2A_SMP$L.POR.X..px.==0)} )

#REPLACE w. dplyr call
for (i in 1:nrow(MV2A_SMP)) {
  #Target
  if ((MV2A_SMP[i, c("L.POR.X..px.")] == 0) & (MV2A_SMP[i, c("L.POR.Y..px.")] == 0)) {
    MV2A_SMP[i, c("SampleLost")] <- 1
  } else {
    MV2A_SMP[i, c("SampleLost")] <- 0
  }
  if ((i %% 1000) == 0) {
    print(i)
  }
}

write.csv(MV2A_SMP, file = "MV2A_EyetrackingR_Adults_1400s.csv")
#Rename "Time" to TimeSMI
#Rename TIME to TimeMS

#read in file for experimental trials
MV2A <- read.csv(file = "MV2A_EyetrackingR_Adults_1400s.csv", header=TRUE)
colnames(MV2A)
head(MV2A, 10)

youngKids = filter (MV2A, AgeYrs == 3)
middleKids = filter (MV2A, AgeYrs == 4)
oldKids = filter (MV2A, AgeYrs == 5)

####create file for eyetrackingR
MV2A_ETR_Adults_1400ms <- make_eyetrackingr_data(MV2A, 
                               participant_column = "Subject",
                               trial_column = "Trial",
                               time_column = "SampleTime",
                               trackloss_column = "SampleLost",
                               aoi_columns = c('Target','Complement', 'DiffExemplar', 'Distractor'),
                               treat_non_aoi_looks_as_missing = TRUE)
head(M2A_ETR_Adults_1400ms)

###########################################################################################
#######Look at track loss, and remove participants and trials with too much track loss#####
###########################################################################################
#Create a file removing any time points after 1600ms 
MV2A_ETR_1600 <- subset_by_window(MV2A_ETR_Adults_1400ms,
                                        window_start_time = 0, 
                                        window_end_time = 1600, 
                                        rezero = TRUE, 
                                        remove = TRUE)

#How many subjects
length(unique(MV2A_ETR_1600[, c('Subject')]))
unique(MV2A_ETR_1600$Subject)

###Look at track loss data (data where eyetracker lost the eye)
trackloss <- trackloss_analysis(data = MV2A_ETR_1600)
trackloss

write.csv(trackloss, file = "MV2A_EyetrackingR_Adults_1400ms_TRACKLOSS.csv")

#Check for looks into each AOI (low is not necessarily a problem; justs suggests overall prop of looks)
table(MV2A_ETR_1600$Target)
table(MV2A_ETR_1600$Complement)
table(MV2A_ETR_1600$DiffExemplar)
table(MV2A_ETR_1600$Distractor)

###Remove subjects and trials with over a proportion certain proportion of track loss
MV2A_ETR_1600_Clean <- clean_by_trackloss(MV2A_ETR_1600, participant_prop_thresh = 0.25, 
                                          trial_prop_thresh = 0.25)
#MV2A_ETR_1600_Clean <- clean_by_trackloss(MV2A_ETR_1600, trial_prop_thresh = 0.6)
unique(MV2A_ETR_1600_Clean$Subject)
write.csv(MV2A_ETR_1600_Clean, file = "MV2A_EyetrackingR_Adults_1400ms_CLEANED.csv")

###get the proportion of track loss per subject
trackloss_clean <- trackloss_analysis(data = MV2A_ETR_1600_Clean)
trackloss_clean_subjects <- unique(trackloss_clean[, c('Subject','TracklossForParticipant')])

# get mean proportion of samples contributed per trials, with SD
mean(1 - trackloss_clean_subjects$TracklossForParticipant)
sd(1- trackloss_clean_subjects$TracklossForParticipant)

#head ( MV2A_ETR_1600_Clean, 3)

# look at the number of trials contributed by each subject
final_summary <- describe_data(MV2A_ETR_1600_Clean, 'Target', 'Subject')
final_summary <- describe_data(MV2A_ETR_1600_Clean, 'Complement', 'Subject')

#  mySummary <- MV2A_ETR_1600_Clean %>% 
#  group_by(Subject, Condition, Trial) %>% 
#  summarise(MaxSampleTime=max(SampleTime), MaxTimeStamp=max(TimeStamp))

mySummary <- ddply(MV2A_ETR_1600_Clean, c("Subject", "Condition", "Trial"), summarise, 
                   MaxTime=max(SampleTime), SEMaxTime= sd(MaxTime)/sqrt(length(MaxTime)) )

mySummary <- ddply(MV2A_ETR_1600_Clean, c("Subject", "Condition", "Trial"), 
                   summarise, MaxTStamp=max(TimeStamp), SEMaxTStamp=sd(MaxTStamp)/sqrt(length(MaxTStamp)))

write.csv(mySummary, file = "SampleTimeSummary.csv")

########################################################3
###Descriptives and figure time!!!!!!!!!!!!!!###########
########################################################

###  NOTE!! ####
#  Need to run this for each level of Kids!  Or not - just run separately  #
MV2A_Exp_Clean <- MV2A_ETR_1600_Clean

###Plot proportion of looking to the target for parts and features, for each participant
MeanLooksToTarget <- describe_data(MV2A_Exp_Clean, 
                              describe_column='Target', group_columns=c('Condition','Subject'))
plot(MeanLooksToTarget)

###Plot proportion of looking to the compliment for parts and features, for each participant
MeanLooksToComplement <- describe_data(MV2A_Exp_Clean, 
                                   describe_column='Complement', group_columns=c('Condition','Subject'))
plot(MeanLooksToComplement)

###Plot proportion of looking to the different exemplar for parts and features, for each participant
MeanLooksToDiffExemplar <- describe_data(MV2A_Exp_Clean, 
                                       describe_column='DiffExemplar', group_columns=c('Condition','Subject'))
plot(MeanLooksToDiffExemplar)

###Plot proportion of looking to the distractor for parts and features, for each participant
MeanLooksToDistractor <- describe_data(MV2A_Exp_Clean, 
                                         describe_column='Distractor', group_columns=c('Condition','Subject'))
plot(MeanLooksToDistractor)



# Create proportion of looks in each AOI across features and parts, in time bins of 
LookProp <- make_time_sequence_data(MV2A_ETR_1600_Clean, time_bin_size = 25, 
                                         predictor_columns = c("Condition", "TrialType"),
                                         aois = c("Target", "Complement", "DiffExemplar", "Distractor")
)

###Cluster analysis, or hopefully plotting all AOIs together
LookProp$AOI<-as.factor(LookProp$AOI)

####Remove NaN's from proportion
LookProp_NoNA <- LookProp[LookProp$Prop != 'NaN', ]

###Create file with only Experimental Trials
LookProp_Exp <- LookProp_NoNA[LookProp_NoNA$TrialType == 'Experimental', ]

###Create file with only Filler Trials
LookProp_Filler <- LookProp_NoNA[LookProp_NoNA$TrialType == 'Filler', ]


###########################################
###########################################
##########Experimental Trials##############
########Figures and Analyses###############
###########################################
###########################################

###Fixation proportion###
#  CHANGE: Calc SD per Subject FIRST?
FixProp <- ddply(LookProp_Exp, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
#FixProp

###################################
###Parts x AOI Figure##############
###################################
FixProp_Parts <- FixProp[FixProp$Condition == 'Parts', ]

###Reorder bars###
FixProp_Parts <- FixProp_Parts[c(202:268, 1:67, 68:134, 135:201), ]
FixProp_Parts$AOI <- factor(FixProp_Parts$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))

myGGplot(FixProp_Parts)

ggplot(data=FixProp_Parts, aes(x=TimeBin, y=MeanProp, colour=AOI, ymin = MeanProp - SEProp, ymax = MeanProp + SEProp)) + 
  geom_smooth(stat="identity", position = "dodge", size=2.0)+
  xlab("Time-bin Number (25 ms bins)") + ylab("Proportion of Looks") + ggtitle("5-Yr-Olds: Parts Exp. Trials") +
  coord_cartesian(xlim = c(0, 66), ylim = c(0.0, 0.90)) +        
  scale_fill_brewer(palette="Dark2", name = "Orientation") +
  theme(plot.title = element_text(colour="Black", size=28, face="bold"), 
        legend.title = element_text(colour="Black", size=20, face="bold"), 
        legend.text = element_text(colour="Black", size=18, face="bold"),
        axis.title.x = element_text(colour="Black", size=20, face="bold"), 
        axis.title.y = element_text(colour="Black", size=20, face="bold"),
        axis.text.x = element_text(colour="Gray48", size=18, face="bold"),
        axis.text.y = element_text(colour="Gray48", size=18, face="bold"))

////// 
  
###################################
###Features x AOI Figure##############
###################################

FixProp_Features <- FixProp[FixProp$Condition == 'Features', ]

###Reorder bars###
FixProp_Features <- FixProp_Features[c(202:268, 1:67, 68:134, 135:201), ]
FixProp_Features$AOI <- factor(FixProp_Features$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))


ggplot(data=FixProp_Features, aes(x=TimeBin, y=MeanProp, colour=AOI, ymin = MeanProp - SEProp, ymax = MeanProp + SEProp)) + 
  geom_smooth(stat="identity", position= "dodge", size=2.0)+
  xlab("Time-bin Number (25 ms bins)") + ylab("Proportion of Looks") + ggtitle("5-Yr-Olds: Features Exp. Trls") +
  coord_cartesian(xlim = c(0, 66), ylim = c(0.0, 0.90)) +        
  scale_fill_brewer(palette="Dark2", name = "Orientation") +
  theme(plot.title = element_text(colour="Black", size=28, face="bold"), 
        legend.title = element_text(colour="Black", size=20, face="bold"), 
        legend.text = element_text(colour="Black", size=18, face="bold"),
        axis.title.x = element_text(colour="Black", size=20, face="bold"), 
        axis.title.y = element_text(colour="Black", size=20, face="bold"),
        axis.text.x = element_text(colour="Gray48", size=18, face="bold"),
        axis.text.y = element_text(colour="Gray48", size=18, face="bold"))


#############################
######LM proportion --
##   Time bin collapsed across trials
######  Stats ######
############################

###Fixation proportion###
Prop_Sub_Trial <- ddply(LookProp_Exp, c("Condition", "AOI", "Subject", "Trial"), summarise, TrialProp=mean(Prop), TrialSEProp=sd(Prop)/sqrt(length(Prop)))
Prop_Sub_Trial

####quick and dirty lm of average proportion in each AOI
Prop_Sub_Trial$DistractorRef<- relevel(Prop_Sub_Trial$AOI, ref = "Distractor")

LMProp <- lmer(TrialProp~ DistractorRef*Condition +(1|Subject), data=Prop_Sub_Trial)
summary(LMProp)

#Parts only
Prop_Sub_Trial_Parts <- Prop_Sub_Trial[Prop_Sub_Trial$Condition == 'Parts', ]

LMProp_Parts <- lmer(TrialProp~ DistractorRef +(1|Subject), data=Prop_Sub_Trial_Parts)
summary(LMProp_Parts)


#Features only
Prop_Sub_Trial_Features <- Prop_Sub_Trial[Prop_Sub_Trial$Condition == 'Features', ]

LMProp_Features <- lmer(TrialProp~ DistractorRef +(1|Subject), data=Prop_Sub_Trial_Features)
summary(LMProp_Features)



#################################################################################
###################Trying to do the bootstrap cluster analysis!!!!!!#############
#################################################################################
######################
######################
######Parts first#####
######################
######################

###Create file with only Experimental Trials
LookProp_Exp_Parts <- LookProp_Exp[LookProp_Exp$Condition == 'Parts', ]

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_TvsC <- LookProp_Exp_Parts[(LookProp_Exp_Parts$AOI == 'Target' | LookProp_Exp_Parts$AOI == 'Complement') , ]

LookProp_Exp_Parts_TvsC$TargetRef<- relevel(LookProp_Exp_Parts_TvsC$AOI, ref = "Target")

##################################

#   T-VALUE!! #
# Find our threshold t value 
# (2.03 or so - df determined) 
# Consider df and the overall approach!

#################################
threshold_t = qt(p = 1 - .05/2, 
                 df = length(unique(LookProp_Exp_Parts_TvsC$Subject)) - 1)

length(unique(LookProp_Exp_Parts_TvsC$Subject))
threshold_t  #see calc - will be used next

#Bootstrapped cluster-based permutation analysis
time_cluster_data_parts_TvsC <- make_time_cluster_data(data = LookProp_Exp_Parts_TvsC,
                                            predictor_column = "AOI",
                                            test = "lmer",
                                            threshold = threshold_t,
                                            treatment_level = 'Target',
                                            formula = Prop ~ AOI + (1|Subject) + (1|Trial) )


plot(time_cluster_data_parts_TvsC) +  ylab("T-Statistic") + theme_light()


clust_analysis_Parts_TvsC <- analyze_time_clusters(time_cluster_data_parts_TvsC, samples=100, within_subj=TRUE)
plot(clust_analysis_Parts_TvsC) + theme_light()
summary(clust_analysis_Parts_TvsC)

############################################
######## Complement vs, different exemplar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_CvsD <- LookProp_Exp_Parts[(LookProp_Exp_Parts$AOI == 'Complement' | LookProp_Exp_Parts$AOI == 'DiffExemplar') , ]

#find our threshold t value (2.03)
#threshold_t = qt(p = 1 - .05/2, 
#                 df = 35-1)


#Bootsrapped cluster-based permutation analysis

time_cluster_data_parts_CvDiff <- make_time_cluster_data(data = LookProp_Exp_Parts_CvsD,
                                            predictor_column = "AOI",
                                            test = "lmer",
                                            threshold = threshold_t,
                                            treatment_level = 'DiffExemplar',
                                            formula = Prop ~ AOI + (1|Subject) + (1|Trial) )


plot(time_cluster_data_parts_CvDiff) +  ylab("T-Statistic") + theme_light()
summary(time_cluster_data_parts_CvDiff)

clust_analysis_Parts_CvsDiff <- analyze_time_clusters(time_cluster_data_parts_CvDiff, samples=100, within_subj=TRUE)
plot(clust_analysis_Parts_CvsDiff) + theme_light()
summary(clust_analysis_Parts_CvsDiff)


############################################
######## distractor vs, different exemplar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_DvsDiff <- LookProp_Exp_Parts[(LookProp_Exp_Parts$AOI == 'Distractor' | LookProp_Exp_Parts$AOI == 'DiffExemplar') , ]

#find our threshold t value (2.03)
threshold_t = qt(p = 1 - .05/2, 
                 df = length(unique(LookProp_Exp_Parts_TvsC$Subject)) - 1)

#Bootsrapped cluster-based permutation analysis

time_cluster_data_parts_DvDiff <- make_time_cluster_data(data = LookProp_Exp_Parts_DvsDiff,
                                                         predictor_column = "AOI",
                                                         test = "lmer",
                                                         threshold = 2,
                                                         treatment_level = 'Distractor',
                                                         formula = Prop ~ AOI + (1|Subject) + (1|Trial) )


plot(time_cluster_data_parts_DvDiff) +  ylab("T-Statistic") + theme_light()
summary(time_cluster_data_parts_DvDiff)

clust_analysis_Parts_DvsDiff <- analyze_time_clusters(time_cluster_data_parts_DvDiff, samples=100, within_subj=TRUE)
plot(clust_analysis_Parts_DvsDiff) + theme_light()
summary(clust_analysis_Parts_DvsDiff)
######################
######################
######Features########
######################
######################

###Create file with only Experimental Trials
LookProp_Exp_Features <- LookProp_Exp[LookProp_Exp$Condition == 'Features', ]

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Features_TvsC <- LookProp_Exp_Features[(LookProp_Exp_Features$AOI == 'Target' | LookProp_Exp_Features$AOI == 'Complement') , ]

LookProp_Exp_Parts_TvsC$TargetRef<- relevel(LookProp_Exp_Parts_TvsC$AOI, ref = "Target")

#find our threshold t value (2.03)
threshold_t = qt(p = 1 - .05/2, 
                 df = 35-1)


#Bootsrapped cluster-based permutation analysis

time_cluster_data_features_TvsC <- make_time_cluster_data(data = LookProp_Exp_Features_TvsC,
                                                       predictor_column = "AOI",
                                                       test = "lmer",
                                                       threshold = 2,
                                                       treatment_level = 'Target',
                                                       formula = Prop ~ AOI + (1|Subject) + (1|Trial) )


plot(time_cluster_data_features_TvsC) +  ylab("T-Statistic") + theme_light()
summary(time_cluster_data_features_TvsC)


clust_analysis_Features_TvsC <- analyze_time_clusters(time_cluster_data_features_TvsC, samples=100, within_subj=TRUE)
plot(clust_analysis_Features_TvsC) + theme_light()
summary(clust_analysis_Features_TvsC)


############################################
######## Complement vs, different exemplar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Features_CvsD <- LookProp_Exp_Features[(LookProp_Exp_Features$AOI == 'Complement' | LookProp_Exp_Features$AOI == 'DiffExemplar') , ]

#find our threshold t value (2.03)
threshold_t = qt(p = 1 - .05/2, 
                 df = 35-1)


#Bootsrapped cluster-based permutation analysis

time_cluster_data_Features_CvDiff <- make_time_cluster_data(data = LookProp_Exp_Features_CvsD,
                                                         predictor_column = "AOI",
                                                         test = "lmer",
                                                         threshold = 2,
                                                         treatment_level = 'DiffExemplar',
                                                         formula = Prop ~ AOI + (1|Subject) + (1|Trial) )


plot(time_cluster_data_Features_CvDiff) +  ylab("T-Statistic") + theme_light()
summary(time_cluster_data_Features_CvDiff)

clust_analysis_Features_CvsDiff <- analyze_time_clusters(time_cluster_data_Features_CvDiff, samples=100, within_subj=TRUE)
plot(clust_analysis_Features_CvsDiff) + theme_light()
summary(clust_analysis_Features_CvsDiff)

############################################
######## distractor vs, different exemplar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Features_DvsDiff <- LookProp_Exp_Features[(LookProp_Exp_Features$AOI == 'Distractor' | LookProp_Exp_Features$AOI == 'DiffExemplar') , ]

#find our threshold t value (2.03)
threshold_t = qt(p = 1 - .05/2, 
                 df = 35-1)


#Bootsrapped cluster-based permutation analysis

time_cluster_data_features_DvDiff <- make_time_cluster_data(data = LookProp_Exp_Features_DvsDiff,
                                                            predictor_column = "AOI",
                                                            test = "lmer",
                                                            threshold = 2,
                                                            treatment_level = 'Distractor',
                                                            formula = Prop ~ AOI + (1|Subject) + (1|Trial) )


plot(time_cluster_data_features_DvDiff) +  ylab("T-Statistic") + theme_light()
summary(time_cluster_data_features_DvDiff)

clust_analysis_Features_DvsDiff <- analyze_time_clusters(time_cluster_data_features_DvDiff, samples=1000, within_subj=TRUE)
plot(clust_analysis_Features_DvsDiff) + theme_light()
summary(clust_analysis_Features_DvsDiff)


###########################################
###########################################
###########################################
##########Filler Trials####################
########Figures and Analyses###############
###########################################
###########################################
###########################################

###Make Variable so I can average across all the distractors
LookProp_Filler$AOI_Filler <- ifelse((LookProp_Filler$AOI == 'Target'), 'Target', 'Distractor')

###Fixation proportion###
FixProp_Filler <- ddply(LookProp_Filler, c("Condition", "AOI_Filler", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
#FixProp_Filler

###Reorder bars###
FixProp_Filler <- FixProp_Filler[c(68:134, 1:67, 202:268 , 135:201), ]
FixProp_Filler$AOI_Filler <- factor(FixProp_Filler$AOI_Filler, levels = c("Target", "Distractor"))

ggplot(data=FixProp_Filler, aes(x=TimeBin, y=MeanProp, colour=AOI_Filler, ymin = MeanProp - SEProp, ymax = MeanProp + SEProp)) + 
  geom_smooth(stat="identity",aes(linetype=Condition), position= "dodge", size=2.0)+
  xlab("Time-bin Number (25 ms bins)") + ylab("Proportion of Looks") + ggtitle("5-Yrs: Targ-Only Trials") +
  coord_cartesian(xlim = c(0, 66), ylim = c(0.0, 0.90)) +        
  scale_fill_brewer(palette="Dark2", name = "Orientation") +
  theme(plot.title = element_text(colour="Black", size=28, face="bold"), 
        legend.title = element_text(colour="Black", size=20, face="bold"), 
        legend.text = element_text(colour="Black", size=18, face="bold"),
        axis.title.x = element_text(colour="Black", size=20, face="bold"), 
        axis.title.y = element_text(colour="Black", size=20, face="bold"),
        axis.text.x = element_text(colour="Gray48", size=18, face="bold"),
        axis.text.y = element_text(colour="Gray48", size=18, face="bold"))


###############################
######LM proportion, time bin collapsed across trials

###Fixation proportion###
Prop_Sub_Trial_Filler <- ddply(LookProp_Filler, c("Condition", "AOI_Filler", "Subject", "Trial"), summarise, TrialProp=mean(Prop), TrialSEProp=sd(Prop)/sqrt(length(Prop)))
Prop_Sub_Trial_Filler

####quick and dirty lm of average proportion in each AOI
LMProp <- lmer(TrialProp~ AOI_Filler*Condition +(1|Subject), data=Prop_Sub_Trial_Filler)
summary(LMProp)

###Target only
Prop_Sub_Trial_Filler_Target <- Prop_Sub_Trial_Filler[Prop_Sub_Trial_Filler$AOI_Filler == 'Target', ]

LMProp <- lmer(TrialProp~ Condition +(1|Subject), data=Prop_Sub_Trial_Filler_Target)
summary(LMProp)

Prop_Filler <- ddply(LookProp_Filler, c("Condition", "AOI_Filler"), summarise, TrialProp=mean(Prop), TrialSEProp=sd(Prop)/sqrt(length(Prop)))
Prop_Filler

###Target only
Prop_Sub_Trial_Filler_Distractor <- Prop_Sub_Trial_Filler[Prop_Sub_Trial_Filler$AOI_Filler == 'Distractor', ]

LMProp <- lmer(TrialProp~ Condition +(1|Subject), data=Prop_Sub_Trial_Filler_Distractor)
summary(LMProp)



#########################################################
#########################################################
#############Trying to get average type of responses#####
#########################################################
#########################################################
colnames(MV2A_ETR_1600_Clean)
ForResponseTypeAverages <- unique(MV2A_ETR_1600_Clean[, c('Subject','Trial', 'Condition', 'TrialType', 'Array.ACC', 'Array.RESP')])

###Make Variable so I can average across all the distractors
ForResponseTypeAverages$InAnyQuad <- ifelse((ForResponseTypeAverages$Array.RESP == 'NoResponse'), 0, 1)

Average_Response <- ddply(ForResponseTypeAverages, c("TrialType", "Condition"), summarise, TargetResponse=mean(Array.ACC), TargetResponseSE=sd(Array.ACC)/sqrt(length(Array.ACC)), AnyQuadResponse=mean(InAnyQuad), AnyQuadResponseSE=sd(InAnyQuad)/sqrt(length(InAnyQuad)))
Average_Response

#LabelPhase_Age <- LabelPhase_Age[c(2,1,3,5,4,6,8,7,9,11,10,12,14,13,15,17,16), ]
#LabelPhase_Age$Age <- factor(LabelPhase_Age$Age, levels = c("24m", "30m", "36m"))
#LabelPhase_Age$Label_Phase <- factor(LabelPhase_Age$Label_Phase, levels = c("Demo", "Test", "Post Test"))

Average_Response$AgeYr<-as.factor(Average_Response$AgeYr)

####Any response
ggplot(data=Average_Response, aes(x=TrialType, y=AnyQuadResponse, fill=Condition, ymin = AnyQuadResponse - AnyQuadResponseSE, ymax = AnyQuadResponse + AnyQuadResponseSE)) +
  geom_bar(stat="identity", position= "dodge") + 
  geom_errorbar(position=position_dodge(.9), width=.4, colour = 'black') +
  xlab("Condition") + ylab("Proportion who looked any quadrant for allotted time") +
  coord_cartesian(xlim = NULL, ylim = c(0,0.8)) +        
  scale_fill_manual(values=c("#8c96c6", "#8856a7", "#810f7c"), name="Age") +
  theme(strip.text = element_text(colour="Black", size=20, face="bold"), 
        legend.title = element_text(colour="Black", size=20, face="bold"), 
        legend.text = element_text(colour="Black", size=18, face="bold"),
        axis.title.x = element_text(colour="Black", size=20, face="bold"), 
        axis.title.y = element_text(colour="Black", size=14, face="bold"),
        axis.text.x = element_text(colour="Gray48", size=18, face="bold"),
        axis.text.y = element_text(colour="Gray48", size=18, face="bold"))


####Target response
ggplot(data=Average_Response, aes(x=TrialType, y=TargetResponse, fill=Condition, ymin = TargetResponse - TargetResponseSE, ymax = TargetResponse + TargetResponseSE)) +
  geom_bar(stat="identity", position= "dodge") + 
  geom_errorbar(position=position_dodge(.9), width=.4, colour = 'black') +
  xlab("Condition") + ylab("Proportion who looked at target for allotted time") +
  coord_cartesian(xlim = NULL, ylim = c(0,0.8)) +        
  scale_fill_manual(values=c("#8c96c6", "#8856a7", "#810f7c"), name="Age") +
  theme(strip.text = element_text(colour="Black", size=20, face="bold"), 
        legend.title = element_text(colour="Black", size=20, face="bold"), 
        legend.text = element_text(colour="Black", size=18, face="bold"),
        axis.title.x = element_text(colour="Black", size=20, face="bold"), 
        axis.title.y = element_text(colour="Black", size=14, face="bold"),
        axis.text.x = element_text(colour="Gray48", size=18, face="bold"),
        axis.text.y = element_text(colour="Gray48", size=18, face="bold"))




