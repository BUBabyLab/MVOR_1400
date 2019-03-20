setwd("/Users/laurieazoulai/Desktop/MVOR")

library(reshape2)
library(gtools) 
library(lme4)
library(plyr)
library(dplyr)
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
install.packages('doParallel')
library(doParallel)

library(doMC)
library(doSNOW)

registerDoSNOW(makeCluster(4, type = "SOCK"))

doMC::registerDoMC(cores=4) # or however many cores you have access to
cores <- detectCores()
cores

########################################
###
###   Load E-Prime Data for Trials   ###
###
########################################
#
#   Read in Samples Master File
#   Made by MVOR_MakeMaster_Samples.r
#   S. Olsen
##  Modified PG - June 29, 2018
#   Meshes SMI with EMerge files
#   Builds AOI data for EyeTrackingR
##
########################################

#This prevents Scientific Notation from being invoked
#  Should not be an issue - SMI merging script recalculates times
options(scipen=999)

# STOP!! Open and remove random "x" column here!
#  Delete "UseSample column as well unless done in prev. script
Samples <- read.csv(file = "MVOR2_1400ms_Full_Samples_Adults.csv", header=TRUE)

#unique(Samples$ExperimentName)
unique(Samples$Trial)
unique(Samples$Type)
unique(Samples$Subject)
nrow(Samples)

### Eprime File ###
### Read in the csv file 
#EPrime <- read.csv(file = "MVOR Adults EPrime Data Fall 2017 TRIMMED.csv", header=TRUE)
EPrime <- read.csv(file = "MVOR2 Adults Data Fall 2018 1400ms.csv", header=TRUE)
nrow(EPrime)

colnames(EPrime)

## For KIDS ONLY- AGE RANGE ANALYSIS Load age qualifiers and merge into SMI file
AgeList <- read.csv(file = "Age_List.csv", header = TRUE)
Samples <- merge(AgeList, Samples, by=c("Subject"), all=T)

###Reduce the number of variables
EPrime <- EPrime[ , c("ExperimentName", "Subject", "Trial", "ArrayFile", "ComplementFile", "ComplementLocation",
                      "Condition", "DiffExemplarFile", "DiffExemplarLocation", "DistractorFile", "DistractorLocation", "FvPMatch",
                      "ItemNumber", "ListNumber", "Q1File", "Q1Item", "Q2File", "Q2Item", "Q3File", "Q3Item", "Q4File", "Q4Item", 
                      "TargetFile", "TargetLocation", "TrialType"   )]


#Remove practice trials
EPrime <- EPrime[EPrime$Trial > 4, ]

unique(EPrime$ExperimentName)
unique(EPrime$Subject)
unique(Samples$Subject)
unique(EPrime$Trial)
nrow(EPrime)

###Merge the eye tracking and EPrime files###
MVOR2_MERG <- merge(EPrime, Samples, by=c("Subject", "Trial"), all=T)
nrow(MVOR2_MERG)
unique(MVOR2_MERG$Subject)
unique(MVOR2_MERG$ExperimentName)
nrow(MVOR2_MERG)

# Sort MVOR2 (using 'arrange' from ddply)
newData <- arrange(MVOR2_MERG, MVOR2_MERG$Subject, MVOR2_MERG$Trial, MVOR2_MERG$Time)
MVOR2_MERG <- newData
write.csv(MVOR2_MERG, file = "MVOR2_1400ms_Adults_Full_EPRIME-SMI.csv")

###########################################

##########Experimental Trials##############
########AOIS and Timestamp Columns#########

###########################################

##Read in the file


unique(MVOR2_MERG$Subject)
unique(MVOR2_MERG$ExperimentName)
unique(MVOR2_MERG$TrialType)

#create file with only messages
MSG_Only <- MVOR2_MERG[MVOR2_MERG$Type == 'MSG', ]
head(MSG_Only, n=15)

#create file with only file name messages (timestamp for start time of array presentation)
MSG_STIM <- MSG_Only[(MSG_Only$L.Raw.X..px. != '# Message: ArrayOnset') & (MSG_Only$L.Raw.X..px. != '# Message: Correct') & (MSG_Only$L.Raw.X..px. != '# Message: Incorrect') , ]

#remove unnessesary variables
MSG_STIM <- MSG_STIM[ , c("Subject","Trial", "Time", "L.Raw.X..px.")]
colnames(MSG_STIM)

###Rename variables
colnames(MSG_STIM) <- c("Subject","Trial", "ArrayOnset", "FileMsg")  # Change to FileMsg
colnames(MSG_STIM)

###Merge files
MV2_FULL <- merge(MVOR2_MERG, MSG_STIM, by=c("Subject", "Trial"), all=T)

#################################
### STOP AND SORT!! #############
#################################
# 'order' is the quick-sort function
newData2 <- MV2_FULL[
  order( MV2_FULL[,"Subject"], MV2_FULL[,"Trial"] ),
  ]

MV2_FULL <- newData2

#Invoke only if needed
write.csv(MV2_FULL, file = "MV2_Full_MSGSTIM_Sorted.csv")

unique (MV2_FULL$ExperimentName)
unique (MV2_FULL$Subject)
unique (MV2_FULL$TrialType)

###remove messages from data set
MV2_SMP <- MV2_FULL[MV2_FULL$Type != 'MSG', ]
unique (MV2_FULL$TrialType)
unique (MV2_SMP$TrialType)
write.csv(MV2_SMP, file = "MV2_SMP_Only.csv")

##################################

##  ADD AOIs                    ##
##  Might be able to            ##
##  Use Eyetracking R function  ##

##################################

#######Formula for marking 1 if the sample is in the AOI, and 0 if it is not.######
# First, create file with ONLY experimental trials
MV2_SMP_EXP <- MV2_SMP[MV2_SMP$TrialType == "Experimental", ]

##
##  Quandrant ID = QuadID
##  Sarah used binary flags - use 1-4 instead?  Faster...
##  Quad indicates location of the look
##  Set to 999 first to indicate NO AOI
##
MV2_SMP_EXP$Target <- 0   #  These first 4 are the AOI columns - set to 1 if present, 0 otherwise
MV2_SMP_EXP$Complement <- 0
MV2_SMP_EXP$DiffExemplar <- 0
MV2_SMP_EXP$Distractor <- 0
MV2_SMP_EXP$InQuad <- 0

#  When doing a data check:
#  DO NOT Check with all rows - change for a test to about 1000!  (Takes a long time)
#  Future - convert to lapply
for (i in 1:nrow(MV2_SMP_EXP)) { 
  x = MV2_SMP_EXP[i, c("L.POR.X..px.")]
  y = MV2_SMP_EXP[i, c("L.POR.Y..px.")]
  quad = 999
  
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 0) & (y <= 425)) ) { quad = 1 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 0) & (y <= 425)) ) {  quad = 2 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 655) & (y <= 1080)) ) { quad = 3 }
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 655) & (y <= 1080)) ) { quad = 4 }
  
  MV2_SMP_EXP[i,"InQuad"] <- quad  # Write value back into dataframe!
  
  if (quad == 1) {
    if (MV2_SMP_EXP[i,"TargetLocation"] == 'Q1') { MV2_SMP_EXP[i,"Target"] = 1}
    if (MV2_SMP_EXP[i,"ComplementLocation"] == 'Q1') { MV2_SMP_EXP[i,"Complement"] = 1}
    if (MV2_SMP_EXP[i,"DiffExemplarLocation"] == 'Q1') { MV2_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_SMP_EXP[i,"DistractorLocation"] == 'Q1') { MV2_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if (quad == 2) {
    if (MV2_SMP_EXP[i,"TargetLocation"] == 'Q2') { MV2_SMP_EXP[i,"Target"] = 1}
    if (MV2_SMP_EXP[i,"ComplementLocation"] == 'Q2') { MV2_SMP_EXP[i,"Complement"] = 1}
    if (MV2_SMP_EXP[i,"DiffExemplarLocation"] == 'Q2') { MV2_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_SMP_EXP[i,"DistractorLocation"] == 'Q2') { MV2_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if (quad == 3) {
    if (MV2_SMP_EXP[i,"TargetLocation"] == 'Q3') { MV2_SMP_EXP[i,"Target"] = 1}
    if (MV2_SMP_EXP[i,"ComplementLocation"] == 'Q3') { MV2_SMP_EXP[i,"Complement"] = 1}
    if (MV2_SMP_EXP[i,"DiffExemplarLocation"] == 'Q3') { MV2_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_SMP_EXP[i,"DistractorLocation"] == 'Q3') { MV2_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if (quad == 4) {
    if (MV2_SMP_EXP[i,"TargetLocation"] == 'Q4') { MV2_SMP_EXP[i,"Target"] = 1}
    if (MV2_SMP_EXP[i,"ComplementLocation"] == 'Q4') { MV2_SMP_EXP[i,"Complement"] = 1}
    if (MV2_SMP_EXP[i,"DiffExemplarLocation"] == 'Q4') { MV2_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_SMP_EXP[i,"DistractorLocation"] == 'Q4') { MV2_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if ((i %% 1000) == 0) {  print(i) }
}

unique(MV2_SMP_EXP$InQuad)
unique(MV2_SMP_EXP$TrialType)

write.csv(MV2_SMP_EXP,file="MVOR2_1400ms_Adults_SMP_EXP.csv")

#MV2_SMP_EXP<-read.csv(file="MVOR2_Test.csv",header=TRUE)

#######################################
#####Make time stamp number column#####
#######################################

###make column with time of first sample
MV2_SMP_EXP <- ddply(MV2_SMP_EXP, c("Subject", "Trial"), transform, TimeOfFirstSample=min(Time))

###Provide a time stamp once this command is done running!

###make column with sample time starting from the time of the first sample recorded after appearance of the array image
MV2_SMP_EXP$SampleTime <- (MV2_SMP_EXP$Time - MV2_SMP_EXP$TimeOfFirstSample) 

###Create a column (TimeStamp) which labels each incrementing sample time as a time stamp number starting with one, and increasing sequentially by one 
MV2_SMP_EXP <- ddply(MV2_SMP_EXP, c("Subject", "Trial"), transform, TimeStamp=seq(1, (length(SampleTime)), by=1))

###Write file with timestamp number###
##DO NOT REWRITE FILE!!!##
write.csv(MV2_SMP_EXP, file = "MVOR2_1400ms_Adults_TSTAMP.csv")



####we ran until here, 10/29/18, 5:40~####
###########################################
##
####        Filler Trials           #######
####  AOIS and Timestamp Columns    #######
##
###########################################
MV2_SMP <- read.csv(file = "MVOR2_1400ms_Adults_TSTAMP.csv", header=TRUE)

###remove messages from data set
MV2_SMP <- MV2_FULL[MV2_FULL$Type != 'MSG', ]

#Create file with only filler trials
MV2_SMP_FILLER <- MV2_SMP[MV2_SMP$TrialType == 'Filler', ]
unique (MV2_SMP_FILLER$TrialType)

##
##  Quandrant ID = QuadID
##  Sarah used binary flags - use 1-4 instead?  Faster...
##  Quad indicates location of the look
##  Set to 999 first to indicate NO AOI
##
MV2_SMP_FILLER$Target <- 0   #  These first 4 are the AOI columns - set to 1 if present, 0 otherwise
MV2_SMP_FILLER$Complement <- 0
MV2_SMP_FILLER$DiffExemplar <- 0
MV2_SMP_FILLER$Distractor <- 0
MV2_SMP_FILLER$InQuad <- 0

for (i in 1:nrow(MV2_SMP_FILLER)) { #  DO NOT Check with all rows - change for a test to about 1000!
  x = MV2_SMP_FILLER[i, c("L.POR.X..px.")]
  y = MV2_SMP_FILLER[i, c("L.POR.Y..px.")]
  quad = 999
  
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 0) & (y <= 425)) ) { quad = 1 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 0) & (y <= 425)) ) {  quad = 2 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 655) & (y <= 1080)) ) { quad = 3 }
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 655) & (y <= 1080)) ) { quad = 4 }
  
  MV2_SMP_FILLER[i,"InQuad"] <- quad  # Write value back into dataframe!
  
  if (quad == 1) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q1') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q1') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q1') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q1') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if (quad == 2) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q2') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q2') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q2') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q2') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if (quad == 3) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q3') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q3') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q3') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q3') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if (quad == 4) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q4') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q4') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q4') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q4') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if ((i %% 1000) == 0) {  print(i) }
}

unique(MV2_SMP_FILLER$InQuad)
unique(MV2_SMP_FILLER$TrialType)

#######################################
#####Make time stamp number column#####
#######################################

###make column with time of first sample
system.time(MV2_SMP_FILLER <- ddply(MV2_SMP_FILLER, c("Subject", "Trial"), transform, TimeOfFirstSample=min(Time)))

###make column with sample time starting from the time of the first sample recorded after appearance of the array image
MV2_SMP_FILLER$SampleTime <- (MV2_SMP_FILLER$Time - MV2_SMP_FILLER$TimeOfFirstSample) 

###Create a column (TimeStamp) which labels each incrementing sample time as a time stamp number starting with one, and increasing sequentially by one 
MV2_SMP_FILLER <- ddply(MV2_SMP_FILLER, c("Subject", "Trial"), transform, TimeStamp=seq(1, (length(SampleTime)), by=1))

###Write file with timestamp number###
##DO NOT REWRITE FILE!!!##
write.csv(MV2_SMP_FILLER, file = "MV2_SMP_FILLER_TSTAMP_ADULTS.csv")

MV2A_SMP = rbind(MV2_SMP_EXP, MV2_SMP_FILLER)

#  Next File is MV2Adults_eyetrackingR_PGv2.r