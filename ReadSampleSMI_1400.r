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
  #library(multicore)
  library(tibble)
  library(parallel)
  
  install.packages('doParallel')
  library(doParallel)
  
  cores <- detectCores()
  cores
  
  ##########################################
  ##########################################
                                  
  ##  General Samples Cleaning script  ####
  ##  Opens all files in a folder       ####
  ##  Adds subject and Exp Name         ####
  ##  Alters time to start at zero      ####
  ##  July 17, 2018 - PG                ####

  ##########################################
  ##########################################

  #This prevents Scientific Notation RUN EVERY TIME!
  options(scipen=999)
  
  require(plyr)
  
  filenames <- list.files("MVOR2 2018 Subject Data/MVOR2 SMI", pattern="*.txt", full.names=TRUE )  ## Rename samples files as CSV or alter to TXT
  filenames
  
 # for (r in 1:length(filenames)) {
 #   myTest <- read.csv(filenames[r], header = T, skip = 38, sep = "\t")
 # }

    # OR
  
  #  NOTE!! 'comment.char = "#"' cannot be used (MSGs start w. #)
  
  system.time(ldf <- mclapply(filenames, read.csv, header = T, skip = 38, sep = "\t") )
  
  system.time(ldf <- lapply(filenames, read.csv, header = T, skip = 38, sep = "\t") )
  
  #  For Adults MVOR
  subjNums <- paste ("18299", substr(filenames, 35, 37), sep = "" ) # Need extra to match to E-Merge file
  filenames
  subjNums
  
  #  For Kids MVOR
  #subjNums <- paste ("1629", substr(filenames, 21, 23), sep = "" ) # Need extra to match to E-Merge file
  #ExpNames <- substr(filenames, 20, 25)
  
  
  #head (ldf, 3) # Prints a set of rows from first 3 data frames in list
  #head( ldf[[1]], 5)
    
  cat ("Number of Files Opened = ", length(ldf)()

  for (i in 1:length(ldf)) {
    
    print ("Start of Subj ")
    print (subjNums[i])
    
    #Remove all practice and SMI excess trials
    ldf[[i]] <- ldf[[i]] [(ldf[[i]]$Trial > 4 & ldf[[i]]$Trial < 53), ]
    
    # Drop all unnecessary columns from SMI files
    ldf[[i]] <- ldf[[i]] [ , c("Time", "Type", "Trial", "L.Raw.X..px.", 
                          "L.POR.X..px.", "L.POR.Y..px.", "Latency" )]
    
    ###  CONVERT TO MSEC & Add time index ## 
    ldf[[i]]$Time <- ldf[[i]]$Time/1000
    startTime <- ldf[[i]] [1, ("Time")]
    ldf[[i]]$Time <- ldf[[i]]$Time - startTime
    #ldf[[i]] <- myTimeConvert (ldf[[i]])  #  Function to call lapply with (does not work)
    
    #Figure out how to rename the columns to RawX, X and Y!!!
    
    ldf[[i]]$UseSample <- 0
    #ldf[[i]]$ExperimentName <- ExpNames[1]
    ldf[[i]]$Subject <- subjNums[i]  
   # add_column(ldf[[i]], Subject = 1919, .before = 1)
   # add_column(ldf[[i]], ExperimentName = "Test", .before = 1)
    
    #  Add flag for "in-stimulus-array" marker
    for (r in 2:nrow(ldf[[i]])) {
      k= r-1
      if ((ldf[[i]] [k, ("L.Raw.X..px.")]== '# Message: Correct') | (ldf[[i]] [k, ("L.Raw.X..px.")]== '# Message: Incorrect')) {
        ldf[[i]] [r, ("UseSample")] <- 0
        #print("Msg encountered")
      } else if (ldf[[i]] [ r, ("L.Raw.X..px.")] == '# Message: ArrayOnset') {
        ldf[[i]] [r, ("UseSample")] <- 1
        #print("Array Onset")
      } else {
        ldf[[i]] [r, ("UseSample")] <- ldf[[i]] [k, ("UseSample")]
      }
      
      if ((r %% 5000) == 0) { 
        cat ("Iteration = ", r)
      }
    }
    
    # Eliminate the non-array samples
    ldf[[i]] <- ldf[[i]] [(ldf[[i]]$UseSample == 1), ]
  }

  fulldf <- ldf[[1]]
  for (i in 2:length(ldf)) { fulldf <- rbind(fulldf, ldf[[i]]) }
  write.csv(fulldf, "MVOR2_1400ms_FULL_Samples_Adults.csv")
  
  #Print number of samples per file for check on over-run
  for (i in 1:length(ldf)) { 
    print (subjNums[i]) 
    print (nrow (ldf[[i]])) 
  }