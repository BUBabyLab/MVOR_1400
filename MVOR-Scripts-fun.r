#  my GGPLOT for MVOR 
#  Written 10-8-18 (PG)
#
myGGplot <- function( theData) 
{
  theData <- theData[c(202:268, 1:67, 68:134, 135:201), ]
  theData$AOI <- factor(theData$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))

  p <- ggplot(data=theData, aes(x=TimeBin, y=MeanProp, colour=AOI, ymin = MeanProp - SEProp, ymax = MeanProp + SEProp)) + 
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
  print (p)  # To force ggplot to show graph in R-Studio output window
}

myFillerGGplot <- function (theData)
{
  p <- ggplot(data=FixProp_Filler, aes(x=TimeBin, y=MeanProp, colour=AOI_Filler, ymin = MeanProp - SEProp, ymax = MeanProp + SEProp)) + 
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
  print (p)  # To force ggplot to show graph in R-Studio output window
}

myCluster_LMER_Compare <- function (ClusterDF)
{
    ClusterDF$TargetRef<- relevel(ClusterDF$AOI, ref = "Target")
  
    TCData <- make_time_cluster_data(data = ClusterDF, 
                                    predictor_column = "AOI",
                                    test = "lmer", 
                                    threshold = 2,
                                    treatment_level = "Target",
                                    formula = Prop~AOI+(1 | Subject)
                                    #aoi = c("Target"), 
                                    #p_adjust_method = "holm",
    )
  
    summary(TCData)
  
    p <- plot(TCData, type = "estimate") + theme_light()
    print (p)  # To force ggplot to show graph in R-Studio output window
    
    tb_bootstrap <- analyze_time_clusters(TCData, 
                                        samples = 25,
                                        within_subj = TRUE
    )
  
    p <- plot(tb_bootstrap) + theme_light()
    print (p)  # To force ggplot to show graph in R-Studio output window
    
    summary(tb_bootstrap)
}

myClusterTTest <- function (ClusterDF)
{
  # visualize timecourse
  p <- plot(ClusterDF, predictor_column = "AOI") + theme_light() +
      coord_cartesian(ylim = c(0,1))
  print (p)  # To force ggplot to show graph in R-Studio output window
  
  ClusterDF$TargetRef<- relevel(ClusterDF$AOI, ref = "Target")
  #sink("myRfile.txt", append=FALSE, split=FALSE)
  
  Time_Cluster_Data_Kids_TvsC <- make_time_cluster_data(ClusterDF, 
                                                        test = "t.test", paired = TRUE,
                                                        predictor_column = "AOI",
                                                        threshold = 2.06#,
                                                        #formula = Prop~AOI+(1 | Subject)
  )
  summary(Time_Cluster_Data_Kids_TvsC)
  
  p<- plot(Time_Cluster_Data_Kids_TvsC, type = "estimate") + theme_light()
  print (p)  # To force ggplot to show graph in R-Studio output window
  
  tb_bootstrap <- analyze_time_clusters(Time_Cluster_Data_Kids_TvsC, within_subj=TRUE, paired=TRUE,
                                        samples = 250)
  
  p <- plot(tb_bootstrap) + theme_light()
  print (p)  # To force ggplot to show graph in R-Studio output window
  
  summary(tb_bootstrap)
}  #  End T-Test Cluster Analysis

#Not DONE
myCluster_LM_Test <- function (ClusterDF)
{
  ClusterDF$TargetRef<- relevel(ClusterDF$AOI, ref = "Target")
  
  response_time_between <- make_time_sequence_data(response_window_clean,
                                                 time_bin_size = 25, 
                                                 predictor_columns = ClusterDF$AOI,
                                                 aois = ClusterDF$AOI,
                                                 summarize_by = aois )

  df_timeclust_between <- make_time_cluster_data(response_time_between, 
                                               test= "lm",
                                               predictor_column = "AOI", 
                                               threshold = threshold_t) 
  plot(df_timeclust_between) +  ylab("T-Statistic") + theme_light()
}
