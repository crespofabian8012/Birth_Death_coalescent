#R code accompanying the paper "Coalescent Models Derived from
#Birth-Death Processes" by Crespo and Wiuf
#post_processing
rm(list=ls())
#####################################################################
#assemble_output
#Reads and assemble the output of the simulator ans produce a single data frame
#number_sim: number of simulations
#sample_size: sample size
#Delta_list: vector of Delta  values
#Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#path_scenario_A: path to a file containing the list of data frames generated for the simulator for scenario A 
#path_scenario_B: path to a file containing the list of data frames generated for the simulator for scenario B 
#path_scenario_C: path to a file containing the list of data frames generated for the simulator for scenario C 
#####################################################################
assemble_output<-function(number_sim, sample_size, Delta_list,Gamma_list,
                          path_scenario_A, path_scenario_B, path_scenario_C, path_scenario_hybrid, path_mean_times_origin,
                          include_scenario_C=FALSE, include_scenario_hybrid=FALSE, shifted_scale,nB
                          ){
  
  stopifnot(length(Delta_list)==length(Gamma_list))
  
  list_data_frames_A<-readRDS( path_scenario_A)
  list_data_frames_B<-readRDS( path_scenario_B)
  
  if (include_scenario_hybrid){
      list_data_frames_hybrid<-readRDS( path_scenario_hybrid)
  }
  
  if (include_scenario_C){
    list_data_frames_C<-readRDS( path_scenario_C)
  }
  
  
  delta_labs <- lapply(Delta_list, function(x, name) paste(name, "=", x, sep=""), "Delta")
  names(delta_labs)<-factor(Delta_list)
  gamma_labs <- lapply(Gamma_list, function(x, name) paste(name, "=", x, sep=""), "Gamma")
  names(gamma_labs)<- factor(Gamma_list)
  
  
  simulationColumn=rep(1:length(Delta_list), rep(sample_size, length(Delta_list)))
  
  
  Delta_column=rep(Delta_list, rep(sample_size-1, length(Delta_list)))
  
  Gamma_column=rep(Gamma_list, rep(sample_size-1, length(Gamma_list)))
  
  
  data_frame_A<-do.call("rbind", list_data_frames_A)
  data_frame_B<-do.call("rbind", list_data_frames_B)
  
  if (include_scenario_hybrid)
    data_frame_hybrid<-do.call("rbind", list_data_frames_hybrid)
  
  if (include_scenario_C){
    
    data_frame_C<-do.call("rbind", list_data_frames_C)
    
    data_frame_C$Delta=Delta_column
    data_frame_C$GammaVar=Gamma_column
    
    data_frame_C$x = rep(((sample_size-1):1), length(list_data_frames_hybrid))
    
    #data_frame_A$Gamma <- factor(data_frame_A$Gamma, levels = factor(Gamma_list[commonIndexes]),
       #                         labels = gamma_labs[commonIndexes])
    #data_frame_B$Gamma <- factor(data_frame_B$Gamma, levels = factor(Gamma_list[commonIndexes]),
      #                          labels = gamma_labs[commonIndexes])
    #data_frame_C$Gamma <- factor(data_frame_C$Gamma, levels = factor(Gamma_list[commonIndexes]),
       #                         labels = gamma_labs[commonIndexes])
    
  }
    #data_frame_A$Delta=factor(Delta_column)
    #data_frame_B$Delta=factor(Delta_column)
    #data_frame_A$Gamma=factor(Gamma_column)
    #data_frame_B$Gamma=factor(Gamma_column)
    
    data_frame_A$Delta=Delta_column
    data_frame_B$Delta=Delta_column
    
    if (include_scenario_hybrid)
       data_frame_hybrid$Delta=Delta_column
    
    data_frame_A$GammaVar=Gamma_column
    data_frame_B$GammaVar=Gamma_column
    
    if (include_scenario_hybrid)
       data_frame_hybrid$GammaVar=Gamma_column
    
    #data_frame_A$Gamma <- factor(data_frame_A$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
    #data_frame_B$Gamma <- factor(data_frame_B$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
    
    data_frame_A$x=rep(((sample_size-1):1), length(list_data_frames_A))
    
    if (include_scenario_hybrid)
       data_frame_hybrid$x=rep(((sample_size-1):1), length(list_data_frames_hybrid))
    
    interval_B<-(2*sample_size-2):sample_size
    if (shifted_scale)
        data_frame_B$x=rep(interval_B, length(list_data_frames_B))
    else
      data_frame_B$x=rep(((sample_size-1):1), length(list_data_frames_B))
  
  
  data_frame_A$Scenario="A"
  data_frame_B$Scenario="B"
  
  if (include_scenario_hybrid)
     data_frame_hybrid$Scenario=paste("H",nB, sep="")
  
  if (include_scenario_C)
    data_frame_C$Scenario="C"
  
  if (include_scenario_C & include_scenario_hybrid)
    combined_data_frame<-do.call("rbind", list(data_frame_A,data_frame_hybrid,  data_frame_B, data_frame_C))
  else if (include_scenario_hybrid)
    combined_data_frame<-do.call("rbind", list(data_frame_A, data_frame_hybrid, data_frame_B))
  else 
    combined_data_frame<-do.call("rbind", list(data_frame_A, data_frame_B))
  
  return(combined_data_frame)
  
}
#####################################################################
#get_ratio_data
#Reads and assemble the output of the simulator for scenarios A and B
#number_sim: number of simulations
#sample_size: sample size
#Delta_list: vector of Delta  values
#Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#path_scenario_A: path to a file containing the list of data frames generated for the simulator for scenario A 
#path_scenario_B: path to a file containing the list of data frames generated for the simulator for scenario B 
#####################################################################
get_ratio_data<-function(number_sim, sample_size, Delta_list,Gamma_list,
                         path_scenario_A, path_scenario_B){
  
  stopifnot(length(Delta_list)==length(Gamma_list))
  
  list_data_fames_A<-readRDS( path_scenario_A)
  list_data_fames_B<-readRDS( path_scenario_B)
 
  
  data_frame_A<-do.call("rbind", list_data_fames_A)
  data_frame_B<-do.call("rbind", list_data_fames_B)
  
  delta_labs <- lapply(Delta_list, function(x, name) paste(name, "=", x, sep=""), "Delta")
  names(delta_labs)<-factor(Delta_list)
  gamma_labs <- lapply(Gamma_list, function(x, name) paste(name, "=", x, sep=""), "Gamma")
  names(gamma_labs)<- factor(Gamma_list)
  
  
  simulationColumn=rep(1:length(Delta_list), rep(sample_size, length(Delta_list)))
  Delta_column=rep(Delta_list, rep(sample_size-1, length(Delta_list)))
  Gamma_column=rep(Gamma_list, rep(sample_size-1, length(Gamma_list)))
  
  # data_frame_A$Delta=factor(Delta_column)
  # data_frame_B$Delta=factor(Delta_column)
  # data_frame_A$Gamma=factor(Gamma_column)
  # data_frame_B$Gamma=factor(Gamma_column)
  
  data_frame_A$Delta=Delta_column
  data_frame_B$Delta=Delta_column
  data_frame_A$GammaVar=Gamma_column
  data_frame_B$GammaVar=Gamma_column
  
  #data_frame_A$Gamma <- factor(data_frame_A$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
  #data_frame_B$Gamma <- factor(data_frame_B$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
  
  data_frame_A$x=rep((sample_size-1):1, length(list_data_fames_A))
  data_frame_B$x=rep((sample_size-1):1, length(list_data_fames_B))
  
  mean_ratio_colum<-data_frame_A$mean /  data_frame_B$mean
  median_ratio_colum<-data_frame_A$median /  data_frame_B$median
  ratio_data_frame<-data.frame(RatioMeanA1overA2=mean_ratio_colum,RatioMedianA1overA2 =median_ratio_colum, Delta=Delta_column,
                               GammaVar=Gamma_column, x=rep((sample_size-1):1, length(list_data_fames_A)))
  
  #ratio_data_frame$Gamma <- factor(ratio_data_frame$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
  #ratio_data_frame$Delta <- factor(ratio_data_frame$Delta, levels = factor(Delta_list), labels=gamma_lbs)
  return(ratio_data_frame)
  
}
#####################################################################
#get_mean_time_origin_data
#Reads and assemble the output of the simulator for the mean time origin for scenarios 
#number_sim: number of simulations
#sample_size: sample size
#Delta_list: vector of Delta  values
#Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#path_mean_times_origin: path to a file containing the list of mean times of origins generated for the simulator for scenario A
#####################################################################
get_mean_time_origin_data<-function(number_sim, sample_size, Delta_list,Gamma_list,path_mean_times_origin){
  
  
  mean_times_origin<-readRDS(path_mean_times_origin)
  delta_labs <- lapply(Delta_list, function(x, name) paste(name, "=", x, sep=""), "Delta")
  names(delta_labs)<-factor(Delta_list)
  gamma_labs <- lapply(Gamma_list, function(x, name) paste(name, "=", x, sep=""), "Gamma")
  names(gamma_labs)<- factor(Gamma_list)
  
  mean_time_origin_df<-data.frame(Delta=factor(Delta_list), mean_time_origin= mean_times_origin, Gamma_list)
 # mean_time_origin_df<-data.frame(Delta=factor(Delta_list), mean_time_origin= mean_times_origin, Gamma=factor(Gamma_list))
#mean_time_origin_df$Gamma <- factor(mean_time_origin_df$Gamma, levels = factor(Gamma_list),  labels=gamma_lbs)
 # mean_time_origin_df$Delta <- factor(mean_time_origin_df$Delta, levels = factor(Delta_list), labels=gamma_lbs)
  
  return(mean_time_origin_df)
}
#####################################################################
#get_CI_time_origin_data
#Reads and assemble the output of the simulator for the mean time origin for scenarios 
#number_sim: number of simulations
#sample_size: sample size
#Delta_list: vector of Delta  values
#Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#path_CI_times_origin: path to a file containing the table of CI ot time of origins generated for the simulator for scenario A
#path_mean_times_origin: path to a file containing the list of mean times of origins generated for the simulator for scenario A
#####################################################################
get_CI_time_origin_data<-function(number_sim, sample_size, Delta_list,Gamma_list,path_CI_times_origin, path_mean_times_origin){
  
  mean_times_origin<-readRDS(path_mean_times_origin)
  mean_time_origin_df<-as.data.frame(readRDS(path_CI_times_origin))
  colnames(mean_time_origin_df) <- c("LI", "median", "UI")
  
  delta_labs <- lapply(Delta_list, function(x, name) paste(name, "=", x, sep=""), "Delta")
  names(delta_labs)<-factor(Delta_list)
  gamma_labs <- lapply(Gamma_list, function(x, name) paste(name, "=", x, sep=""), "Gamma")
  names(gamma_labs)<- factor(Gamma_list)
  
  mean_time_origin_df$Delta<-factor(Delta_list) 
  mean_time_origin_df$mean<- mean_times_origin
  mean_time_origin_df$GammaVar<-factor(Gamma_list)
  #mean_time_origin_df$Gamma <- factor(mean_time_origin_df$Gamma, levels = factor(Gamma_list),  labels=gamma_lbs)
  #mean_time_origin_df$Delta <- factor(mean_time_origin_df$Delta, levels = factor(Delta_list),  labels=gamma_lbs)
  
  return(mean_time_origin_df)
}
plot_comparison_Gumbel<-function(number_sim, Delta_list, path_to_save){
  
  list.sampled.Trans.Gumbel<-lapply(Delta_list, function(x, n) sampleTOriginFormEVD8(n, x,x), n=number_sim)
  means.sampled.Trans.Gumbel<-unlist(lapply(list.sampled.Trans.Gumbel, function(x) mean(x)))


  list.sampled.TOrigins<-lapply(Delta_list, function(x, n) drawRandomTimeOrigin2(n, x), n=number_sim)
  means.sampled.TOrigins<-unlist(lapply(list.sampled.Trans.Gumbel, function(x) mean(x)))

  comparisonTOriginDataFrame<- data.frame(logDelta=log(Delta_list),
                                          TOriginTransGumbel=means.sampled.Trans.Gumbel,
                                          TOrigin=means.sampled.TOrigins)
  
  p4<- ggplot(comparisonTOriginDataFrame, aes(logDelta))+
      geom_point(aes(y=TOriginTransGumbel), colour="red")+
      geom_point(aes(y=TOrigin), colour="green")+
      ggtitle(paste("Comparison of average times of origin vs Delta"))+
      labs(y="Time of origin", x ="log of Delta")+
      theme(plot.margin = unit(c(1,1,2,1), "lines")) +
      theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(paste(path_to_save, "/",  "ComparisonTransGumbel.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "ComparisonTransGumbel.png", sep=""))
  ggsave(paste(path_to_save, "/", "ComparisonTransGumbel.jpg", sep=""))
  
}
#####################################################################
#mean_ratio_waiting_times_A_respect_2th_coal_time
#build a list of ratios of mean waiting times per TMRCA per value pair (Delta, Gamma)
#number_sim: number of simulations
#sample_size: sample size
#Delta_list: vector of Delta  values
#Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#data_frame_A: the data frame for scenario A 
#####################################################################
mean_ratio_waiting_times_A_respect_2th_coal_time<-function(number_sim, sample_size, Delta_list,Gamma_list,
                                                      data_frame_A){
  
  stopifnot(length(Delta_list)==length(Gamma_list))
  
  
  
  list_result<-lapply(1:length(Gamma_list), function(i,data_frame_A, Gamma_list  ){
    
    dataFrame_for_gamma<-data_frame_A[data_frame_A$GammaVar==Gamma_list[i],]
    
    
    int_list<-lapply(1:(sample_size-1), function(k,dataFrame_for_gamma, Gamma_list ){
      if (k==(sample_size-1))
        last_term=0
      else{
        last_term<- dataFrame_for_gamma[dataFrame_for_gamma$x==(k+1),]$mean
      }
      return((dataFrame_for_gamma[dataFrame_for_gamma$x==k,]$mean-last_term)/ dataFrame_for_gamma[dataFrame_for_gamma$x==1,]$mean)
    }, dataFrame_for_gamma, Gamma_list)
    
    return(unlist(int_list))
  }, data_frame_A, Gamma_list)
  
  return(list_result)
}
#####################################################################
#mean_ratio_waiting_times_A_respect_time_origin
#build a list of ratios of mean waiting times per TMRCA per value pair (Delta, Gamma)
#number_sim: number of simulations
#sample_size: sample size
#Delta_list: vector of Delta  values
#Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#data_frame_A: the data frame for scenario A 
#meanTOriginDataFrame: the data frame with the mean TOrigin
#####################################################################
mean_ratio_waiting_times_A_respect_time_origin<-function(number_sim, sample_size, Delta_list,Gamma_list,
                                                        data_frame_A, meanTOriginDataFrame){
  
  stopifnot(length(Delta_list)==length(Gamma_list))
  
  
  
  list_result<-lapply(1:length(Gamma_list), function(i,data_frame_A, Gamma_list, meanTOriginDataFrame  ){
    
    dataFrame_for_gamma<-data_frame_A[data_frame_A$GammaVar==Gamma_list[i],]
    
    meanTorigin<-meanTOriginDataFrame[meanTOriginDataFrame$Gamma_list==Gamma_list[i],]$mean_time_origin
    
    first_term<-(meanTorigin -dataFrame_for_gamma[dataFrame_for_gamma$x==1,]$mean)/meanTorigin
    
    int_list<-lapply(1:(sample_size-1), function(k,dataFrame_for_gamma,  Gamma_list, meanTorigin ){
      if (k==(sample_size-1))
        last_term=0
      else{
        last_term<- dataFrame_for_gamma[dataFrame_for_gamma$x==(k+1),]$mean
      } 
     
      return((dataFrame_for_gamma[dataFrame_for_gamma$x==k,]$mean-last_term)/meanTorigin)
    }, dataFrame_for_gamma, Gamma_list, meanTorigin)
    
  
    return(c(first_term,unlist(int_list)))
  }, data_frame_A, Gamma_list, meanTOriginDataFrame)
  
  return(list_result)
}
getCurrentFileLocation <-  function()
{
  require(tibble)
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
#####################################################################
#compare scenarios A1(A), A2(B), A3(C)

sim=10000
sample_size=10
sample.size=sample_size
runA3Scenario=FALSE
doTestA1=TRUE
varyDelta=FALSE
runA3Scenario=FALSE
includeA3 =FALSE

Gamma_list=c(0.001,  0.1,  10)
Delta_list=rep(1,length(Gamma_list))
#####################################################################
require(ggplot2)
require(dplyr)
require(tidyr)

path_to_save = getCurrentFileLocation()

K=0.8
sample_size=100
nB=70
pathDataFramesA<-paste(path_to_save,"/", "listDataFramesA_",sample.size,"_", sim,".rds", sep="")
pathDataFramesB<-paste(path_to_save,"/","listDataFramesB_",sample.size,"_", sim,".rds", sep="")
pathDataFramesC<-paste(path_to_save,"/","listDataFramesC_",sample.size,"_", sim,".rds", sep="")
pathDataFramesTOrigin<-paste(path_to_save,"/", "meansToriginA_",sample_size,"_", sim,".rds", sep="")
pathDataFramesTOriginCI<-paste(path_to_save,"/", "listToriginCI_A_",sample_size,"_", sim,".rds", sep="")

pathDataFramesB_K<-paste(path_to_save,"/","listDataFramesB_",sample.size,"_K=",K, "_", sim,".rds", sep="")
pathDataFramesHybrid<-paste(path_to_save,"/","listDataFramesHybrid_",sample.size,"_",nB,"_", sim,".rds", sep="")  


all.data.frame_100_70<-assemble_output(sim, sample_size, Delta_list,Gamma_list,
                                path_scenario_A = pathDataFramesA,
                                path_scenario_B = pathDataFramesB,
                                path_scenario_hybrid=pathDataFramesHybrid,
                                path_scenario_C=pathDataFramesC , 
                                pathDataFramesTOrigin,
                          include_scenario_C=FALSE, include_scenario_hybrid=FALSE, shifted_scale=FALSE, nB )

ratio.data.frame_100_70<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                        pathDataFramesA, 
                                        pathDataFramesB)

ratio.data.frame_100_A_Hybrid_70<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                                 pathDataFramesA, 
                                                 pathDataFramesHybrid)

meanTOriginDataFrame_100_70<-get_mean_time_origin_data(sim, sample_size, Delta_list,Gamma_list,
                                                       pathDataFramesTOrigin)

CI_TOriginDataFrame_100_70<-get_CI_time_origin_data(sim, sample_size, Delta_list,Gamma_list,pathDataFramesTOriginCI,
                                                    pathDataFramesTOrigin)

nB=50
pathDataFramesHybrid<-paste(path_to_save,"/","listDataFramesHybrid_",sample.size,"_",nB,"_", sim,".rds", sep="") 

all.data.frame_100_50<-assemble_output(sim, sample_size, Delta_list,Gamma_list,
                                    path_scenario_A= pathDataFramesA,
                                    path_scenario_B=pathDataFramesB,
                                    path_scenario_hybrid=pathDataFramesHybrid,
                                    path_scenario_C=pathDataFramesC, 
                                    pathDataFramesTOrigin,
                                    include_scenario_C=FALSE,include_scenario_hybrid=TRUE, shifted_scale=FALSE, nB )

ratio.data.frame_100_50<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                        pathDataFramesA, 
                                        pathDataFramesB)

ratio.data.frame_100_A_Hybrid_50<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                                 pathDataFramesA, 
                                                 pathDataFramesHybrid)

meanTOriginDataFrame_100_50<-get_mean_time_origin_data(sim, sample_size, Delta_list,Gamma_list,
                                                       pathDataFramesTOrigin)

CI_TOriginDataFrame_100_50<-get_CI_time_origin_data(sim, sample_size, Delta_list,Gamma_list,pathDataFramesTOriginCI,
                                                    pathDataFramesTOrigin)

nB=90
pathDataFramesHybrid<-paste(path_to_save,"/","listDataFramesHybrid_",sample.size,"_",nB,"_", sim,".rds", sep="") 
all.data.frame_100_90<-assemble_output(sim, sample_size, Delta_list,Gamma_list,
                                       path_scenario_A= pathDataFramesA,
                                       path_scenario_B=pathDataFramesB,
                                       path_scenario_hybrid=pathDataFramesHybrid,
                                       path_scenario_C=pathDataFramesC, 
                                       pathDataFramesTOrigin,
                                       include_scenario_C=FALSE, include_scenario_hybrid=FALSE, shifted_scale=FALSE, nB)

ratio.data.frame_100_90<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                        pathDataFramesA, 
                                        pathDataFramesB)

ratio.data.frame_100_A_Hybrid_90<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                                 pathDataFramesA, 
                                                 pathDataFramesB)

meanTOriginDataFrame_100_90<-get_mean_time_origin_data(sim, sample_size, Delta_list,Gamma_list,
                                                       pathDataFramesTOrigin)

CI_TOriginDataFrame_100_90<-get_CI_time_origin_data(sim, sample_size, Delta_list,Gamma_list,pathDataFramesTOriginCI,
                                                    pathDataFramesTOrigin)

nB=99
pathDataFramesHybrid<-paste(path_to_save,"/","listDataFramesHybrid_",sample.size,"_",nB,"_", sim,".rds", sep="") 
all.data.frame_100_99<-assemble_output(sim, sample_size, Delta_list,Gamma_list,
                                       path_scenario_A= pathDataFramesA,
                                       path_scenario_B=pathDataFramesB,
                                       path_scenario_hybrid=pathDataFramesHybrid,
                                       path_scenario_C=pathDataFramesC, 
                                       pathDataFramesTOrigin,
                                       include_scenario_C=FALSE,  include_scenario_hybrid=FALSE, shifted_scale=FALSE, nB )

ratio.data.frame_100_99<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                        pathDataFramesA, 
                                        pathDataFramesB)

ratio.data.frame_100_A_Hybrid_99<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                                 pathDataFramesA, 
                                                 pathDataFramesHybrid)

meanTOriginDataFrame_100_99<-get_mean_time_origin_data(sim, sample_size, Delta_list,Gamma_list,
                                                       pathDataFramesTOrigin)

CI_TOriginDataFrame_100_99<-get_CI_time_origin_data(sim, sample_size, Delta_list,Gamma_list,pathDataFramesTOriginCI,
                                                    pathDataFramesTOrigin)


sample_size=10
sample.size=10
nB=5
pathDataFramesA<-paste(path_to_save,"/", "listDataFramesA_",sample.size,"_", sim,".rds", sep="")
pathDataFramesB<-paste(path_to_save,"/","listDataFramesB_",sample.size,"_", sim,".rds", sep="")
pathDataFramesC<-paste(path_to_save,"/","listDataFramesC_",sample.size,"_", sim,".rds", sep="")
pathDataFramesTOrigin<-paste(path_to_save,"/", "meansToriginA_",sample_size,"_", sim,".rds", sep="")
pathDataFramesTOriginCI<-paste(path_to_save,"/", "listToriginCI_A_",sample_size,"_", sim,".rds", sep="")

pathDataFramesB_K<-paste(path_to_save,"/","listDataFramesB_",sample.size,"_K=",K, "_", sim,".rds", sep="")
pathDataFramesHybrid<-paste(path_to_save,"/","listDataFramesHybrid_",sample.size,"_",nB,"_", sim,".rds", sep="")  

all.data.frame_10<-assemble_output(sim, sample_size, Delta_list,Gamma_list,
                                    path_scenario_A= pathDataFramesA,
                                    path_scenario_B=pathDataFramesB,
                                    path_scenario_hybrid=pathDataFramesHybrid,
                                   path_scenario_C=pathDataFramesC,  
                                   pathDataFramesTOrigin,
                                    include_scenario_C=FALSE, include_scenario_hybrid=TRUE, shifted_scale=FALSE, nB)

ratio.data.frame_10<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                    pathDataFramesA, 
                                    pathDataFramesB)

ratio.data.frame_10_A_Hybrid<-get_ratio_data(sim, sample_size, Delta_list,Gamma_list,
                                             pathDataFramesA, 
                                             pathDataFramesHybrid)


meanTOriginDataFrame_10<-get_mean_time_origin_data(sim, sample_size, Delta_list,Gamma_list,
                                                   pathDataFramesTOrigin)

CI_TOriginDataFrame_10<-get_CI_time_origin_data(sim, sample_size, Delta_list,Gamma_list,pathDataFramesTOriginCI,
                                                pathDataFramesTOrigin)

result_mean_ratio_waiting_times_respect_TMRCA<-mean_ratio_waiting_times_A_respect_2th_coal_time(number_sim, 10, Delta_list,Gamma_list,data_frame_A =all.data.frame_10[all.data.frame_10$Scenario=="A",] )

result_mean_ratio_waiting_times_respect_TOrigin<-mean_ratio_waiting_times_A_respect_time_origin(number_sim, 10, Delta_list,Gamma_list,data_frame_A =all.data.frame_10[all.data.frame_10$Scenario=="A",] , meanTOriginDataFrame_10)

sample_size=10
all.data.frame<-all.data.frame_10
meanTOriginDataFrame<-meanTOriginDataFrame_10
CI_TOriginDataFrame<-CI_TOriginDataFrame_10
ratio.data.frameAB<-ratio.data.frame_10

sample_size=100
all.data.frame<-all.data.frame_100_50
meanTOriginDataFrame<-meanTOriginDataFrame_100_50
CI_TOriginDataFrame<-CI_TOriginDataFrame_100_50
ratio.data.frameAB<-ratio.data.frame_100_50

lab1 <- c(expression(Log~median~T[1]), expression(Log~CI95~'%'~T[1]))
lab2 <- c(expression(Log~mean~T[1]), expression(Log~CI95~'%'~T[1]))

ratio.data.frame_10$Sample<-as.factor(10)
ratio.data.frame_10$x<-10*ratio.data.frame_10$x
ratio.data.frame_100_99$Sample<-as.factor(100)
ratio.data.frame.10.and.100 <- rbind(ratio.data.frame_10, ratio.data.frame_100_99)

all.data.frameA_Hybrid<-all.data.frame[substr(all.data.frame$Scenario, start = 1, stop = 1) %in%  c("A" ,"H") ]
all.data.frameB_Hybrid<-all.data.frame[substr(all.data.frame$Scenario, start = 1, stop = 1) %in%  c("B" ,"H") ]
all.data.frameB_C<-all.data.frame[all.data.frame$Scenario %in% c("B", "C"),]
############################## plots
require(ggplot2)

if (doTestA1){
  
  listDataFramesA_test<-readRDS(paste(path_to_save,"/", "listDataFramesA_test.rds", sep=""))
  
  expected.test.A1<-unlist(lapply((sample_size-1):1, function(k, n) 2*(1.0/k-1.0/n), n=sample_size))
  dataFrameA_test<-do.call("rbind", listDataFramesA1_test)
  dataFrameA_test$scenario="A"
  dataFrameA_test$Delta=rep(Delta_list, rep(sample_size-1, length(Delta_list)))
  dataFrameA_test$Gamma=rep(Gamma_list, rep(sample_size-1, length(Gamma_list)))
  dataFrameA_test$Simulation=rep(1:length(Delta_list), rep(sample_size-1, length(Delta_list)))
  dataFrameA_test$trueValues<- rep(expected.test.A1, length(Delta_list))
  dataFrameA_test$x=rep((sample_size-1):1, length(Delta_list))
  
  p1<-ggplot(dataFrameA_test, aes(x = x)) +
    geom_point(aes(y=log(mean)), colour="red")+
    geom_point(aes(y=log(trueValues)), colour="green")+
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    #geom_errorbar(aes(ymax = log(UI), ymin = log(LI)), colour="red")+
    ggtitle(paste("Comparison of simulated event times(red) vs expected(green)"))+
    labs(y="Log of event times", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    facet_grid(~ factor(Simulation))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  ggsave(paste(path_to_save, "/",  "plotTestA.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "plotTestA.png", sep=""))
  ggsave(paste(path_to_save, "/", "plotTestA.jpg", sep=""))
}




#####################################################################
require(ggplot2)

if (varyDelta){
  
  p1<-ggplot(all.data.frame, aes(x = x, y = log(mean), colour = Scenario)) +
    geom_point(size=1) +
    geom_hline(data = data.frame(Delta=factor(DeltaList), MeanTOrigin= meansTorigin, Gamma=factor(GammaList)),
               aes(yintercept = log(MeanTOrigin)))+
    geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
    ggtitle(paste("Means with 95%CI of event times"))+
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    facet_grid(~ factor(Delta))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  ggsave(paste(path_to_save, "/",  "TestPlotA.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "TestPlotA.png", sep=""))
  ggsave(paste(path_to_save, "/", "TestPlotA.jpg", sep=""))
  
  
}else{
  
  
  tiff(paste(path_to_save, "/",  "MeansCIn=",sample_size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  
  p111<-ggplot(all.data.frame, aes(x = x, y = log(mean), colour=Scenario))+
    geom_point(size=1) +
    geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    facet_grid(~ paste0("\u0393=",GammaVar))+
    ylim(min(log(all.data.frame$median)) -5,max(log(all.data.frame$median)) +1.5)+
    # theme(axis.text.x=element_text(hjust=1,vjust=.9,size=12),
    #       axis.text.y=element_text(hjust=1,vjust=.8,size=12),
    #       axis.title.x=element_text(size=12,face="bold",vjust=-0.5,hjust=0.5),
    #       axis.title.y=element_text(size=12,face="bold",vjust=1.5,hjust=0.5))+
    theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p111+ 
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(LI), linetype ="Log CI T1"),colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(UI), linetype ="Log CI T1"), colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(mean),  linetype ="Log mean T1"), colour="black")+
    scale_linetype_manual(name = "Lines", values = c(2, 1), guide = guide_legend(override.aes = list(color = c("black", "gray"))), labels = lab1)+
    theme( legend.position="none")+
             #c(0.9, 0.31))+
    theme(legend.text=element_text(size=6), legend.title=element_text(size=8))+
    theme(legend.key.size = unit(0.1, "cm"))+
    scale_fill_brewer(palette="Accent")
  
  ggsave(paste(path_to_save, "/",  "MeansCIn=",sample_size,".pdf", sep=""))
  ggsave(paste(path_to_save, "/", "MeansCIn=",sample_size,".png", sep=""))
  ggsave(paste(path_to_save, "/", "MeansCIn=",sample_size,".jpg", sep=""))
  
  dev.off()
  
  ann_text <- data.frame(x = 6,median = exp(-7),lab = paste0("\u0394=",1),
                         GammaVar =  factor(paste0("\u0393=",10), levels =unlist(lapply(Gamma_list, function(x) paste0("\u0393=",x))))
  )
  
  ann_text <- data.frame(x = 5,median = exp(-7),lab = paste0("\u0394=",1),
                         GammaVar =  factor(10, levels = Gamma_list)
  )
  
  tiff(paste(path_to_save, "/",  "MediansCIn=",sample_size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  
  p222<-ggplot(all.data.frame)+
    geom_point(size=0.5,  aes(x = x, y = log(median), colour=Scenario)) +
    geom_errorbar(data=all.data.frame, aes(x=x, y=log(median),ymax = log(UI), ymin = log(LI),  colour=Scenario))+
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/4)))+
    facet_grid(~ paste0("\u0393=",GammaVar))+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    ylim(min(log(all.data.frame$median)) -5,max(log(all.data.frame$median)) +1.5)+
    theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
     
  p222+ 
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(LI), linetype ="Log CI T1"),colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(UI), linetype ="Log CI T1"), colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(median),  linetype ="Log median T1"), colour="black")+
    scale_linetype_manual(name = "Lines", values = c(2, 1), guide = guide_legend(override.aes = list(color = c("black", "gray"))), labels = lab1)+
    theme( legend.position="none")+
    theme(legend.text=element_text(size=6), legend.title=element_text(size=8))+
    theme(legend.key.size = unit(0.1, "cm"))+
    scale_fill_brewer(palette="Accent")
    #geom_text(data = ann_text,mapping = aes(x = x, y = log(median), label = lab), size=4)
 
  
  ggsave(paste(path_to_save, "/",  "MediansCIn=",sample_size,".pdf", sep=""))
  ggsave(paste(path_to_save, "/", "MediansCIn=",sample_size,".png", sep=""))
  ggsave(paste(path_to_save, "/", "MediansCIn=",sample_size,".jpg", sep=""))

  dev.off()
  #facet_grid(~ paste0(GammaVar, "= \u0393"))+
  
}



if (varyDelta){
  
  p2<-ggplot(all.data.frame, aes(x = x, y = log(mean), colour =Scenario)) +
    geom_point(size=1) +
    geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    ylim(min(log(all.data.frame$mean)) -5,max(log(all.data.frame$mean)) +5)+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    facet_grid(~ factor(Delta))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p2+ 
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(LI)), colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(UI)), colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(mean),  linetype = "Log mean T1"), colour="black")+
    scale_linetype_manual(name = "Lines", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("black"))))
  
  
  ggsave(paste(path_to_save, "/",  "plot2TOriginGamma=1.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "plot2TOriginGamma=1.png", sep=""))
  ggsave(paste(path_to_save, "/", "plot2TOriginGamma=1.jpg", sep=""))
  
  p3<-ggplot(all.data.frame, aes(x = x, y = log(median), colour = Scenario)) +
    geom_point(size=1) +
    geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    ylim(min(log(all.data.frame$median)) -5,max(log(all.data.frame$median)) +5)+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    facet_grid(~ factor(Delta))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p3+ 
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(LI)), colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(UI)), colour="gray")+
    geom_hline(data = CI_TOriginDataFrame, aes(yintercept = log(median),  linetype = "Log median T1"), colour="black")+
    scale_linetype_manual(name = "Lines", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("black"))))
  
  
  ggsave(paste(path_to_save, "/",  "plot3TOriginUpdated.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "plot3TOriginUpdated.png", sep=""))
  ggsave(paste(path_to_save, "/", "plot3TOriginUpdated.jpg", sep=""))
}else{
  
  p2<-ggplot(all.data.frame, aes(x = x, y = log(mean), colour = Scenario)) +
    geom_point(size=1) +
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    facet_grid(~ paste0("\u0393=",GammaVar))+
    ylim(min(log(all.data.frame$mean)) -5,max(log(all.data.frame$mean)) +5)+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p2+  geom_hline(data = meanTOriginDataFrame, aes(yintercept = log(mean_time_origin),  linetype = "Log mean TOrigin"), colour="black")+
    scale_linetype_manual(name = "Lines", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("black"))))
  
  
  ggsave(paste(path_to_save, "/",  "MeansDelta=1n=100_2.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "MeansDelta=1n=100_2.png", sep=""))
  ggsave(paste(path_to_save, "/", "MeansDelta=1n=100_2.jpg", sep=""))
  
  p3<-ggplot(all.data.frame, aes(x = x, y = log(median), colour = Scenario)) +
    geom_point(size=1)+
    labs(y="Log of coalescent times", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    facet_grid(~ paste0("\u0393=",GammaVar))+
    ylim(min(log(all.data.frame$median)) -5,max(log(all.data.frame$median)) +5)+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  
  p3+  geom_hline(data = meanTOriginDataFrame, aes(yintercept = log(mean_time_origin),  linetype = "Log mean TOrigin"), colour="black")+
    scale_linetype_manual(name = "Lines", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("black"))))
  
  ggsave(paste(path_to_save, "/",  "Medians100.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "Median100.png", sep=""))
  ggsave(paste(path_to_save, "/", "Median100.jpg", sep=""))
  
  
  tiff(paste(path_to_save, "/",  "RatioMeansScenariosn=",sample_size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  
  p4<-ggplot(ratio.data.frame_10_A_Hybrid, aes(x = x, y = RatioMeanA1overA2)) +
    geom_point(size=1) +
    geom_hline(data = NULL, aes(yintercept = 1))+
    labs(y="", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    facet_grid(~ paste0("\u0393=",GammaVar))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  ggsave(paste(path_to_save, "/",  "RatioMeansDelta=1n=10.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMeansDelta=1n=10.png", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMeansDelta=1n=10.jpg", sep=""))
  
  dev.off()
  
  tiff(paste(path_to_save, "/",  "RatioMediansScenariosAHybrid50n=",sample_size,"_",50,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  
  p5<-ggplot(ratio.data.frame_100_A_Hybrid_50, aes(x = x, y = RatioMedianA1overA2)) +
    geom_point(size=0.5) +
    geom_hline(data = NULL, aes(yintercept = 1))+
    ylim(0.93,1.05)+
    labs(y="", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    facet_grid(~ paste0("\u0393=",GammaVar))+
    theme(text = element_text(size=15),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                   axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank())
  
  p5
  ggsave(paste(path_to_save, "/",  "RatioMediansA_Hybrid50n=100_2.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansA_Hybrid50n=100_2.png", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansA_Hybrid50n=100_2.jpg", sep=""))
  
  dev.off()
  
  
  tiff(paste(path_to_save, "/",  "RatioMediansA_Bn=",sample_size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  
  p5<-ggplot(ratio.data.frameAB, aes(x = x, y = RatioMedianA1overA2)) +
    geom_point(size=0.5) +
    geom_hline(data = NULL, aes(yintercept = 1))+
    scale_x_continuous(breaks=seq(0,sample_size, by=floor(sample_size/10)))+
    ylim(0.95,1.20)+
    labs(y="", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    facet_grid(~ paste0("\u0393=",GammaVar))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p5
  ggsave(paste(path_to_save, "/",  "RatioMediansA_Bn=10.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansA_Bn=10.png", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansA_Bn=10.jpg", sep=""))
  
  dev.off()
  
  ann_text <- data.frame(x = 60,RatioMedianA1overA2 = 0.95,lab = paste0("\u0394=",1),
                         GammaVar =  factor(10, levels = Gamma_list), Sample=as.factor(100))
  
  tiff(paste(path_to_save, "/",  "RatioMediansn=10Vsn=",sample_size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  
  p5<-ggplot(ratio.data.frame.10.and.100, aes(x = x, y = RatioMedianA1overA2, colour=Sample)) +
    geom_point(size=0.5) +
    geom_hline(data = NULL, aes(yintercept = 1))+
    scale_x_continuous(breaks=NULL)+
    ylim(0.95,1.15)+
    labs(y="", x ="Number of sample ancestors ")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    facet_grid(~ paste0("\u0393=",GammaVar))+
    theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme( legend.position="none")+
    theme(legend.text=element_text(size=6), legend.title=element_text(size=8))+
    theme(legend.key.size = unit(0.1, "cm"))+
    scale_fill_brewer(palette="Accent")
    #geom_text(data = ann_text,mapping = aes(x = x, y = RatioMedianA1overA2 , label = lab), size=4, colour="black")
  
  p5
  
  ggsave(paste(path_to_save, "/",  "RatioMediansn=10Vsn=100.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansn=10Vsn=100.png", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansn=10Vsn=100.jpg", sep=""))
  
  dev.off()
  
}


p4<- ggplot(comparisonTOriginDataFrame, aes(logDelta))+
  geom_point(aes(y=TOriginTransGumbel), colour="red")+
  geom_point(aes(y=TOrigin), colour="green")+
  ggtitle(paste("Comparison of average times of origin vs Delta"))+
  labs(y="Time of origin", x ="log of Delta")+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
        axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
        plot.title = element_text(size=15, face="bold", vjust=2),
        axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
        axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
        legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(path_to_save, "/",  "plot4.pdf", sep=""))
ggsave(paste(path_to_save, "/", "plot4.png", sep=""))
ggsave(paste(path_to_save, "/", "plot4.jpg", sep=""))
# ####################################################################################
#Distribution of TMRCA when n=2
TMRCA_A1_n_2<-readRDS(paste(path_to_save,"/", "simul_coal_times_A1_",1,"_",2,".rds", sep=""))
TMRCA_A2_n_2<-readRDS(paste(path_to_save,"/", "simul_coal_times_A2_",1,"_",2,".rds", sep=""))

TMRCA_A1_n_2_Gamma_1<-TMRCA_A1_n_2[[1]]
TMRCA_A2_n_2_Gamma_1<-TMRCA_A2_n_2[[1]]

TMRCA_n_2<- data.frame( TMRCA=c(TMRCA_A1_n_2_Gamma_1$coal_times_1, TMRCA_A2_n_2_Gamma_1$coal_times_1), 
              Scenario=c(rep("A", nrow(TMRCA_A1_n_2_Gamma_1)), rep("B", nrow(TMRCA_A1_n_2_Gamma_1)))
              )


library(RColorBrewer)

tiff(paste(path_to_save, "/",  "TMRCAn=",".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)


p8 <- ggplot(TMRCA_n_2, aes(x = TMRCA, fill = Scenario)) +
  geom_density(position="identity", alpha=0.6) +
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "Density of tmrca") +
  ggtitle("") +
  xlim(0,2.3)+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="none")+
  theme(legend.text=element_text(size=6), legend.title=element_text(size=8))+
  theme(legend.key.size = unit(0.1, "cm"))+
  scale_fill_brewer(palette="Accent")
p8

ggsave(paste(path_to_save, "/",  "tmrca_n_2_2.pdf", sep=""))
ggsave(paste(path_to_save, "/", "tmrca_n_2_2.png", sep=""))
ggsave(paste(path_to_save, "/", "tmrca_n_2_2.jpg", sep=""))

dev.off()