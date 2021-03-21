#R code accompanying the paper "Coalescent Models Derived from
#Birth-Death Processes" by Crespo and Wiuf
#post_processing
rm(list=ls())
################################################################################################################
#' get_ratio_data
#' Reads and assemble the output of the simulator for scenarios A and B
#' @param number_sim: number of simulations
#' @param sample_size: sample size
#' @param Delta_list: vector of Delta  values
#' @param Gamma_list: vector of Gamma(must be of the same length of Delta_list)
#' @param path_scenario_A: path to a file containing the list of data frames generated for the simulator for scenario A 
#' @param path_scenario_B: path to a file containing the list of data frames generated for the simulator for scenario B 
#' @output the list of ratio dataframe(one per Gamma value)
################################################################################################################
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
  ratio_data_frame<-data.frame(RatioMeanA1overA2=mean_ratio_colum,RatioMedianA1overA2 =median_ratio_colum, 
                               Delta=Delta_column,
                               GammaVar=Gamma_column, x=rep((sample_size-1):1, length(list_data_fames_A)))
  
  #ratio_data_frame$Gamma <- factor(ratio_data_frame$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
  #ratio_data_frame$Delta <- factor(ratio_data_frame$Delta, levels = factor(Delta_list), labels=gamma_lbs)
  return(ratio_data_frame)
  
}
################################################################################################################
#' getCurrentFileLocation
#' get the path of the current R script
#' @output the path of the current R script
################################################################################################################
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
################################################################################################################
#' saveRatioPlot
#' save the ratio of coalescent times between model M_K and BD model
#' @param ratio.data.frameAB
#' @param path_to_save
#' @param sample_size
#' @output None
################################################################################################################
saveRatioPlot<-function(ratio.data.frameAB, path_to_save, sample_size){
  
  tiff(paste(path_to_save, "/",  "RatioMediansA_Bn=",sample_size,"_K=",K,".tiff", sep=""), units="in", width=3.25, 
       height=3.25, res=300)
  
  p5<-ggplot(ratio.data.frameAB, aes(x = x, y = RatioMedianA1overA2)) +
    geom_point(size=0.5) +
    geom_hline(data = NULL, aes(yintercept = 1))+
    scale_x_continuous(breaks=seq(0,sample_size,  length.out
                                  = 6))+
    ylim(0.8,1.20)+
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
  ggsave(paste(path_to_save, "/",  "RatioMediansA_Bn=",sample_size, "_K=",K,".pdf", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansA_Bn=",sample_size,"_K=",K, ".png", sep=""))
  ggsave(paste(path_to_save, "/", "RatioMediansA_Bn=",sample_size,"_K=",K, ".jpg", sep=""))
  
  dev.off()
  
}
################################################################################################################
sim=10000
sample_size=100
Gamma_list=c(0.001,  0.1,  10)
Delta_list=rep(1,length(Gamma_list))
################################################################################################################
require(ggplot2)
require(dplyr)
require(tidyr)

path_to_save = getCurrentFileLocation()
           
################################################################################################################
# Computes the better K
################################################################################################################
K.list= seq(0,2, by=0.1)
sample.size=10
sample_size=sample.size
nB_10 <- 9
min_norm <- 100000
best_K <- -1
for(i in  2:length(K.list)){
  
  K = K.list[i]
  fileB=paste(path_to_save,"/","listDataFramesB_",sample.size,"_",nB_10,"_K=",K, ".rds", sep="")
  
  ratio.data.frameAB<-get_ratio_data(sim, 10, Delta_list,Gamma_list,
                                      paste(path_to_save,"/", "listDataFramesA_10_",nB_10,".rds", sep=""), 
                                      fileB)

  #print(unlist(ratio.data.frameAB$RatioMedianA1overA2))
  norm_2 <- sum((ratio.data.frameAB$RatioMedianA1overA2-rep(1, length(ratio.data.frameAB$RatioMedianA1overA2)))^2)
  print(paste("For K=",K," the 2-norm is ", norm_2 ))
  if (norm_2 < min_norm)
  {
    min_norm= norm_2
    best_K = K
  }
  saveRatioPlot(ratio.data.frameAB, path_to_save, sample_size)
}

print(paste("The K that gives the minimum norm is  K=",best_K, " the 2-norm is ", min_norm, sep="" ))
################################################################################################################

