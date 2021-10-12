rm(list = ls())
library(tibble)
library(tidyr)
library(dplyr)
library(pryr)
library(ape)
library(tidyverse)
library(phytools)
getCurrentFileLocation <-  function()
{
  require(tibble)
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value,
      into = c("key", "value"),
      sep = "=",
      fill = 'right'
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
path = getCurrentFileLocation()
source("/Users/faustofabiancrespofernandez/Downloads/BD/MBE/NewApprox/M_asterVsM_2/CorrelationsHybrid/Birth_Death_coalescent/src/CoalSimulationBirthDeath.R")

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
                         list_data_frames_A, list_data_frames_B){
  
  stopifnot(length(Delta_list)==length(Gamma_list))
  
  data_frame_A<-do.call("rbind", list_data_frames_A)
  data_frame_B<-do.call("rbind", list_data_frames_B)
  
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
  
  data_frame_A$x=rep((sample_size-1):1, length(list_data_frames_A))
  data_frame_B$x=rep((sample_size-1):1, length(list_data_frames_B))
  
  mean_ratio_colum<-data_frame_A$mean /  data_frame_B$mean
  median_ratio_colum<-data_frame_A$median /  data_frame_B$median
  ratio_data_frame<-data.frame(RatioMeanA1overA2=mean_ratio_colum,RatioMedianA1overA2 =median_ratio_colum, Delta=Delta_column,
                               GammaVar=Gamma_column, x=rep((sample_size-1):1, length(list_data_frames_A)))
  
  #ratio_data_frame$Gamma <- factor(ratio_data_frame$Gamma, levels = factor(Gamma_list), labels=gamma_lbs)
  #ratio_data_frame$Delta <- factor(ratio_data_frame$Delta, levels = factor(Delta_list), labels=gamma_lbs)
  return(ratio_data_frame)
  
}

list_packages <-
  c(
    "parallel",
    "doFuture",
    "bigstatsr",
    "doRNG",
    "microbenchmark",
    "matrixStats",
    "cmna",
    "NLRoot",
    "ggplot2",
    "future.apply",
    "ape",
    "tibble",
    "ggtree",
    "geiger",
    "dplyr",
    "pracma",
    "foreach",
    "doParallel",
    "stats",
    "doSNOW",
    "tcltk",
    "ggplot2",
    "tidyverse",
    "pryr",
    "lineup"
  )


install_required_packages(list_packages)
lapply(list_packages,
       library,
       character.only = TRUE)
#########################################################################
#### Simulate data from BD coalescent 
#########################################################################
path_to_save = getCurrentFileLocation()

GammaList = c(0.001,0.1,10)
DeltaList = rep(1, length(GammaList) )
#########################################################################
#BD coal for n=10
sim = 10000
sample.size = 10
K = 0.8

Time.Origin.STD <-
  bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
coal.events.times.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))
number.ancestors.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                   1))
number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), sim)

listDataFramesBDCoal10 <-
  simulateA.parallel(
    DeltaList,
    GammaList,
    sim,
    sample.size,
    Time.Origin.STD,
    coal.events.times.simA,
    number.ancestors.simA,
    number.ancestors.Transition
  )
saveRDS(
  listDataFramesBDCoal10,
  file = paste(
    path_to_save,
    "/",
    "listDataFramesBDCoal10_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

saveRDS(
  Time.Origin.STD,
  file = paste(
    path_to_save,
    "/",
    "Time.Origin.STD_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

meansTorigin <-
  bigstatsr::big_apply(
    Time.Origin.STD,
    a.FUN = function(Time.Origin.STD, ind)
      rowMeans(Time.Origin.STD[ind,]),
    ind = rows_along(Time.Origin.STD),
    a.combine = 'c',
    block.size = 500
  )

saveRDS(
  meansTorigin,
  file = paste(
    path_to_save,
    "/",
    "meansToriginA_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

listToriginCI <-
  bigstatsr::big_apply(
    Time.Origin.STD,
    a.FUN = function(Time.Origin.STD, ind) {
      quants <- c(0.025, 0.50, 0.975)
      matrixStats::rowQuantiles(Time.Origin.STD[ind,],  probs = quants)
    } ,
    ind = rows_along(Time.Origin.STD),
    a.combine = "rbind",
    block.size = 500
  )
saveRDS(
  listToriginCI,
  file = paste(
    path_to_save,
    "/",
    "listToriginCI_A_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

simulationColumn=rep(1:length(GammaList), rep(sample.size, length(GammaList)))

Delta_column=rep(DeltaList, rep(sample.size-1, length(DeltaList)))

Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))

data_frame_BDCoal10<-do.call("rbind", listDataFramesBDCoal10)
data_frame_BDCoal10$GammaVar=Gamma_column

data_frame_BDCoal10$x=10*rep(((sample.size-1):1), length(listDataFramesBDCoal10))

data_frame_BDCoal10$method="BDCoal"
###############################################################################
#########################################################################
#BD coal for n=100
sim = 10000
sample.size = 100


Time.Origin.STD <-
  bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
coal.events.times.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))
number.ancestors.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                   1))
number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), sim)

listDataFramesBDCoal100 <-
  simulateA.parallel(
    DeltaList,
    GammaList,
    sim,
    sample.size,
    Time.Origin.STD,
    coal.events.times.simA,
    number.ancestors.simA,
    number.ancestors.Transition
  )
saveRDS(
  listDataFramesBDCoal100,
  file = paste(
    path_to_save,
    "/",
    "listDataFramesBDCoal100_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

saveRDS(
  Time.Origin.STD,
  file = paste(
    path_to_save,
    "/",
    "Time.Origin.STD_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

meansTorigin <-
  bigstatsr::big_apply(
    Time.Origin.STD,
    a.FUN = function(Time.Origin.STD, ind)
      rowMeans(Time.Origin.STD[ind,]),
    ind = rows_along(Time.Origin.STD),
    a.combine = 'c',
    block.size = 500
  )

saveRDS(
  meansTorigin,
  file = paste(
    path_to_save,
    "/",
    "meansToriginA_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

listToriginCI <-
  bigstatsr::big_apply(
    Time.Origin.STD,
    a.FUN = function(Time.Origin.STD, ind) {
      quants <- c(0.025, 0.50, 0.975)
      matrixStats::rowQuantiles(Time.Origin.STD[ind,],  probs = quants)
    } ,
    ind = rows_along(Time.Origin.STD),
    a.combine = "rbind",
    block.size = 500
  )
saveRDS(
  listToriginCI,
  file = paste(
    path_to_save,
    "/",
    "listToriginCI_A_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

simulationColumn=rep(1:length(GammaList), rep(sample.size, length(GammaList)))

Delta_column=rep(DeltaList, rep(sample.size-1, length(DeltaList)))

Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))

data_frame_BDCoal100<-do.call("rbind", listDataFramesBDCoal100)
data_frame_BDCoal100$GammaVar=Gamma_column

data_frame_BDCoal100$x=rep(((sample.size-1):1), length(listDataFramesBDCoal100))

data_frame_BDCoal100$method="BDCoal"
###############################################################################
#Kingman Coal for n=10
sample.size = 10
Delta_column=rep(DeltaList, rep(sample.size-1, length(DeltaList)))
Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))

coal.events.times.simKingman = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))

listDataFramesKingman10<-simulateKingmanCoal.parallel(DeltaList, GammaList,
                                        sim,
                                        sample.size,
                                        coal.events.times.simKingman)

data_frame_Kingman10<-do.call("rbind", listDataFramesKingman10)
data_frame_Kingman10$GammaVar=Gamma_column
data_frame_Kingman10$x=10*rep(((sample.size-1):1), length(listDataFramesKingman10))
data_frame_Kingman10$method="Kingman"

#################################################################
##################################################################
#Kingman Coal for n=100
sample.size = 100
Delta_column=rep(DeltaList, rep(sample.size-1, length(DeltaList)))
Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))

coal.events.times.simKingman100 = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                          1))

listDataFramesKingman100<-simulateKingmanCoal.parallel(DeltaList, GammaList,
                                                      sim,
                                                      sample.size,
                                                      coal.events.times.simKingman100)

data_frame_Kingman100<-do.call("rbind", listDataFramesKingman100)
data_frame_Kingman100$GammaVar=Gamma_column
data_frame_Kingman100$x=rep(((sample.size-1):1), length(listDataFramesKingman100))
data_frame_Kingman100$method="Kingman"
#################################################################
#Exponential coal for n=10
sample.size = 10
Delta_column=rep(DeltaList, rep(sample.size-1, length(DeltaList)))
Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))

coal.events.times.sim.Exp10 = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))
listDataFramesExpCoal10<-simulateExponentialCoal.parallel(DeltaList, GammaList,
                                            sim,
                                            sample.size,
                                            coal.events.times.sim.Exp10)

data_frame_ExpCoal10<-do.call("rbind", listDataFramesExpCoal10)
data_frame_ExpCoal10$GammaVar=Gamma_column
data_frame_ExpCoal10$x=10*rep(((sample.size-1):1), length(listDataFramesExpCoal10))
data_frame_ExpCoal10$method="ExpCoal"
######################################################
#################################################################
#Exponential coal for n=100
sample.size = 100
Delta_column=rep(DeltaList, rep(sample.size-1, length(DeltaList)))
Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))

coal.events.times.sim.Exp100 = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                       1))
listDataFramesExpCoal100<-simulateExponentialCoal.parallel(DeltaList, GammaList,
                                                          sim,
                                                          sample.size,
                                                          coal.events.times.sim.Exp100)

data_frame_ExpCoal100<-do.call("rbind", listDataFramesExpCoal100)
data_frame_ExpCoal100$GammaVar=Gamma_column
data_frame_ExpCoal100$x=rep(((sample.size-1):1), length(listDataFramesExpCoal100))
data_frame_ExpCoal100$method="ExpCoal"
###############################################################################################
all.data.frame10<-rbind(data_frame_BDCoal10,data_frame_Kingman10 ,data_frame_ExpCoal10)
all.data.frame100<-rbind(data_frame_BDCoal100,data_frame_Kingman100 ,data_frame_ExpCoal100)
all.data.frame10$n=10
all.data.frame100$n=100
all.data.frame=rbind(all.data.frame10, all.data.frame100)
#########################################################################################

lab1 <- c(expression(Log~median~T[1]), expression(Log~CI95~'%'~T[1]))
lab2 <- c(expression(Log~mean~T[1]), expression(Log~CI95~'%'~T[1]))
library(ggplot2)

tiff(paste(path_to_save, "/",  "Comparison=",sample.size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)

p111<-ggplot(all.data.frame, aes(x = x, y = log(median), colour=method))+
  geom_point(size=1) +
  geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
  labs(y="", x ="Number of sample ancestors ")+
  scale_x_continuous(breaks=seq(1,sample.size, by=floor(sample.size/5)))+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  facet_grid(~ paste0("\u0393=",GammaVar))+
  ylim(min(log(all.data.frame$mean)) -5,max(log(all.data.frame$mean)) +1.5)+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="none")+
  theme(legend.text=element_text(size=6), legend.title=element_text(size=6), axis.title.y=element_text(size=6), axis.title.x=element_text(size=6))+
  #theme(legend.key.size = unit(0.1, "cm"))+
  scale_fill_brewer(palette="Accent")

  p111


ggsave(paste(path_to_save, "/",  "Comparison=",sample.size,".pdf", sep=""))
ggsave(paste(path_to_save, "/", "Comparison=",sample.size,".png", sep=""))
ggsave(paste(path_to_save, "/", "Comparison=",sample.size,".jpg", sep=""))

dev.off()



all.data.frame1<-all.data.frame[all.data.frame$method %in% c("BDCoal" ,"ExpCoal"),  ]
tiff(paste(path_to_save, "/",  "ComparisonBDCoalvsExp=",sample.size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)


p111<-ggplot(all.data.frame1, aes(x = x, y = log(median), colour=method))+
  geom_point(size=1) +
  geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
  labs(y="", x ="Number of sample ancestors ")+
  scale_x_continuous(breaks=seq(1,sample.size, by=floor(sample.size/5)))+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  facet_grid(~ paste0("\u0393=",GammaVar))+
  ylim(min(log(all.data.frame1$median)) -5,max(log(all.data.frame1$median)) +1.5)+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="none")+
  theme(legend.text=element_text(size=6), legend.title=element_text(size=6), axis.title.y=element_text(size=6), axis.title.x=element_text(size=6))+
  #theme(legend.key.size = unit(0.1, "cm"))+
  scale_fill_brewer(palette="Accent")

p111

ggsave(paste(path_to_save, "/",  "ComparisonBDCoalvsExp=",sample.size,".pdf", sep=""))
ggsave(paste(path_to_save, "/", "ComparisonBDCoalvsExp=",sample.size,".png", sep=""))
ggsave(paste(path_to_save, "/", "ComparisonBDCoalvsExp=",sample.size,".jpg", sep=""))

dev.off()


all.data.frame2<-all.data.frame[all.data.frame$method %in% c("Kingman" ,"ExpCoal"),  ]
tiff(paste(path_to_save, "/",  "ComparisonKingmanvsExp=",sample.size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)


p111<-ggplot(all.data.frame2, aes(x = x, y = log(median), colour=method))+
  geom_point(size=1) +
  geom_errorbar(aes(ymax = log(UI), ymin = log(LI)))+
  labs(y="", x ="Number of sample ancestors ")+
  scale_x_continuous(breaks=seq(1,sample.size, by=floor(sample.size/5)))+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  facet_grid(~ paste0("\u0393=",GammaVar))+
  ylim(min(log(all.data.frame2$median)) -5,max(log(all.data.frame2$median)) +1.5)+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="none")+
  theme(legend.text=element_text(size=6), legend.title=element_text(size=6), axis.title.y=element_text(size=6), axis.title.x=element_text(size=6))+
  #theme(legend.key.size = unit(0.1, "cm"))+
  scale_fill_brewer(palette="Accent")

p111

ggsave(paste(path_to_save, "/",  "ComparisonKingmanvsExp=",sample.size,".pdf", sep=""))
ggsave(paste(path_to_save, "/", "ComparisonKingmanvsExp=",sample.size,".png", sep=""))
ggsave(paste(path_to_save, "/", "ComparisonKingmanvsExp=",sample.size,".jpg", sep=""))

dev.off()
###################################################################3
#ratios
ratio.data.frame_BDCoal_vs_ExpCoal_n_10<-get_ratio_data(sim, 10, DeltaList,GammaList,
                                        listDataFramesBDCoal10, 
                                        listDataFramesExpCoal10)

ratio.data.frame_BDCoal_vs_Kingman_n_10<-get_ratio_data(sim, 10, DeltaList,GammaList,
                                                   listDataFramesBDCoal10, 
                                                   listDataFramesKingman10)

ratio.data.frame_BDCoal_vs_Kingman_n_10$n= 10
ratio.data.frame_BDCoal_vs_ExpCoal_n_10$n= 10

ratio.data.frame_BDCoal_vs_Kingman_n_10$x= 10*ratio.data.frame_BDCoal_vs_Kingman_n_10$x
ratio.data.frame_BDCoal_vs_ExpCoal_n_10$x= 10*ratio.data.frame_BDCoal_vs_ExpCoal_n_10$x

ratio.data.frame_BDCoal_vs_Kingman_n_10$model_ratio= "BDCoal_Kingman"
ratio.data.frame_BDCoal_vs_ExpCoal_n_10$model_ratio= "BDCoal_ExpCoal"
  
ratio.data.frame_BDCoal_vs_ExpCoal_n_100<-get_ratio_data(sim, 100, DeltaList,GammaList,
                                                        listDataFramesBDCoal100, 
                                                        listDataFramesExpCoal100)

ratio.data.frame_BDCoal_vs_Kingman_n_100<-get_ratio_data(sim, 100, DeltaList,GammaList,
                                                        listDataFramesBDCoal100, 
                                                        listDataFramesKingman100)

ratio.data.frame_BDCoal_vs_Kingman_n_100$n= 100
ratio.data.frame_BDCoal_vs_ExpCoal_n_100$n= 100

ratio.data.frame_BDCoal_vs_Kingman_n_100$model_ratio= "BDCoal_Kingman"
ratio.data.frame_BDCoal_vs_ExpCoal_n_100$model_ratio= "BDCoal_ExpCoal"

all.ratio.data.frame=rbind(ratio.data.frame_BDCoal_vs_Kingman_n_10, ratio.data.frame_BDCoal_vs_ExpCoal_n_10,
                           ratio.data.frame_BDCoal_vs_Kingman_n_100, ratio.data.frame_BDCoal_vs_ExpCoal_n_100)
library(tidyverse)
require(tidyr)
library(lemon)
all.ratio.data.frame<-unite(all.ratio.data.frame, ModelSample, model_ratio:n, sep='', remove=F)

all.ratio.data.frameBDCoal_ExpCoal<-all.ratio.data.frame[all.ratio.data.frame$model_ratio == "BDCoal_ExpCoal",]

all.ratio.data.frameBDCoal_ExpCoal$x[all.ratio.data.frameBDCoal_ExpCoal$n == 10] <- all.ratio.data.frameBDCoal_ExpCoal$x[all.ratio.data.frameBDCoal_ExpCoal$n == 10] / 10
count <- 0
labels_fun <- function(x) {
  count <<- count + 1L
  print(paste0("count=",count))
  switch(
    count,
    c("1","3","6","9"),
    c("1", "3","6","9"),
    c( "1","3","6","9"),
    c("1","30",  "60", "90"),
    c("1","30",  "60", "90"),
    c("1", "30",  "60", "90")
  )
}
breaks_fun2 <- function(limits) {
  limits = c(floor(limits[1]), ceiling(limits[2]))
  print(limits)
  d = diff(limits)
  print(d)
  if (max(d) < 95) {
    c(1,3,6,9)
  } else {
    c(1,30,60,90)
  }
}

count <- 0
all.ratio.data.frameBDCoal_ExpCoal$x[all.ratio.data.frameBDCoal_ExpCoal$Sample == 10] <- all.ratio.data.frameBDCoal_ExpCoal$x[all.ratio.data.frameBDCoal_ExpCoal$Sample == 10] / 10
tiff(paste(path_to_save, "/",  "Ratio_Medians_BDCoal_over_ExpCoal",".tiff", sep=""), units="in", width=5, height=5, res=300)

#quartz(type = 'pdf', file = paste(path_to_save, "/",  "Ratio_Medians_BDCoal_over_ExpCoal.pdf", sep=""))

p5<-ggplot(all.ratio.data.frameBDCoal_ExpCoal, aes(x = x, y = RatioMedianA1overA2, colour=ModelSample)) +
  geom_point(size=0.5) +
  geom_hline(data = NULL, aes(yintercept = 1))+
  scale_x_continuous(breaks=breaks_fun2, labels =labels_fun )+
  ylim(min(all.ratio.data.frameBDCoal_ExpCoal$RatioMedianA1overA2),max(all.ratio.data.frameBDCoal_ExpCoal$RatioMedianA1overA2))+
  labs(y="", x ="Number of sample ancestors ")+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  #facet_grid(paste0("n=",n) ~ paste0("\u0393=",GammaVar))+
  facet_wrap(~interaction(paste0("\u0393=",GammaVar), paste0(" n=",n)), scales ="free_x")+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="none")+
  theme(legend.text=element_text(size=6), legend.title=element_text(size=8))+
  theme(legend.key.size = unit(0.1, "cm"))+
  scale_fill_brewer(palette="Accent")
#geom_text(ata = ann_text,mapping = aes(x = x, y = RatioMedianA1overA2 , label = lab), size=4, colour="black")

p5
count <- 0
#ggsave(paste(path_to_save, "/",  "Ratio_Medians_BDCoal_over_ExpCoal.pdf", sep=""))
ggsave(paste(path_to_save, "/", "Ratio_Medians_BDCoal_over_ExpCoal.png", sep=""), device ="png", units="in", width=5, height=5, dpi=300 )
#ggsave(paste(path_to_save, "/", "Ratio_Medians_BDCoal_over_ExpCoal.jpg", sep=""))

dev.off()


tiff(paste(path_to_save, "/",  "Ratio_Medians_Models",".tiff", sep=""), units="in", width=5, height=5, res=300)


quartz(type = 'pdf', file = paste(path_to_save, "/",  "Ratio_Medians_Models.pdf", sep=""))

p5<-ggplot(all.ratio.data.frame, aes(x = x, y = RatioMedianA1overA2, colour=ModelSample)) +
  geom_point(size=0.5) +
  geom_hline(data = NULL, aes(yintercept = 1))+
  scale_x_continuous(breaks=c(1,  20,  40, 60, 80, 100))+
  ylim(min(all.ratio.data.frame$RatioMedianA1overA2),max(all.ratio.data.frame$RatioMedianA1overA2))+
  labs(y="", x ="Number of sample ancestors ")+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  facet_grid(paste0("n=",n) ~ paste0("\u0393=",GammaVar))+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="bottom")+
  theme(legend.text=element_text(size=6), legend.title=element_text(size=8))+
  theme(legend.key.size = unit(0.1, "cm"))+
  scale_fill_brewer(palette="Accent")
#geom_text(data = ann_text,mapping = aes(x = x, y = RatioMedianA1overA2 , label = lab), size=4, colour="black")

p5

#ggsave(paste(path_to_save, "/",  "Ratio_Medians_Models.pdf", sep=""))
ggsave(paste(path_to_save, "/", "Ratio_Medians_Models.png", sep=""))
ggsave(paste(path_to_save, "/", "Ratio_Medians_Models.jpg", sep=""))

dev.off()

#####################################################################################
ratio.data.frame_BDCoal_over_KingmanCoal_n_10<-get_ratio_data(sim, 10, DeltaList,GammaList,
                                                        listDataFramesBDCoal10, 
                                                        listDataFramesKingman10)

ratio.data.frame_ExpCoal_over_KingmanCoal_n_10<-get_ratio_data(sim, 10, DeltaList,GammaList,
                                                        listDataFramesExpCoal10, 
                                                        listDataFramesKingman10)

ratio.data.frame_BDCoal_over_KingmanCoal_n_10$n= 10
ratio.data.frame_ExpCoal_over_KingmanCoal_n_10$n= 10

ratio.data.frame_BDCoal_over_KingmanCoal_n_10$x= 10*ratio.data.frame_BDCoal_over_KingmanCoal_n_10$x
ratio.data.frame_ExpCoal_over_KingmanCoal_n_10$x= 10*ratio.data.frame_ExpCoal_over_KingmanCoal_n_10$x

ratio.data.frame_BDCoal_over_KingmanCoal_n_10$model_ratio= "BDCoal/KingmanCoal"
ratio.data.frame_ExpCoal_over_KingmanCoal_n_10$model_ratio= "ExpCoal/KingmanCoal"

ratio.data.frame_BDCoal_over_KingmanCoal_n_100<-get_ratio_data(sim, 100, DeltaList,GammaList,
                                                         listDataFramesBDCoal100, 
                                                         listDataFramesKingman100)

ratio.data.frame_ExpCoal_over_KingmanCoal_n_100<-get_ratio_data(sim, 100, DeltaList,GammaList,
                                                                listDataFramesExpCoal100, 
                                                         listDataFramesKingman100)

ratio.data.frame_BDCoal_over_KingmanCoal_n_100$n= 100
ratio.data.frame_ExpCoal_over_KingmanCoal_n_100$n= 100

ratio.data.frame_BDCoal_over_KingmanCoal_n_100$model_ratio= "BDCoal/KingmanCoal"
ratio.data.frame_ExpCoal_over_KingmanCoal_n_100$model_ratio= "ExpCoal/KingmanCoal"

all.ratio.data.frame=rbind(ratio.data.frame_BDCoal_over_KingmanCoal_n_10, ratio.data.frame_ExpCoal_over_KingmanCoal_n_10,
                           ratio.data.frame_BDCoal_over_KingmanCoal_n_100, ratio.data.frame_ExpCoal_over_KingmanCoal_n_100)
library(tidyverse)
require(tidyr)
all.ratio.data.frame<-unite(all.ratio.data.frame, ModeloverModel, model_ratio:n, sep='', remove=F)

count <- 0
labels_fun <- function(x) {
  count <<- count + 1L
  switch(
    count,
    c("1","3","6","9"),
    c("1","3","6","9"),
    c("1","3","6","9"),
    c("1", "30",  "60", "90"),
    c("1", "30",  "60", "90"),
    c("1", "30",  "60", "90")
  )
}


tiff(paste(path_to_save, "/",  "Ratio_Medians_Models_over_KingmanCoal",".tiff", sep=""), units="in", width=5, height=5, res=300)

count <- 0
#quartz(type = 'pdf', file = paste(path_to_save, "/",  "Ratio_Medians_Models_over_KingmanCoal.pdf", sep=""))

p5<-ggplot(all.ratio.data.frame, aes(x = x, y = RatioMedianA1overA2, colour=ModeloverModel)) +
  geom_point(size=0.5) +
  geom_hline(data = NULL, aes(yintercept = 1))+
  scale_x_continuous(breaks=c(1, 30, 60, 90),limits=c(1, NA), labels= labels_fun)+
  #ylim(min(all.ratio.data.frame$RatioMedianA1overA2),max(all.ratio.data.frame$RatioMedianA1overA2))+
  labs(y="", x ="Number of sample ancestors ")+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  facet_rep_wrap(~interaction(paste0("\u0393=",GammaVar), paste0("n=",n)), scales ="free")+
  #facet_grid(paste0("n=",n) ~ paste0("\u0393=",GammaVar))+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme( legend.position="bottom")+
  theme(legend.text=element_text(size=6), legend.title=element_blank())+
  theme(legend.key.size = unit(0.1, "cm"))+
  scale_colour_brewer(palette = "Set1")
  #scale_fill_brewer(palette="Set1")
#geom_text(data = ann_text,mapping = aes(x = x, y = RatioMedianA1overA2 , label = lab), size=4, colour="black")

p5


ggsave(paste(path_to_save, "/", "Ratio_Medians_Models_over_KingmanCoal.png", sep=""))
ggsave(paste(path_to_save, "/", "Ratio_Medians_Models_over_KingmanCoal.jpg", sep=""))

dev.off()
