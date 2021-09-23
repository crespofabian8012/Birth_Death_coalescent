#code for testing M_k  as a race process

library(tibble)
library(tidyr)
library(dplyr)
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
print(path)
source(paste0(path,"/CoalSimulationBirthDeath.R"))
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
sim = 1000
sample.size = 10
runScenarioA = TRUE
runScenarioB = TRUE

K = 0.8

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



path_to_save = getCurrentFileLocation()

GammaList = c(10, 100)
DeltaList = GammaList

Time.Origin.STD <-
  bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
coal.events.times.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))

number.ancestors.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                   1))

coal.events.times.simB = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))




if (runScenarioA) {
  number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), sim)
  
  listDataFramesA <-
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
    listDataFramesA,
    file = paste(
      path_to_save,
      "/",
      "listDataFramesA_",
      sample.size,
      "_",
      sim,
      ".rds",
      sep = ""
    )
  )
  
  number.ancestors.Transition.Matrix <-
    as.matrix(number.ancestors.Transition[])
  saveRDS(
    number.ancestors.Transition.Matrix,
    file = paste(
      path_to_save,
      "/",
      "number.ancestors.Transition_",
      sample.size,
      "_",
      sim,
      ".rds",
      sep = ""
    )
  )
  
  
  
  listDataFramesAncestors <-
    compute_stats_number_ancestors(sample.size, sim, DeltaList, number.ancestors.simA)
  saveRDS(
    listDataFramesAncestors,
    file = paste(
      path_to_save,
      "/",
      "listDataFramesA_NumberAncestors",
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
  
  if (length(DeltaList)==1){
    meansTorigin <-mean(Time.Origin.STD[1,])
  }
  else{
    meansTorigin <-
      bigstatsr::big_apply(
        Time.Origin.STD,
        a.FUN = function(Time.Origin.STD, ind)
          rowMeans(Time.Origin.STD[ind,]),
        ind = rows_along(Time.Origin.STD),
        a.combine = 'c',
        block.size = 500
      )
    
  }
  
  
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
  
  if (length(DeltaList)==1){
    quants <- c(0.025, 0.50, 0.975)
    listToriginCI<-quantile(Time.Origin.STD[1,],  probs = quants)
  }
  else{
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
    
  }
  
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
  
  pos = sample.size - 1
  simul_coal_times_A <-
    get_simulated_kth_coal_times(coal.events.times.simA, pos, sim, sample.size)
  
  saveRDS(
    simul_coal_times_A,
    file = paste(
      path_to_save,
      "/",
      "simul_coal_times_A_",
      pos,
      "_",
      sample.size,
      "_",
      sim,
      ".rds",
      sep = ""
    )
  )
  
  list_corr_mat <-
    computeCorrelationMatrixModel(coal.events.times.simA, DeltaList, sample.size, sim)
  
  # pryr::mem_change(rm(coal.events.times.simA))
  # pryr::mem_change(rm(number.ancestors.simA))
  # pryr::mem_change(rm(number.ancestors.Transition))
  
  saveRDS(
    list_corr_mat,
    file = paste(
      path_to_save,
      "/",
      "list_corr_mat_A",
      "_",
      sample.size,
      "_",
      sim,
      ".rds",
      sep = ""
    )
  )
}


######################################################################################
#Scenario M*, M_k(cancer model)
if (runScenarioB) {
  coal.events.times.simB = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                      1))
  
  coal.events.times.simB.race = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                      1))
   K.list<-c(0.8)
  
  print("starting simulation M_k")
  
  for (i in  1:length(K.list)) {
    K = K.list[i]
    listDataFramesB <-
      simulateB_K.parallel(
        DeltaList,
        GammaList,
        sim,
        sample.size,
        Time.Origin.STD,
        coal.events.times.simB,
        K
      )
    saveRDS(
      listDataFramesB,
      file = paste(
        path_to_save,
        "/",
        "listDataFramesB_",
        sample.size,
        "_K=",
        K,
        "_",
        sim,
        ".rds",
        sep = ""
      )
    )
    
    
    list_corr_mat <-
      computeCorrelationMatrixModel(coal.events.times.simB, DeltaList, sample.size, sim)
    saveRDS(
      list_corr_mat,
      file = paste(
        path_to_save,
        "/",
        "list_corr_mat_B_",
        "_",
        sample.size,
        "_K=",
        K,
        "_",
        sim,
        ".rds",
        sep = ""
      )
    )
    
    
  }
  


  pryr::mem_change(rm(coal.events.times.simB))
 
  
  print("starting simulation M_k as a race process")
  
  for (i in  1:length(K.list)) {
    K = K.list[i]
    listDataFramesBrace <-
      simulateB_K_Race.parallel(
        DeltaList,
        GammaList,
        sim,
        sample.size,
        Time.Origin.STD,
        coal.events.times.simB.race,
        K
      )
    saveRDS(
      listDataFramesBrace,
      file = paste(
        path_to_save,
        "/",
        "listDataFramesB_race_",
        sample.size,
        "_K=",
        K,
        "_",
        sim,
        ".rds",
        sep = ""
      )
    )
    
  
  }
  
  pryr::mem_change(rm(coal.events.times.simB.race))
}
