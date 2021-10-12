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
sim = 100
sample.size = 100
runScenarioA = FALSE
runScenarioC = FALSE
runScenarioB = FALSE
doTestA = FALSE
varyDelta = FALSE
doJustBenchMark = TRUE
runHybridScenario = FALSE
number.sim.trees = 5
simulate.Trees = FALSE
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

GammaList = c(0.001,  0.1,  10)
DeltaList = rep(1, length(GammaList))

Time.Origin.STD <-
  bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
coal.events.times.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))

number.ancestors.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                   1))

Time.Origin.STD.trees <-
  bigstatsr::FBM(length(DeltaList),
                 number.sim.trees ,
                 type = "double",
                 init = 0)
coal.events.times.simA.trees = bigstatsr::FBM(length(DeltaList), number.sim.trees *
                                                (sample.size - 1))

coal.events.times.simB = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))

coal.events.times.simHybrid = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                         1))

nB = floor(sample.size / 2)


nB_List <- c(floor(sample.size / 2), floor(sample.size / 4))


if (doTestA)
{
  Time.Origin.STD.test <-
    bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
  coal.events.times.sim.testA = bigstatsr::FBM(length(DeltaList), sim *
                                                 (sample.size - 1))
  
  listDataFramesA_test <-
    simulateA.test.parallel(
      DeltaList,
      GammaList,
      sim,
      sample.size,
      Time.Origin.STD.test,
      coal.events.times.sim.testA
    )
  saveRDS(
    listDataFramesA_test,
    paste(
      path_to_save,
      "/",
      "listDataFramesA_",
      sample.size,
      "_test.rds",
      sep = ""
    )
  )
  if (length(DeltaList)==1){
    meansTorigin.test <-mean(Time.Origin.STD.test[1,])
  }
  else{
    meansTorigin.test <-
      bigstatsr::big_apply(
        Time.Origin.STD.test,
        a.FUN = function(Time.Origin.STD.test, ind)
          rowMeans(Time.Origin.STD.test[ind,]),
        ind = rows_along(Time.Origin.STD.test),
        a.combine = 'c',
        block.size = 100
      )
    
  }
  
  saveRDS(
    meansTorigin.test,
    file = paste(
      path_to_save,
      "/",
      "meansToriginA_",
      sample.size,
      "_test.rds",
      sep = ""
    )
  )
  
  
  pryr::mem_change(rm(Time.Origin.STD.test))
  pryr::mem_change(rm(coal.events.times.sim.testA))
  
}
if (!doJustBenchMark) {
  ######################################################################################
  # Scenario BD coalescent with stochastic population size
  
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
  else{
    print(paste(
      path_to_save,
      "/",
      "listDataFramesA_",
      sample.size,
      ".rds",
      sep = ""
    ))
    listDataFramesA <-
      readRDS(
        paste(
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
    
    print(
      paste(
        path_to_save,
        "/",
        "Time.Origin.STD_",
        sample.size,
        "_",
        nB,
        ".rds",
        sep = ""
      )
    )
    Time.Origin.STD <-
      readRDS(
        paste(
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
    
    print(paste(
      path_to_save,
      "/",
      "meansToriginA_",
      sample.size,
      ".rds",
      sep = ""
    ))
    meansTorigin <-
      readRDS(paste(
        path_to_save,
        "/",
        "meansToriginA_",
        sample.size,
        "_",
        sim,
        ".rds",
        sep = ""
      ))
  }
  
  
  #######################################################################################
  #hybrid scenario
  
  if (runHybridScenario) {
    number.Fail.hybrid = bigstatsr::FBM(length(nB_List), length(GammaList))
    coal.events.times.simHybrid = bigstatsr::FBM(length(DeltaList), sim *
                                                   (sample.size - 1))
    number.ancestors.simHybrid = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                            1))
    number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), sim)
    
    K = 2
    
    listDataFramesHybrid <-
      simulateHybrid.parallel(
        DeltaList,
        GammaList,
        sim,
        sample.size,
        Time.Origin.STD,
        coal.events.times.simHybrid,
        nB,
        number.Fail.hybrid,
        1,
        K,
        number.ancestors.simHybrid,
        number.ancestors.Transition
      )
    
    saveRDS(
      listDataFramesHybrid,
      file = paste(
        path_to_save,
        "/",
        "listDataFramesHybrid_",
        "K=",
        K,
        "_",
        sample.size,
        "_",
        nB,
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
        "K=",
        K,
        "_",
        sample.size,
        "_",
        nB,
        "_",
        sim,
        ".rds",
        sep = ""
      )
    )
    
    
    print(number.Fail.hybrid[])
    
    listDataFramesAncestors <-
      compute_stats_number_ancestors(sample.size, sim, DeltaList, number.ancestors.simHybrid)
    saveRDS(
      listDataFramesAncestors,
      file = paste(
        path_to_save,
        "/",
        "listDataFramesHybrid_NumberAncestors",
        sample.size,
        "_",
        sim,
        ".rds",
        sep = ""
      )
    )
    
    list_corr_mat <-
      computeCorrelationMatrixModel(coal.events.times.simHybrid, DeltaList, sample.size, sim)
    
    saveRDS(
      list_corr_mat,
      file = paste(
        path_to_save,
        "/",
        "list_corr_mat_Hybrid",
        "_",
        "K=",
        K,
        "_",
        sample.size,
        "_",
        nB,
        "_",
        sim,
        ".rds",
        sep = ""
      )
    )
    
    pryr::mem_change(rm(coal.events.times.simHybrid))
    pryr::mem_change(rm(number.ancestors.simHybrid))
    pryr::mem_change(rm(number.ancestors.Transition))
  }
  
  ######################################################################################
  #Scenario M*, M_k(cancer model)
  if (runScenarioB) {
    coal.events.times.simB = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                        1))
   # print("starting simulation M*")
    # #K.list= seq(0.5,2, by=0.1)
    # #K.list= seq(2,2, by=1)
     K.list = seq(0.5,1.2, by=0.1)
    # K.list<-c(0.8)
    
    # listDataFramesB <-
    #   simulateB_MAster.parallel(DeltaList,
    #                             GammaList,
    #                             sim,
    #                             sample.size,
    #                             Time.Origin.STD,
    #                             coal.events.times.simB)
    # saveRDS(
    #   listDataFramesB,
    #   file = paste(
    #     path_to_save,
    #     "/",
    #     "listDataFramesB_",
    #     sample.size,
    #     "_",
    #     sim,
    #     ".rds",
    #     sep = ""
    #   )
    # )
  
   
    
    
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
      
      # if (i < length(K.list)) {

      #   listDataFramesA <-
      # simulateA.parallel(
      #   DeltaList,
      #   GammaList,
      #   sim,
      #   sample.size,
      #   Time.Origin.STD,
      #   coal.events.times.simA,
      #   number.ancestors.simA,
      #   number.ancestors.Transition
      # )
    
      # }
      
    }
    
    pos = sample.size - 1
    simul_coal_times_B <-
      get_simulated_kth_coal_times(coal.events.times.simB, pos, sim, sample.size)
    
    saveRDS(
      simul_coal_times_B,
      file = paste(
        path_to_save,
        "/",
        "simul_coal_times_B_",
        sample.size,
        ".rds",
        sep = ""
      )
    )
    pryr::mem_change(rm(coal.events.times.simB))
    
    
  }
  
  ######################################################################################
  #Scenario C
  
  if (runScenarioC)
  {
    coal.events.times.simC = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                        1))
    
    listDataFramesC <-
      simulateC.parallel(DeltaList,
                         GammaList,
                         sim,
                         sample.size,
                         coal.events.times.simC)
    saveRDS(
      listDataFramesC,
      file = paste(
        path_to_save,
        "/listDataFramesC_",
        sample.size,
        "_",
        nB,
        ".rds",
        sep = ""
      )
    )
    
    pryr::mem_change(rm(coal.events.times.simC))
    
  }
  
  ##############################################################################################
  #Draw some trees: number.sim.trees trees per Delta and Gamma combination
  
  if (simulate.Trees) {
    packageVersion("phangorn")
    packageVersion("phytools")
    
    print("Simulating some sample trees...")
    
    listDataFramesA.trees <-
      simulateA.parallel(
        DeltaList,
        GammaList,
        number.sim.trees,
        sample.size,
        Time.Origin.STD.trees,
        coal.events.times.simA.trees
      )
    
    lapply(1:length(GammaList), function(param_index,
                                         coal.events.times.simA.trees,
                                         number.sim.trees,
                                         sample.size,
                                         path_to_save) {
      random_simulated_trees <-
        get_list_simulated_trees(coal.events.times.simA.trees,
                                 param_index,
                                 number.sim.trees,
                                 sample.size)
      
      lapply(1:number.sim.trees, function(idx,
                                          param_index,
                                          GammaList,
                                          path_to_save,
                                          random_simulated_trees) {
        tree = random_simulated_trees[[idx]]
        tree <- ape::as.phylo(tree)
        tiff(
          paste(
            path_to_save,
            "/",
            "tree_n=",
            sample.size,
            "_Gamma=",
            GammaList[param_index],
            "_",
            idx,
            ".tiff",
            sep = ""
          ),
          units = "in",
          width = 3.25,
          height = 3.25,
          res = 300
        )
        plot(tree, show.tip.label = FALSE)
        dev.off()
      }, param_index, GammaList, path_to_save, random_simulated_trees)
      
      
    }, coal.events.times.simA.trees, number.sim.trees, sample.size, path_to_save)
    
    pryr::mem_change(rm(coal.events.times.simA.trees))
    pryr::mem_change(rm(Time.Origin.STD.trees))
    
  }
  
} else{
  #################################################################################################
  #Benchmarking
  print("Starting simulation benchmark...")

   number_runs=100
   sample.size=10
   
   sim=1
   GammaList = c( 10)
   DeltaList = rep(1, length(GammaList))


  Time.Origin.STD <-
  bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
  

   coal.events.times.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))

    number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), sim)

    number.ancestors.simHybrid = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                        1))
    
    nB_List <- c(floor(sample.size / 2), floor(sample.size / 4))
  
    number.Fail.hybrid = bigstatsr::FBM(length(nB_List), length(GammaList))
    coal.events.times.simHybrid = bigstatsr::FBM(length(DeltaList), sim *
                                                   (sample.size - 1))
    number.ancestors.simHybrid = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                            1))
 

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


  print(Time.Origin.STD[])
  K = 1

 
  number.Fail.hybrid = bigstatsr::FBM(length(nB_List), length(GammaList))

  mbm <-
    microbenchmark::microbenchmark(
      "BD" = listDataFramesA <-
       simulateA ( DeltaList[1],
                  GammaList[1],
                     sim,
                     sample.size,
                     Time.Origin.STD,
                     1,
                     coal.events.times.simA,
                     number.ancestors.simA,
                     number.ancestors.Transition),
        "M*" = listDataFramesB <-
             simulateB_MAster (1,  DeltaList[1],
                                   GammaList[1],
                                    sim,
                                    sample.size,
                                    Time.Origin.STD,
                                     coal.events.times.simB),

        "M0.8" = listDataFramesB <-simulateB_K (1,
                       DeltaList[1],
                       GammaList[1],
                       sim,
                       sample.size,
                       Time.Origin.STD,
                       coal.events.times.simB,
                       0.8),

        "M2" = listDataFramesB <-simulateB_K (1,
                       DeltaList[1],
                       GammaList[1],
                       sim,
                       sample.size,
                       Time.Origin.STD,
                       coal.events.times.simB,
                       2),
        "Hn/2" =simulateHybrid(1,DeltaList[1],
                                GammaList[1],
                                   sim,
                                   sample.size,
                                   Time.Origin.STD,
                                   coal.events.times.simHybrid,
                                   nB_List[1],
                                   number.Fail.hybrid,
                                   1,
                                   2,
                                   number.ancestors.simHybrid,
                                   number.ancestors.Transition),
        "Hn/4" =simulateHybrid(1,DeltaList[1],
                                GammaList[1],
                                   sim,
                                   sample.size,
                                   Time.Origin.STD,
                                   coal.events.times.simHybrid,
                                   nB_List[2],
                                   number.Fail.hybrid,
                                   2,
                                   2,
                                   number.ancestors.simHybrid,
                                   number.ancestors.Transition)

        # simulateA.parallel(
        #   DeltaList,
        #   GammaList,
        #   sim,
        #   sample.size,
        #   Time.Origin.STD,
        #   coal.events.times.simA,
        #   number.ancestors.simA,
        #   number.ancestors.Transition
        # )
        #,
      # "M*" = listDataFramesB <-
      #   simulateB_MAster.parallel(
      #     DeltaList,
      #     GammaList,
      #     sim,
      #     sample.size,
      #     Time.Origin.STD,
      #     coal.events.times.simB
      #   )
      #   ,
      # "M0.8" = listDataFramesB <-
      #   simulateB_K.parallel(
      #     DeltaList,
      #     GammaList,
      #     sim,
      #     sample.size,
      #     Time.Origin.STD,
      #     coal.events.times.simB,
      #     0.8
      #   ),
      # "M2" = listDataFramesB <-
      #   simulateB_K.parallel(
      #     DeltaList,
      #     GammaList,
      #     sim,
      #     sample.size,
      #     Time.Origin.STD,
      #     coal.events.times.simB,
      #     2
      #   )
      # ,
      # "Hn/2" = listDataFramesHybrid <-
      #   simulateHybrid.parallel(
      #     DeltaList,
      #     GammaList,
      #     sim,
      #     sample.size,
      #     Time.Origin.STD,
      #     coal.events.times.simHybrid,
      #     nB_List[1],
      #     number.Fail.hybrid,
      #     1,
      #     2,
      #     number.ancestors.simHybrid,
      #     number.ancestors.Transition
      #   )
      #   ,
      #   "Hn/4" = listDataFramesHybrid <-
      #   simulateHybrid.parallel(
      #     DeltaList,
      #     GammaList,
      #     sim,
      #     sample.size,
      #     Time.Origin.STD,
      #     coal.events.times.simHybrid,
      #     nB_List[2],
      #     number.Fail.hybrid,
      #     2,
      #     2,
      #     number.ancestors.simHybrid,
      #     number.ancestors.Transition
      #   )
        ,
      times = number_runs
    )
  print(mbm)
  
  
  tiff(
    paste(path_to_save, "/",  "PlotBenchmark", ".tiff", sep = ""),
    units = "in",
    width = 3.25,
    height = 3.25,
    res = 300
  )
  plot_mbm <-
    ggplot2::autoplot(mbm, log = TRUE) +  theme_bw() + theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
  plot_mbm
  ggsave(paste(path_to_save, "/",  "plot_mbm_n=", sample.size,"_runs=", number_runs, "_sim=", sim, ".pdf", sep = ""))
  ggsave(paste(path_to_save, "/",  "plot_mbm_n=", sample.size,"_runs=", number_runs, "_sim=", sim, ".png", sep = ""), dpi=300, width=4, height=4)

  
  dev.off()
  
 if (length(nB_List)>1){
   total.failed <-
     bigstatsr::big_apply(
       number.Fail.hybrid,
       a.FUN = function(X, ind) {
         rowSums(X[, ind, drop = FALSE])
       },
       a.combine = 'plus',
       block.size = 10
     )
  
 }
  else{
    
    total.failed <-sum(X[, 1, drop = FALSE])
    
  }
  print(paste0("Total number failures in the Hybrid per k0 value in ", nB_List))
  print(total.failed)
  print(paste0("Number failures, rows correspond to k0 values ", nB_List," and columns Gamma values ", GammaList))
  print(number.Fail.hybrid[])
  
}
