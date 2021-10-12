#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tibble)
library(tidyr)
library(dplyr)
require(foreach)
require(doParallel)
require(doSNOW)
require(parallel)
require(doFuture)
require(bigstatsr)
require(stats)
require(rgenoud)
require(doRNG)
require(doMC)
library(shinyjs)
library(ggpubr)
library(ggtree)
library(ggplotify)
library(nLTT)
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
source(paste0(path,"/CoalSimulationBirthDeath.R"))
get_ratio_data<-function(number_sim, sample_size, Delta_list,Gamma_list,
                         list_data_fames_A, list_data_fames_B){
    
    stopifnot(length(Delta_list)==length(Gamma_list))

    
    data_frame_A<-do.call("rbind", list_data_fames_A)
    data_frame_B<-do.call("rbind", list_data_fames_B)
    
    delta_labs <- lapply(Delta_list, function(x, name) paste(name, "=", x, sep=""), "Delta")
    names(delta_labs)<-factor(Delta_list)
    gamma_labs <- lapply(Gamma_list, function(x, name) paste(name, "=", x, sep=""), "Gamma")
    names(gamma_labs)<- factor(Gamma_list)
    
    
    simulationColumn=rep(1:length(Delta_list), rep(sample_size, length(Delta_list)))
    Delta_column=rep(Delta_list, rep(sample_size-1, length(Delta_list)))
    Gamma_column=rep(Gamma_list, rep(sample_size-1, length(Gamma_list)))

    data_frame_A$Delta=Delta_column
    data_frame_B$Delta=Delta_column
    data_frame_A$GammaVar=Gamma_column
    data_frame_B$GammaVar=Gamma_column
    

    data_frame_A$x=rep((sample_size-1):1, length(list_data_fames_A))
    data_frame_B$x=rep((sample_size-1):1, length(list_data_fames_B))
    
    mean_ratio_colum<-data_frame_A$mean /  data_frame_B$mean
    median_ratio_colum<-data_frame_A$median /  data_frame_B$median
    ratio_data_frame<-data.frame(RatioMeanA1overA2=mean_ratio_colum,RatioMedianA1overA2 =median_ratio_colum, Delta=Delta_column,
                                 GammaVar=Gamma_column, x=rep((sample_size-1):1, length(list_data_fames_A)))
    return(ratio_data_frame)
    
}
getHTML <- function (frames) {
    innerhtml = '<div class="grid-container">'
    for (row in 1:(nrow(frames))) {
        id <- frames[row, "id"]
        name  <- frames[row, "names"]
        row_html = '<div class="grid-item">'
        row_html = paste(row_html, '<span>Name: ' , name, "id ", row , '</span>')
        row_html = paste(row_html, '</div>')
        
        innerhtml = paste(innerhtml, row_html)
    }
    paste(innerhtml, "</div>")
    return (innerhtml)
}

jsToggleFS <- 'shinyjs.toggleFullScreen = function() {
    var element = document.documentElement,
      enterFS = element.requestFullscreen || element.msRequestFullscreen || element.mozRequestFullScreen || element.webkitRequestFullscreen,
      exitFS = document.exitFullscreen || document.msExitFullscreen || document.mozCancelFullScreen || document.webkitExitFullscreen;
    if (!document.fullscreenElement && !document.msFullscreenElement && !document.mozFullScreenElement && !document.webkitFullscreenElement) {
      enterFS.call(element);
    } else {
      exitFS.call(document);
    }
  }'

build_tree_from_coal_times <-
  function(coal.events.times)
  {
    require(dplyr)
    require(ape)
    
    sample.size <- length(coal.events.times)+1
    
    coalescent_tree = rcoal(sample.size, br = "coalescent")
    x <- dplyr::as_tibble(coalescent_tree)
    
    df <- as.data.frame(x)
    
    node_times <- c(rep(0, sample.size), rev(coal.events.times))
    
    x[, "branch.length"] <-
      unlist(lapply(1:nrow(x), function(idx, node_times, x) {
        if (as.numeric(x[idx, "parent"]) != as.numeric(x[idx, "node"]))
          return(node_times[as.numeric(x[idx, "parent"])] - node_times[as.numeric(x[idx, "node"])])
        else
          return(NA)
        
      }, node_times, x))
    
    tree <- ape::as.phylo(x)
    
    return(tree)
    
  }
nltt_plot2 <- function( phy, xlab = "Normalized Time",
                        ylab = "Normalized Lineages", ...) {
  
  if (!inherits(phy, "phylo")) {
    # Stop imposed by ape::ltt.plot.coords
    stop(
      "phylogeny must be of class 'phylo', ",
      "but was of type '", class(phy), "' instead"
    )
  }
  
  #we use the ltt.plot.coords function from the package ape
  xy <- ape::ltt.plot.coords(phy, backward = TRUE, tol = 1e-6)
  xy[, 2] <- xy[, 2] / max(xy[, 2]) #normalize number lineages
  
  xy[, 1] <- xy[, 1] + abs( min( xy[, 1])) #make sure time runs from 0..T
  xy[, 1] <- xy[, 1] / max( xy[, 1])      #normalize time
  
  graphics::plot.default(xy, xlab = xlab, ylab = ylab, xaxs = "r", yaxs = "r",
                         type = "S", xlim= c(0,1), ylim=c(0,1),...) #type = "S" ensures a stepwise function
  graphics::abline(a=0, b=1,  lty=2)
 
}
ui <- fluidPage(

    useShinyjs(),
    extendShinyjs(text = jsToggleFS, functions = c("winprint")),
    
    titlePanel("Limit Birth Death Coalescent & approximations"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("sample_size",
                        "Sample size:",
                        min = 2,
                        max = 100,
                        value = 10),
            sliderInput("sim",
                        "Number of simulations:",
                        min = 10,
                        max = 1000,
                        value = 10),
            numericInput(inputId = "Gamma", label = HTML("&Gamma;:"), value = 10, step =  1 ),
            numericInput(inputId = "Delta", label = HTML("&Delta;:"), value = 1, step =  1 ),
            sliderInput("K",
                        "K for M_K approximation:",
                        min = 0,
                        max =2,
                        value = 0.8),
            sliderInput("m",
                        "m for Hybrid approximation:",
                        min = 1,
                        max = 100,
                        value = 5),
            tags$head(tags$script(src = "message-handler.js")),
            div(
                HTML("<button type='button'>Toggle Fullscreen</button>"),
                onclick = "shinyjs.toggleFullScreen();"
            ),
            actionButton("go", "Simulate"),
            tableOutput("values")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("nLTTPlot"),
          plotOutput("plotTree"),
           plotOutput("distPlot")
        )
    ),
    
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    v <- reactiveValues(doPlot = FALSE)
    
    observeEvent(input$go, {
        # 0 will be coerced to FALSE
        # 1+ will be coerced to TRUE
        v$doPlot <- input$go
       # session$sendCustomMessage(type = 'testmessage',
       #                           message = 'Simulating, please wait...')
    })
   
    observe({
        val <- input$sample_size
        # Control the value, min, max, and step.
        # Step size is 2 when input value is even; 1 when value is odd.
        updateSliderInput(session, "m", value = floor(val/2),
                          min = 1, max = val, step = 1)
    })
    
    
    
    model <- eventReactive(input$go, {
        print("simulating")
        path = getCurrentFileLocation()
        sim<- input$sim
        K<- input$K
        sample.size<-input$sample_size
        DeltaList = c(input$Delta)
        GammaList = c(input$Gamma)
        nB = input$m
        nB_List = c(nB)
        number.sim.trees = 5
        
      isolate({
        Time.Origin.STD <-
            bigstatsr::FBM(length(DeltaList), input$sim , type = "double", init = 0)
        
        coal.events.times.simA = bigstatsr::FBM(length(DeltaList), input$sim * (input$sample_size -
                                                                                    1))
        
        number.ancestors.simA = bigstatsr::FBM(length(DeltaList), input$sim * (input$sample_size -
                                                                                   1))
        
        Time.Origin.STD.trees <-
            bigstatsr::FBM(length(DeltaList),
                           number.sim.trees ,
                           type = "double",
                           init = 0)
        coal.events.times.simA.trees = bigstatsr::FBM(length(DeltaList), number.sim.trees *
                                                          (input$sample_size- 1))
        
        coal.events.times.simB = bigstatsr::FBM(length(DeltaList), input$sim * (input$sample_size -
                                                                                    1))
        
        coal.events.times.simHybrid = bigstatsr::FBM(length(DeltaList), input$sim * (input$sample_size -
                                                                                         1))
        
        
        
        coal.events.times.simB = bigstatsr::FBM(length(DeltaList), input$sim * (input$sample_size -
                                                                                    1))
        
        number.Fail.hybrid = bigstatsr::FBM(length(nB_List), length(GammaList))
        coal.events.times.simHybrid = bigstatsr::FBM(length(DeltaList), input$sim *
                                                         (input$sample_size - 1))
        number.ancestors.simHybrid = bigstatsr::FBM(length(DeltaList), input$sim * (input$sample_size -
                                                                                        1))
        number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), input$sim)
        
        innerCluster <-
            parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
        doParallel::registerDoParallel(innerCluster)
        
        parallel::clusterExport(innerCluster, varlist = c("path","K", "nB", "nB_List", "sim", "sample.size", "DeltaList","GammaList","number.sim.trees",  "session"), envir = environment())
        
        withProgress(message = 'Simulating....', value = 0, {
            # Number of times we'll go through the loop
            
            n =  3

            
            incProgress(1/n, detail = paste("..LBD coalescent"))
            RNGkind("L'Ecuyer-CMRG")
            a=as.numeric(Sys.time())
            set.seed(floor(runif(1)*a))
            
            rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
            
            opts <- list(chunkSize = 2)
          
            
            
            tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
                foreach::foreach(i = 1:length(DeltaList),
                                 r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                                 .combine = 'c', .export = c("isolate","input")) %dopar% {
                                     rngtools::setRNG(r)
                                     
                                     
                                     source(paste0(path,"/CoalSimulationBirthDeath.R"))
                                     list.number.ancestors.population = simulate.list.number.ancestors.population(sample.size)
                                     
                                     
                                     list.coal.times = simulate.coalescent.times.A(GammaList[i],
                                                                                   DeltaList[i],
                                                                                   sample.size,
                                                                                   list.number.ancestors.population)
                                     
                                     coal.events.times = unlist(list.coal.times)
                                     
                                     Time.Origin = coal.events.times[length(coal.events.times)]
                                     
                                     Time.Origin.STD[i, j] <- Time.Origin
                                     
                                     positions <-
                                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                                                 (sample.size - 1)))
                                     
                                     number.ancestors.simA[positions] <-
                                         list.number.ancestors.population[2:(length(coal.events.times))]
                                     
                                     coal.events.times.simA[positions] <-
                                         coal.events.times[1:(length(coal.events.times) - 1)]
                                     m = list.number.ancestors.population[floor(length(list.number.ancestors.population) /
                                                                                    2)]
                                     number.ancestors.Transition[cbind(i, j)] <- m
                                     
                                     
                                     NULL
                                 }
            
            list.Data.FramesA <-
                parallel::mclapply(
                    1:length(DeltaList),
                    mc.cores = parallel::detectCores() - 1,
                    FUN = function(k,
                                   sample.size,
                                   sim,
                                   coal.events.times.simA) {
                        quants <- c(0.025, 0.50, 0.975)
                        
                        matrixCurrentValue <-
                            matrix(
                                coal.events.times.simA[k, ],
                                nrow = sim ,
                                ncol = sample.size - 1,
                                byrow = TRUE
                            )
                        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
                        quantiles <-
                            apply(matrixCurrentValue , 2 , quantile , probs = quants)
                        meanCoalTimes <- colMeans(matrixCurrentValue)
                        result <-
                            data.frame(
                                mean = meanCoalTimes,
                                LI = quantiles[1, ],
                                median = quantiles[2, ],
                                UI = quantiles[3, ]
                            )
                        result
                    },
                    sample.size,
                    sim,
                    coal.events.times.simA
                )

            
            Sys.sleep(0.1)
            incProgress(1/n, detail = paste("..M_k approximation"))
         
            tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
                foreach::foreach(i = 1:length(DeltaList),
                                 r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                                 .combine = 'c') %dopar% {
                                     
                                     source(paste0(path,"/CoalSimulationBirthDeath.R"))
                                     rngtools::setRNG(r)
                                     
                                     TimeOrigin = Time.Origin.STD[i, j]
                                     
                                     if (K>0){
                                         
                                         coal.events.times = modelCoalB_K(sample.size, TimeOrigin, DeltaList[i], GammaList[i], K)
                                     }
                                     else{#K=0
                                         coal.events.times = modelCoal(sample.size, TimeOrigin, DeltaList[i], GammaList[i])
                                     }
                                     
                                     print(
                                         paste(
                                             "finished B M_k in parallel for Delta",
                                             DeltaList[i],
                                             "Gamma",
                                             GammaList[i],
                                             "K",
                                             K,
                                             " sim ",
                                             j,
                                             sep = " "
                                         )
                                     )
                                     positions <-
                                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                                                 (sample.size - 1)))
                                     coal.events.times.simB[positions] <- rev(coal.events.times)
                                     NULL
                                 }
            list.Data.FramesB <-
                parallel::mclapply(
                    1:length(DeltaList),
                    mc.set.seed = TRUE,
                    mc.cores = parallel::detectCores() - 1,
                    FUN = function(k,
                                   sample.size,
                                   sim,
                                   coal.events.times.simB)
                    {
                        quants <- c(0.025, 0.50, 0.975)
                        matrixCurrentValue <-
                            matrix(
                                coal.events.times.simB[k, ],
                                nrow = sim ,
                                ncol = sample.size - 1,
                                byrow = TRUE
                            )
                        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
                        quantiles <-
                            apply(matrixCurrentValue , 2 , quantile , probs = quants)
                        meanCoalTimes <- colMeans(matrixCurrentValue)
                        result <-
                            data.frame(
                                mean = meanCoalTimes,
                                LI = quantiles[1, ],
                                median = quantiles[2, ],
                                UI = quantiles[3, ]
                            )
                        result
                    },
                    sample.size,
                    sim,
                    coal.events.times.simB
                )

            posMprime = 1
            Sys.sleep(0.1)
            incProgress(1/n, detail = paste("...Hybrid approximation"))
            
          
            tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
                foreach::foreach(i = 1:length(DeltaList),
                                 r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                                 .combine = 'c') %dopar% {
                                     
                                     source(paste0(path,"/CoalSimulationBirthDeath.R"))
                                     rngtools::setRNG(r)
                                     posMprime = 1
                                     
                                     TimeOrigin <- Time.Origin.STD[i, j]
                                     
                                     m = -1
                                     while (m < nB) {
                                         if (m != -1) {
                                             number.Fail.hybrid[cbind(posMprime, i)] = number.Fail.hybrid[cbind(posMprime, i)] +
                                                 1
                                         }
                                         
                                         std.coal.events.times.until.B = rev(modelCoalHybrid(sample.size, TimeOrigin, DeltaList[i], GammaList[i], nB))
                                         
                                         coal.events.times.until.nB = unlist(lapply(
                                             std.coal.events.times.until.B,
                                             FUN = function(x)
                                                 standard2modelB_K(x, TimeOrigin, DeltaList[i], GammaList[i], K)
                                         ))
                                         
                                         t = std.coal.events.times.until.B[length(std.coal.events.times.until.B)]
                                         
                                         m = ceiling(2.0 / t)
                                         
                                         
                                     }
                                     number.ancestors.Transition[cbind(i, j)] <- m
                                     
                                     
                                     list.number.ancestors.population = simulate.list.number.ancestors.population.from(sample.size, m, nB)
                                     
                                     s <-
                                         standard2modelB_K(t, TimeOrigin, DeltaList[i], GammaList[i], K)
                                     
                                     
                                     list.coal.times = simulate.coalescent.times.A.from(GammaList[i],
                                                                                        DeltaList[i],
                                                                                        sample.size,
                                                                                        list.number.ancestors.population,
                                                                                        nB,
                                                                                        s)
                                     
                                     coal.events.times = c(unlist(coal.events.times.until.nB),
                                                           unlist(list.coal.times))
                                     
                                     
                                     positions <-
                                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                                                 (sample.size - 1)))
                                     
                                     coal.events.times.simHybrid[positions] <-
                                         unlist(coal.events.times)
                                     
                                     number.ancestors.simHybrid[positions] <-
                                         c(rep(0, sample.size - nB - 2),
                                           unlist(list.number.ancestors.population),
                                           0)
                                     
                                     
                                     print(
                                         paste(
                                             "finished hybrid scenario in parallel for Delta",
                                             DeltaList[i],
                                             "Gamma",
                                             GammaList[i],
                                             " sim ",
                                             j,
                                             sep = " "
                                         )
                                     )
                                     NULL
                                 }
            
            list.Data.FramesHybrid <-
                parallel::mclapply(
                    1:length(DeltaList),
                    mc.cores = parallel::detectCores() - 1,
                    FUN = function(k,
                                   sample.size,
                                   sim,
                                   coal.events.times.simHybrid) {
                        quants <- c(0.025, 0.50, 0.975)
                        
                        matrixCurrentValue <-
                            matrix(
                                coal.events.times.simHybrid[k, ],
                                nrow = sim ,
                                ncol = sample.size - 1,
                                byrow = TRUE
                            )
                        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
                        quantiles <-
                            apply(matrixCurrentValue , 2 , quantile , probs = quants)
                        meanCoalTimes <- colMeans(matrixCurrentValue)
                        result <-
                            data.frame(
                                mean = meanCoalTimes,
                                LI = quantiles[1, ],
                                median = quantiles[2, ],
                                UI = quantiles[3, ]
                            )
                        result
                    },
                    sample.size,
                    sim,
                    coal.events.times.simHybrid
                )
            registerDoSEQ()
            stopCluster(innerCluster)    
            
            ratio.data.frame_M_K<-get_ratio_data(sim,sample.size, DeltaList,GammaList, list.Data.FramesA,list.Data.FramesB )
            ratio.data.frame_M_K$Model<-paste("M_",K,"/LBD", sep= "")
            ratio.data.frame_M_K$Sample<-as.factor(sample.size)
            
            
            ratio.data.frame_Hybrid<-get_ratio_data(sim,sample.size, DeltaList,GammaList, list.Data.FramesA,list.Data.FramesHybrid )
            ratio.data.frame_Hybrid$Model<-paste("H_",nB,"/LBD", sep= "")
            ratio.data.frame_Hybrid$Sample<-as.factor(sample.size)
            
            ratio.data.frame_all <-rbind(ratio.data.frame_M_K, ratio.data.frame_Hybrid)
            
            ratio.data.frame_M_K$Sample<-as.factor(sample.size)
            
            ratio.data.frame_all<-unite(ratio.data.frame_all, ModelSample, Model:Sample, sep='', remove=F)
            print(ratio.data.frame_all)
            
            list(ratio.data.frame.all =ratio.data.frame_all,LBD_coal_times = as.data.frame(coal.events.times.simA[]), LBD_num_ancestors = as.data.frame(number.ancestors.simA[]), 
                 listDataFramesA = list.Data.FramesA,listDataFramesMK = list.Data.FramesB,listDataFramesHybrid = list.Data.FramesHybrid  )
         })
       })
    })    
     
    
    output$distPlot <- renderPlot({
        #print(v$doPlot)
        #if (v$doPlot == FALSE) return()
       
        sim<- input$sim
        K<- input$K
        sample.size<-input$sample_size
        DeltaList = c(input$Delta)
        GammaList = c(input$Gamma)
        nB = input$m
        nB_List = c(nB)
        number.sim.trees = 5
        
        mod_list=model()
        ratio.data.frame.all=mod_list$ratio.data.frame.all
       
            ggplot(ratio.data.frame.all, aes(x = x, y = RatioMedianA1overA2, colour=ModelSample)) +
                geom_point(size=1) +
                geom_hline(data = NULL, aes(yintercept = 1))+
                scale_x_continuous(breaks =1:sample.size)+
                ylim(min(ratio.data.frame.all$RatioMedianA1overA2),max(ratio.data.frame.all$RatioMedianA1overA2))+
                labs(y="", x ="Number of sample ancestors ")+
                theme(plot.margin = unit(c(1,1,2,1), "lines")) +
                #facet_grid(~ interaction(paste0("\u0393=",GammaVar), paste0("n=",Sample)), scales='free')+
                #facet_wrap(paste0("\u0393=",GammaVar) ~ paste0("n=",Sample), scales='free')+
                facet_wrap(~ interaction(paste0("\u0393=",GammaVar), paste0(" n=",Sample)), scales='free')+
                #coord_capped_cart(bottom='both', left='both', xlim=c(1, 100))+
                theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                theme( legend.position="bottom")+
                theme(legend.text=element_text(size=12), legend.title=element_text(size=12))+
                theme(legend.key.size = unit(0.5, "cm"), axis.title.y=element_text(size=12))+
                scale_fill_brewer(palette="Accent")
       #   })
       # })
    
    })
    
    output$nLTTPlot<-renderPlot({
      
      sim<- input$sim
      K<- input$K
      sample.size<-input$sample_size
      DeltaList = c(input$Delta)
      GammaList = c(input$Gamma)
      nB = input$m
      nB_List = c(nB)
      number.sim.trees = 5
      
      mod_list=model()
      listDataFramesA=mod_list$listDataFramesA
    
      list_trees_mean_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(listDataFramesA[[i]]$mean)
        
      },
      listDataFramesA)
      
      list_trees_median_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(listDataFramesA[[i]]$median)
        
      },
      listDataFramesA)
      
      list_trees_log_mean_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(log(listDataFramesA[[i]]$mean))
        
      },
      listDataFramesA)
      
      list_trees_log_median_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(log(listDataFramesA[[i]]$median))
        
      },
      listDataFramesA)
      
      class(list_trees_mean_coal_times)<-"multiPhylo"
      
      class(list_trees_median_coal_times)<-"multiPhylo"
      
      class(list_trees_log_mean_coal_times)<-"multiPhylo"
      
      class(list_trees_log_median_coal_times)<-"multiPhylo"
      
      tree<- as.phylo(list_trees_mean_coal_times[[1]])
      
      xlab = "Normalized Time"
      ylab = "Normalized Lineages"
      
      xy <- ape::ltt.plot.coords(tree, backward = TRUE, tol = 1e-6)
      xy[, 2] <- xy[, 2] / max(xy[, 2]) #normalize number lineages
      
      xy[, 1] <- xy[, 1] + abs( min( xy[, 1])) #make sure time runs from 0..T
      xy[, 1] <- xy[, 1] / max( xy[, 1])      #normalize time
      
      graphics::plot.default(xy, xlab = xlab, ylab = ylab, xaxs = "r", yaxs = "r",
                             type = "S", xlim= c(0,1), ylim=c(0,1), main = "nLTT with average coalescent times in LBD Coalescent") #type = "S" ensures a stepwise function
      graphics::abline(a=0, b=1,  lty=2)
      
    #   p1<-nltt_plot2(list_trees_mean_coal_times[[1]], lwd=2)
    #   
    #   tree<- as.phylo(tree)
    #   print(tree)
    #   #p2<-ggplotify::as.ggplot(~nLTT::nltt_plot(tree, lwd=2))
    #   
    #   t1<-ggtree(list_trees_mean_coal_times[[1]])
    #  
    # #  p1<-ggplotify::as.ggplot(~nltt_plot2(list_trees_mean_coal_times[[1]], lwd=2))
    #   
    #   
    #   ggpubr::ggarrange(t1,p1,
    #                     labels=c(  "Mean coal times topology", 
    #                                "NLTT"),
    #                     ncol = 2, nrow = 1)
    })
    
    output$plotTree<-renderPlot({
      
      sim<- input$sim
      K<- input$K
      sample.size<-input$sample_size
      DeltaList = c(input$Delta)
      GammaList = c(input$Gamma)
      nB = input$m
      nB_List = c(nB)
      number.sim.trees = 5
      
      mod_list=model()
      listDataFramesA=mod_list$listDataFramesA
      
      list_trees_mean_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(listDataFramesA[[i]]$mean)
        
      },
      listDataFramesA)
      
      list_trees_median_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(listDataFramesA[[i]]$median)
        
      },
      listDataFramesA)
      
      list_trees_log_mean_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(log(listDataFramesA[[i]]$mean))
        
      },
      listDataFramesA)
      
      list_trees_log_median_coal_times <- lapply(1:length(listDataFramesA), function(i, listDataFramesA){
        build_tree_from_coal_times(log(listDataFramesA[[i]]$median))
        
      },
      listDataFramesA)
      
      class(list_trees_mean_coal_times)<-"multiPhylo"
      
      class(list_trees_median_coal_times)<-"multiPhylo"
      
      class(list_trees_log_mean_coal_times)<-"multiPhylo"
      
      class(list_trees_log_median_coal_times)<-"multiPhylo"
      
      tree<- as.phylo(list_trees_mean_coal_times[[1]])
      
      ggtree(tree)+ geom_treescale()+ ggtitle("Random topology with avg coal times of LBD Coalescent ") 
   
    })
    sliderValues <- reactive({
        
        data.frame(
            Name = c("Sample size",
                     "Number simulations",
                     "Gamma",
                     "Delta",
                     "K",
                     "m"
                     ),
            Value = as.character(c(input$sample_size,
                                   input$sim,
                                   input$Gamma,
                                   input$Delta,
                                   input$K,
                                   input$m)),
            stringsAsFactors = FALSE)
        
    })
    
    # Show the values in an HTML table ----
    output$values <- renderTable({
        sliderValues()
        #actionButton('goSimulate', 'Simulate')
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
