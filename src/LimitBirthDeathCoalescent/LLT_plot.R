
rm(list = ls())
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
################################################################################################################
#' build_tree_from_coal_times
#' gets a simulate tree with some  set of parameter value
#' @param coal.events.times is a vector of cola times 
#' @param parameter_index index of parameter value
#' @param number.sim.trees: number of simulations
#' @param sample.size: the sample size used
#' @output list_trees: the list of simulated trees
################################################################################################################
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
##########################################################
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
#############################################################
path = getCurrentFileLocation()


library(phytools)
sim = 10000
sample.size = 10
GammaList = c(0.001,  0.1,  10, 1000)
DeltaList = rep(1, length(GammaList))



listDataFramesA<-readRDS(
   paste(
    path,
    "/",
    "listDataFramesA_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)
 
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


ltt_plot_mean_coal_times <- phytools::ltt(list_trees_mean_coal_times,log=TRUE)
ltt_plot_mean_coal_times

names(list_trees_mean_coal_times) <-c(  "\u0393=0.001", 
                                        "\u0393=0.1",
                                        "\u0393=10",
                                        "\u0393=1000")

names(list_trees_log_mean_coal_times) <-c(  "\u0393=0.001", 
                                        "\u0393=0.1",
                                        "\u0393=10",
                                        "\u0393=1000")

ape::mltt.plot(list_trees_log_mean_coal_times, backward = FALSE, legend = TRUE, dlty=TRUE, dcol=FALSE)

ape::mltt.plot(list_trees_mean_coal_times, backward = FALSE, legend = TRUE, dlty=TRUE, dcol=FALSE)


ape::mltt.plot(list_trees_mean_coal_times, backward = TRUE, legend = TRUE, dlty=TRUE, dcol=FALSE, log="y")


library(ggtree)
layout(matrix(1:8, 2, 4))
plot(list_trees_mean_coal_times[[1]], show.tip.label = FALSE,lwd=2)
title( c(expression(paste(Gamma, " = ", 0.001))), line = -2)
ltt.plot(phy=list_trees_mean_coal_times[[1]], ylab = "n(t)",lwd=2,xlim=c(-6,0))
ltt.lines(phy=list_trees_mean_coal_times[[1]])
title( c(expression(paste(Gamma, " = ", 0.001))), line = -2,lwd=2,)
plot(list_trees_mean_coal_times[[2]], show.tip.label = FALSE)
title( c(expression(paste(Gamma, " = ", 0.1))), line = -2)
ltt.plot(phy=list_trees_mean_coal_times[[2]], ylab = "n(t)", lwd=2,xlim=c(-6,0))
ltt.lines(phy=list_trees_mean_coal_times[[2]])
title( c(expression(paste(Gamma, " = ", 0.1))), line = -2, lwd=2)
plot(list_trees_mean_coal_times[[3]], show.tip.label = FALSE)
title( c(expression(paste(Gamma, " = ", 10))), line = -2)
ltt.plot(phy=list_trees_mean_coal_times[[3]], ylab = "n(t)", lwd=2,xlim=c(-6,0))
ltt.lines(phy=list_trees_mean_coal_times[[3]])
title( c(expression(paste(Gamma, " = ", 10))), line = -2,lwd=2)
plot(list_trees_mean_coal_times[[4]], show.tip.label = FALSE)
title( c(expression(paste(Gamma, " = ", 1000))), line = -2)
ltt.plot(phy=list_trees_mean_coal_times[[4]], ylab = "n(t)", lwd=2,xlim=c(-6,0))
ltt.lines(phy=list_trees_mean_coal_times[[4]])
title( c(expression(paste(Gamma, " = ", 1000))), line = -2, lwd=2)
layout(1)


library(ggtree)
library(ggpubr)
library(treedater)
library(nLTT)
library("grid")
library("ggplotify")
library(cowplot)
library(ggplot2)
library(ggimage)
library(gridExtra)

tiff(file=paste0(path, "/nLLTs_circular.tiff"),
     width=8, height=8, units="in", res=100)

t1<-ggtree(list_trees_mean_coal_times[[1]], layout='circular')
#title( c(expression(paste(Gamma, " = ", 0.001))), line = -2)
r1<-image(nLTT::nltt_plot(list_trees_mean_coal_times[[1]]))

p1<-as.ggplot(~nLTT::nltt_plot(list_trees_mean_coal_times[[1]]))

#title( c(expression(paste(Gamma, " = ", 0.001))), line = -2,lwd=2,)
t2<-ggtree(list_trees_mean_coal_times[[2]], layout='circular')
#title( c(expression(paste(Gamma, " = ", 0.1))), line = -2)
p2<-as.ggplot(~nLTT::nltt_plot(list_trees_mean_coal_times[[2]]))

#title( c(expression(paste(Gamma, " = ", 0.1))), line = -2, lwd=2)
t3<-ggtree(list_trees_mean_coal_times[[3]], layout='circular')
#title( c(expression(paste(Gamma, " = ", 10))), line = -2)
p3<-as.ggplot(~nLTT::nltt_plot(list_trees_mean_coal_times[[3]]))

#title( c(expression(paste(Gamma, " = ", 10))), line = -2,lwd=2)
t4<-ggtree(list_trees_mean_coal_times[[4]], layout='circular')
#title( c(expression(paste(Gamma, " = ", 1000))), line = -2)
p4<-as.ggplot(~nLTT::nltt_plot(list_trees_mean_coal_times[[4]]))
#title( c(expression(paste(Gamma, " = ", 1000))), line = -2, lwd=2)

#treedater::plot(x, t0 = NA, res = 100, ggplot = TRUE,
#               cumulative = FALSE, ...)

ggpubr::ggarrange(t1, t2, t3, t4, p1,p2,p3, p4,
                  labels=c(  "\u0393=0.001", 
                             "\u0393=0.1",
                             "\u0393=10",
                             "\u0393=1000",
                             "\u0393=0.001", 
                             "\u0393=0.1",
                             "\u0393=10",
                             "\u0393=1000"),
                              ncol = 4, nrow = 2)


dev.off()

ggsave(paste0(path, "/nLLTs_circular.png"))

cowplot::plot_grid(plotlist = plist,nrow = 2, ncol=4)

ggsave(paste0(path, "/nLLTs2.tiff"))




ape::ltt.plot(list_trees_mean_coal_times, xlab = "Time", ylab = "N",
         backward = TRUE, tol = 1e-6)
legend(-20, 10, lwd = c(2, 1), lty = c(1, 2), bty = "n",
       legend = c("Bird orders", "Random (coalescent) trees"))
ape::ltt.lines()

#########################################################

tiff(file=paste0(path, "/nLLTs_tiff.tiff"),
     width=8, height=6, units="in", res=100)


p1<-as.ggplot(~nltt_plot2(list_trees_mean_coal_times[[1]], lwd=2))


p2<-as.ggplot(~nltt_plot2(list_trees_mean_coal_times[[2]], lwd=2))

p3<-as.ggplot(~nltt_plot2(list_trees_mean_coal_times[[3]], lwd=2))

p4<-as.ggplot(~nltt_plot2(list_trees_mean_coal_times[[4]], lwd=2))

ggpubr::ggarrange(p1,p2,p3, p4,
                  labels=c(  "\u0393=0.001", 
                             "\u0393=0.1",
                             "\u0393=10",
                             "\u0393=1000"),
                  ncol = 4, nrow = 1)


