rm(list=ls())

library(GGally)
library(ggcorrplot)
library(grid)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tibble)

#Gets the path where the .R file is executed
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
expected_number_ancestors =function(sample.size){
  
  list.number.ancestors<-lapply((sample.size-2):0, FUN=function(k, sample.size)
  {
    
    sample.size*k/(sample.size-k-1)
  }, sample.size)
  
  return(unlist(list.number.ancestors))
}
###########################################################################################
sample.size <- 100
sim = 10000
nB<- floor(sample.size / 2)
GammaList = c(0.001,  0.1,  10)
DeltaList= rep(1, length(GammaList))
K = 0.8
path_to_save<-getCurrentFileLocation()


#read the data
list_corr_mat_hybrid<-readRDS( paste(path_to_save,"/", "list_corr_mat_Hybrid","_","K=", K,"_",sample.size,"_",nB,"_", sim,".rds", sep=""))
list_corr_matA<-readRDS(paste(path_to_save,"/", "list_corr_mat_A_","_",sample.size,"_", sim,".rds", sep=""))
list_corr_mat_B_K<-readRDS(paste(path_to_save,"/", "list_corr_mat_B_","_",sample.size,"_K=",K,"_", sim, ".rds", sep=""))

listDataFramesA<-readRDS(paste(path_to_save,"/", "listDataFramesA_",sample.size,"_", sim,".rds", sep=""))
listDataFramesB<-readRDS(paste(path_to_save,"/","listDataFramesB_",sample.size,"_", sim,".rds", sep=""))

listDataFramesB_K<-readRDS(file=paste(path_to_save,"/","listDataFramesB_",sample.size,"_K=",K, "_", sim,".rds", sep=""))
listDataFramesHybrid<-readRDS(paste(path_to_save,"/","listDataFramesHybrid_",sample.size,"_",nB,"_", sim,".rds", sep=""))

listDataFramesAncestorsA<-readRDS(paste(path_to_save,"/", "listDataFramesA_NumberAncestors",sample.size,"_", sim, ".rds", sep=""))
listDataFramesAncestorsHybrid<-readRDS(paste(path_to_save,"/", "listDataFramesHybrid_NumberAncestors",sample.size,"_", sim,".rds", sep=""))

number.ancestors.Transition.Hybrid<-readRDS(paste(path_to_save,"/","number.ancestors.Transition_","K=", K, "_",sample.size,"_",nB,"_", sim,".rds", sep=""))
number.ancestors.Transition.A<-readRDS(paste(path_to_save,"/","number.ancestors.Transition_", sample.size,"_", sim,".rds", sep=""))


dataframeNumberAncestorsHybrid<-as.data.frame(t(number.ancestors.Transition.Hybrid))
print(colMeans(as.matrix(dataframeNumberAncestorsHybrid)))

number_ancestors= c(dataframeNumberAncestorsHybrid$V1, dataframeNumberAncestorsHybrid$V2, dataframeNumberAncestorsHybrid$V3)
Gamma<- rep(GammaList, rep(sim, length(GammaList)))
dataFrameAllNumbersAncestors<-data.frame(number_ancestors,
                                          Gamma)
dataFrameAllNumbersAncestors$Gamma<-as.factor(dataFrameAllNumbersAncestors$Gamma)
dataFrameAllNumbersAncestors$model <-"Hybrid"
dataFrameAllNumbersAncestors$model<-as.factor(dataFrameAllNumbersAncestors$model)

nB=floor(sample.size/2)
Expected <- rep(sample.size*nB/(sample.size-nB-1), length(GammaList))
Gamma<-levels(dataFrameAllNumbersAncestors$Gamma)
dataTheorModelGamma <- data.frame(Gamma,  Expected)
#dataTheor$Gamma <-as.factor(dataTheor$Gamma)
dataTheorModelGamma$Gamma <-as.factor(dataTheorModelGamma$Gamma)

model<-c("Hybrid", "BD")
Expected<-rep(sample.size*nB/(sample.size-nB-1), length(model))
dataTheorModel <-data.frame(model,  Expected)
#dataTheor$Gamma <-as.factor(dataTheor$Gamma)
dataTheorModel$model <-as.factor(dataTheorModel$model)



dataframeNumberAncestorsBD<-as.data.frame(t(number.ancestors.Transition.A))
print(colMeans(as.matrix(dataframeNumberAncestorsBD)))
number_ancestors= c(dataframeNumberAncestorsBD$V1, dataframeNumberAncestorsBD$V2, dataframeNumberAncestorsBD$V3)
Gamma<- rep(GammaList, rep(sim, length(GammaList)))
dataFrameAllNumbersAncestorsBD<-data.frame(number_ancestors,
                                         Gamma)
dataFrameAllNumbersAncestorsBD$Gamma<-as.factor(dataFrameAllNumbersAncestorsBD$Gamma)
dataFrameAllNumbersAncestorsBD$model <-"BD"

dataFrameAllNumbersAncestors2Models<-rbind(dataFrameAllNumbersAncestors, dataFrameAllNumbersAncestorsBD)
dataFrameAllNumbersAncestors2Models$Model_Gamma <- paste(dataFrameAllNumbersAncestors2Models$model,"-", dataFrameAllNumbersAncestors2Models$Gamma)
dataFrameAllNumbersAncestors2Models$Model_Gamma<-as.factor(dataFrameAllNumbersAncestors2Models$Model_Gamma)

#plot number of ancestor distribution
library(ggplot2)
# density for Hybrid 
p <- ggplot(dataFrameAllNumbersAncestors, aes(x=number_ancestors, color=Gamma)) + 
  geom_density()+
  geom_vline(data=dataTheorModelGamma, aes(xintercept=Expected, color=Gamma),
             linetype="dashed")+
  labs(title="Hybrid model", y="Density", x =paste("Number of ancestors when sample size =", nB))+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
        axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
        plot.title = element_text(size=15, face="bold", vjust=2),
        axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
        axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
        legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p

ggsave(paste(path_to_save, "/",  "DistributionNumberAncestorsHybrid_", sample.size, ".pdf", sep=""))
ggsave(paste(path_to_save, "/", "DistributionNumberAncestorsHybrid_", sample.size, ".png", sep=""))
ggsave(paste(path_to_save, "/", "DistributionNumberAncestorsHybrid_", sample.size, ".jpg", sep=""))
#comparison of Hybrid vs BD
p <- ggplot(dataFrameAllNumbersAncestors2Models, aes(x=number_ancestors, color=model)) + 
  geom_density()+
  geom_vline(data=dataTheorModel, aes(xintercept=Expected, color=model),
             linetype="dashed")+
  labs(title="Hybrid model vs BD", y="Density", x =paste("Number of ancestors when sample size =", nB))+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
        axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
        plot.title = element_text(size=15, face="bold", vjust=2),
        axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
        axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
        legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p

ggsave(paste(path_to_save, "/",  "DistributionNumberAncestorsHybridvsBD_", sample.size, ".pdf", sep=""))
ggsave(paste(path_to_save, "/", "DistributionNumberAncestorsHybridvsBD_", sample.size, ".png", sep=""))
ggsave(paste(path_to_save, "/", "DistributionNumberAncestorsHybridvsBD_", sample.size, ".jpg", sep=""))

Model_Gamma<-levels(dataFrameAllNumbersAncestors2Models$Model_Gamma)
Expected <- rep(sample.size*nB/(sample.size-nB-1), length(Model_Gamma))
dataTheorModel=data.frame(Model_Gamma, Expected)

p <- ggplot(dataFrameAllNumbersAncestors2Models, aes(x=number_ancestors, color=Model_Gamma )) + 
  geom_density()+
  geom_vline(data=dataTheorModel, aes(xintercept=Expected, color=Model_Gamma),
             linetype="dashed")+
  labs(title="Hybrid model vs BD", y="Density", x =paste("Number of ancestors when sample size =", nB))+
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
        axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
        plot.title = element_text(size=15, face="bold", vjust=2),
        axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
        axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
        legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p

ggsave(paste(path_to_save, "/",  "DistributionNumberAncestorsHybridvsBD2_", sample.size, ".pdf", sep=""))
ggsave(paste(path_to_save, "/", "DistributionNumberAncestorsHybridvsBD2_", sample.size, ".png", sep=""))
ggsave(paste(path_to_save, "/", "DistributionNumberAncestorsHybridvsBD2_", sample.size, ".jpg", sep=""))
###############################################################################################################
#make the plot


rownames(list_corr_matA[[1]]) <- as.character(1:9)
colnames(list_corr_matA[[1]]) <- as.character(1:9)
p1<-ggcorrplot(as.matrix(list_corr_matA[[1]]), hc.order = FALSE,
               lab = FALSE, digits = 2, lab_size = 2,show.diag=TRUE, show.legend=FALSE,type = "upper")
                                                                                                   #,
                                                                                                   #ggtheme= theme(axis.title.x=element_blank(),
                                                                                                    #  axis.text.x=element_blank(),
                                                                                                    #  axis.ticks.x=element_blank()))

rownames(list_corr_mat_hybrid[[1]]) <- as.character(1:9)
colnames(list_corr_mat_hybrid[[1]]) <- as.character(1:9)

p2<-ggcorrplot(as.matrix(list_corr_mat_hybrid[[1]]), hc.order = FALSE,
               lab = FALSE, digits = 2, lab_size = 2,show.legend=FALSE ,show.diag=TRUE,type = "upper")

rownames(list_corr_matA[[2]]) <- as.character(1:9)
colnames(list_corr_matA[[2]]) <- as.character(1:9)
p3<-ggcorrplot(as.matrix(list_corr_matA[[2]]), hc.order = FALSE, 
               lab = FALSE, digits = 2, lab_size = 2, show.legend=FALSE,show.diag=TRUE,type = "upper")

rownames(list_corr_mat_hybrid[[2]]) <- as.character(1:9)
colnames(list_corr_mat_hybrid[[2]]) <- as.character(1:9)

p4<-ggcorrplot(as.matrix(list_corr_mat_hybrid[[2]]), hc.order = FALSE, 
               lab = FALSE, digits = 2, lab_size = 2, show.legend=FALSE,show.diag=TRUE,type = "upper")

rownames(list_corr_matA[[3]]) <- as.character(1:9)
colnames(list_corr_matA[[3]]) <- as.character(1:9)

p5<-ggcorrplot(as.matrix(list_corr_matA[[3]]), hc.order = FALSE, 
               lab = FALSE,digits = 2, lab_size = 2, show.legend=FALSE,show.diag=TRUE, type = "upper")

rownames(list_corr_mat_hybrid[[3]]) <- as.character(1:9)
colnames(list_corr_mat_hybrid[[3]]) <- as.character(1:9)

p6<-ggcorrplot(as.matrix(list_corr_mat_hybrid[[3]]), hc.order = FALSE, 
               lab = FALSE,digits = 2, lab_size = 2, show.legend=FALSE,show.diag=TRUE,type = "upper")



figure <- ggarrange(p2, p1, p4,p3,p6, p5,
                    labels = c( text_grob(bquote("H \u0393=10"^-4)),
                                text_grob(bquote("BD \u0393=10"^-4)), 
                                text_grob(bquote("H \u0393=10"^-1)),
                                text_grob(bquote("BD \u0393=10"^-1)),
                                text_grob(bquote("H \u0393=10")),
                                text_grob(bquote("BD \u0393=10"))),
                                ncol = 2, nrow = 3)
figure
annotate_figure(figure,
                top = text_grob(bquote("Correlations "*Hybrid~vs~BD~coal*" "), color = "black", face = "bold", size = 14),
                fig.lab.face = "bold"
)


ggsave(paste(path_to_save, "/",  "CorrelationsHybrid_",sample.size, ".pdf", sep=""), width = 11, height = 9 )
ggsave(paste(path_to_save, "/", "CorrelationsHybrid_",sample.size, ".png", sep=""), width = 11, height = 9)
ggsave(paste(path_to_save, "/", "CorrelationsHybrid_",sample.size, ".jpg", sep=""), width = 11, height = 9)

dev.off()

print(paste0(" Max % of correlation diff between Hybrid and BD ","\u0393=",0.001))
print(max(abs(list_corr_matA[[1]]-list_corr_mat_hybrid[[1]])*100/list_corr_matA[[1]]))
print(paste0(" Max % of correlation diff between Hybrid and BD ","\u0393=",0.1))
print(max(abs(list_corr_matA[[2]]-list_corr_mat_hybrid[[2]])*100/list_corr_matA[[2]]))
print(paste0(" Max % of correlation diff between Hybrid and BD ","\u0393=",10))
print(max(abs(list_corr_matA[[3]]-list_corr_mat_hybrid[[3]])*100/list_corr_matA[[3]]))

########################################################################################################
data_frame_Pop_Size<-do.call("rbind", listDataFramesAncestorsA)
Gamma_column=rep(GammaList, rep(sample.size-1, length(GammaList)))
data_frame_Pop_Size$Gamma= Gamma_column
data_frame_Pop_Size$t<-rep(((sample.size-2):0), length(listDataFramesAncestorsA))

t<-c((sample.size-2):0)
expected<-expected_number_ancestors(sample.size)
#plot(t, N, type = "o", pch = 19, las = 1, ylab= " Expected number of ancestors", xlab="sample size")

N<-matrix(c(data_frame_Pop_Size$mean, expected),nrow=sample.size-1, ncol= length(listDataFramesAncestorsA)+1,byrow=FALSE)

tiff(paste(path_to_save, "/",  "NumberAncestorsPopulation=",sample.size,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)

library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)
library(primer)
library(RColorBrewer)
colors = brewer.pal(5, "Set1")
matplot(t, N, type = "o", las = 1, pch = 1:5, col = colors, lty = 1:5,  ylab= "Number of ancestors", 
        xlab="sample size", cex.axis=0.6,cex.lab=0.6,
        axis.break=seq(0,(sample.size-2)))
legend("topleft", legend = c("0.001", "0.1", "10", "theor."), 
       title = "\u0393", pch = 1:5, lty = 1:5, col = colors, cex = 0.6)
theme(plot.margin = unit(c(1,1,2,1), "lines")) 

ggsave(paste(path_to_save, "/",  "NumberAncestorsn=",sample.size,".pdf", sep=""))
ggsave(paste(path_to_save, "/", "NumberAncestorsn=",sample.size,".png", sep=""))
ggsave(paste(path_to_save, "/", "NumberAncestorsn=",sample.size,".jpg", sep=""))

dev.off()
  