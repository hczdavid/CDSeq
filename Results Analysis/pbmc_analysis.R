#----------------------#
# PBMC results analysis
# Author: David Huang
#----------------------#

library(tidyverse)
library(corrplot)

load_obj <- function(f){
  env <- new.env();nm <- load(f, env)[1]
  env[[nm]]
}

# load data with ground truth
load("readyData_pbmc5.Rdata")

#------------------------------------------------------#
# CDSeq results with different MCMC and Reduce-recovery
#------------------------------------------------------#

# load the data
allfile_mc <- list.files(path="pbmc_result/", full.names = TRUE)[c(6,5,2,3,4)]
alldata_mc <- list()
for(i in 1:5){alldata_mc[[i]] <- load_obj(allfile_mc[i])}

# calcuate the correlation and RMSE
allcor_rmse_mc <- list()
for(jj in 1:5){
  
  allCDSeq <- alldata_mc[[jj]]
  rmse_cor <- matrix(NA, nrow = 10, ncol = 2, dimnames = list(1:10, c("RMSE","CORR")))
  
  for(i in 1:10){
    
    cdseqresult <- allCDSeq[[i]]$estProp
    trueth      <- allP[[i]] %>% as.matrix()
    cdseqresult <- cdseqresult[rownames(trueth),]
    
    eta <- colSums(C)[names(cellnumber)]
    cdseqresult <- CDSeq:::Cell2RNA(eta, cdseqresult)

    mycor  <- c()
    myrmse <- c()
    for(j in 1:nrow(trueth)){
      mycor[j]  <- cor(cdseqresult[j,],trueth[j,]) 
      myrmse[j] <- sqrt(mean((cdseqresult[j,]-trueth[j,])^2)) 
    }
    
    rmse_cor[i,1] <- mean(myrmse) %>% round(3)
    rmse_cor[i,2] <- mean(mycor) %>% round(3)
  }
  allcor_rmse_mc[[jj]] <- rmse_cor
}

allcor_mcmc  <- sapply(allcor_rmse_mc, function(x){x[,2]})
allrmse_mcmc <- sapply(allcor_rmse_mc, function(x){x[,1]})

colnames(allcor_mcmc)  <- c("2000 (RR)",700,1000,2000,3000)
colnames(allrmse_mcmc) <- c("2000 (RR)",700,1000,2000,3000)

# Plot the correlation and RMSE
plot_mcmc <- allcor_mcmc %>% as.data.frame() %>% gather("MCMC", "Corr", 1:5) %>% 
  mutate(MCMC = factor(MCMC, levels = c("2000 (RR)",700,1000,2000,3000)))

# plot_mcmc <- allrmse_mcmc %>% as.data.frame() %>% gather("MCMC", "Corr", 1:5) %>% 
#   mutate(MCMC = factor(MCMC, levels = c("2000 (RR)",700,1000,2000,3000)))


pp <- ggplot(plot_mcmc, aes(x=MCMC,y=Corr, color = MCMC)) + 
  geom_boxplot(width=0.4, position = position_dodge(0.7),outlier.shape = NA) + 
  geom_jitter(aes(color = MCMC), position=position_jitter(width=.1, height=0)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  #geom_hline(yintercept = 1,color = "red", lty= "dashed")+
  ylim(c(0,0.15))+
  labs(y="Correlation",x="MCMC iteration")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
#ggsave(filename = "pbmc_mcmc_rmse.png", pp, height = 3, width = 4.2, dpi = 500)

#--------------------------------------#
# Comparison with other methods Results
#--------------------------------------#

# load the results from other methods
load("pbmc_result/allBulk16.Rdata")
load("pbmc_result/pbmc5_linseed.Rdata")

allbulk <- list()
for(i in 1:length(allbulk_result)){
  
  rmse_cor_rr_bulk <- matrix(NA, nrow = 10, ncol = 2, dimnames = list(1:10, c("RMSE","CORR")))
  for (j in 1:10) {
    
    mytemp <- allbulk_result[[i]][[j]] %>% group_by(CT) %>% 
      dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
                       Pearson=cor(observed_values,expected_values) %>% round(.,4)) 
    
    rmse_cor_rr_bulk[j,1] <- mean(mytemp$RMSE) %>% round(3)
    rmse_cor_rr_bulk[j,2] <- mean(mytemp$Pearson) %>% round(3)
    
  }
  allbulk[[i]] <- rmse_cor_rr_bulk
}

methodname <- c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elastic_net","lasso","ridge","EPIC","ssKL","dtangle",
                "MuSiC","BisqueRNA", "SCDC")

bulkcor <- sapply(allbulk, function(x){x[,2]})
colnames(bulkcor) <- methodname

bulkrmse <- sapply(allbulk, function(x){x[,1]})
colnames(bulkrmse) <- methodname

bulkcor  <- bulkcor[,-c(11:13)] %>% t()
bulkrmse <- bulkrmse[,-c(11:13)] %>% t


plotcor <- rbind(bulkcor,pbmc5_linseed_cor,allcor_mcmc[,4],allcor_mcmc[,1])
#plotcor <- rbind(bulkrmse,pbmc5_linseed_rmse,allrmse_mcmc[,4], allrmse_mcmc[,1])

rownames(plotcor)[nrow(plotcor)-2] <- "Linseed"
rownames(plotcor)[nrow(plotcor)-1] <- "CDSeq"
rownames(plotcor)[nrow(plotcor)] <- "CDSeqRR"

plotfinal <- plotcor %>% t() %>% as.data.frame() %>% gather("Method", "Corr", 1:nrow(plotcor))
plotfinal$Method <- factor(plotfinal$Method, levels = rownames(plotcor))

ct1 <- ggplot(plotfinal, aes(x=Method,y=Corr, color = Method)) + 
  geom_boxplot(width=0.4, position = position_dodge(0.7),outlier.shape = NA) + #avoid plotting outliers twice
  geom_jitter(aes(color = Method), position=position_jitter(width=.1, height=0)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  #geom_hline(yintercept = 1,color = "red", lty= "dashed")+
  labs(y="RMSE",x="Method")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none")

#ggsave(filename = "pbmc_compare_rmse.png", ct1, height = 4, width = 5.5, dpi = 500)
