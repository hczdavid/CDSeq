#--------------------------------#
# Run Other deconvolution methods
#--------------------------------#

library(CDSeq)
library(foreach)
library(doParallel)
library(tidyverse)


load("readyData_pbmc5.Rdata")
source("helper_functions1.R")
source("CIBERSORT.R")


bulk_methods <- c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elastic_net","lasso","ridge","EPIC","ssKL","dtangle",
                   "MuSiC","BisqueRNA", "SCDC")

colnames(pDataC)[2] <- "cellType"

# The effect of Transformation and Scaling is already evaluated in Cobos et al, 2020
# Here we implement one combination: no transformation bulk data and apply column scale
allT_ts <- allT
for(i in 1:10){
  allT_ts[[i]] <- Transformation(allT[[i]], "none")
  allT_ts[[i]] <- Scaling(allT_ts[[i]], "column")
}

C <- Transformation(C, "none")
C <- Scaling(C, "column")

allbulk_result <- list()
for(i in 1:length(bulk_methods)){
  registerDoParallel(detectCores()-1)
  tempresult <- foreach(j = 1:10, .inorder = FALSE)%dopar%{
    print(paste0("Current running method *", bulk_methods[i], "* at ", j, " iteration"))
    Result <- Deconvolution(T =allT_ts[[j]], C = C, phenoDataC = pDataC, method = bulk_methods[i], P = allP[[j]], refProfiles.var = refProfiles.var) 
    return(Result)
  }
  stopImplicitCluster() 
  allbulk_result[[i]] <- tempresult
}
save(allbulk_result, file = "allBulk16.Rdata")
