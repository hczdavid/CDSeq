
#--------------------------------------------#
# Deconvolution Benchmark Pipeline 
# Use PBMC single cell data to simulate data
# Author: David Huang
#--------------------------------------------#

# Note: this code is following the data generating pipeline from 
# https://github.com/favilaco/deconv_benchmark

library(tidyverse)
library(Matrix)

path <- "/Users/huangc9/Documents/CDSeq/new code/"
setwd(path)
load("/Users/huangc9/Documents/CDSeq/raw data/Broad_pbmc_scRNAseq.RData")
source("helper_functions1.R")

# get the raw data
data <- broad_pbmc_read_count_comm_cell
rm(broad_pbmc_read_count_comm_cell)

# get the annotation
pbmc_annotation$SubjectName <- sapply(strsplit(pbmc_annotation[,1],split = "_"), function(x){x[1]})
rownames(pbmc_annotation)   <- pbmc_annotation[,1]

#-------------------#
# Quality Control
#-------------------#

# First: cells with library size further than three MAD away were discarded
filterCells <- function(filterParam){
  cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
  cellsToRemove
}

libSizes      <- colSums(data)
gene_names    <- rownames(data)
cellsToRemove <- filterCells(libSizes)

if(length(cellsToRemove) != 0){
  data <- data[,-cellsToRemove]
  pbmc_annotation <- pbmc_annotation[-cellsToRemove,]
}

# Keep only "detectable" genes: at least 5% of cells (regardless of the group) 
# have a read/UMI count different from 0
keep  <- which(Matrix::rowSums(data > 0) >= round(0.05 * ncol(data)))
data  <- data[keep,]

original_cell_names <- colnames(data)
sum(colnames(data) != pbmc_annotation$cell_id) # cehck is match
colnames(data)      <- pbmc_annotation$cell_type

# Keep major cell types 
cell_counts  <-  table(pbmc_annotation$cell_type)
to_keep      <-  c("Megakaryocyte" , "CD14+ monocyte",  "Dendritic cell", "B cell", "Natural killer cell")
pData        <-  pbmc_annotation[pbmc_annotation$cell_type %in% to_keep,]
data1        <-  data[, colnames(data)%in% to_keep]
original_cell_names <- original_cell_names[colnames(data)%in% to_keep]


#-------------------#
# Creat T and C
#-------------------#

# Data split into train & test  
set.seed(123)
training <- as.numeric(unlist(sapply(unique(colnames(data1)), function(x){
  sample(which(colnames(data1) %in% x), cell_counts[x]/2)})))
testing <- which(!1:ncol(data1) %in% training)

# Generate phenodata for reference matrix C
pDataC  <- pData[training,]
pDataT  <- pData[testing,]
train   <- data1[,training]
test    <- data1[,testing]
cellGEP <- mergecol(test)
cellnumber <- table(colnames(test))

# get the train with original cell name
train_cellID <- train
colnames(train_cellID) <- original_cell_names[training]

# reference matrix (C) + refProfiles.var from training dataset
cellType <- colnames(train)
group    <- list()
for(i in unique(cellType)){ 
  group[[i]] <- which(cellType %in% i)
}

C  <-  lapply(group,function(x) Matrix::rowSums(train[,x])) 
C  <-  round(do.call(cbind.data.frame, C))

refProfiles.var <- lapply(group,function(x) train[,x])
refProfiles.var <- lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
refProfiles.var <- round(do.call(cbind.data.frame, refProfiles.var))
rownames(refProfiles.var) <- rownames(train)

# Generation of 10 datasets with each 100 pseudo-bulk mixtures (T)
cellType         <- colnames(test)
colnames(test)   <- original_cell_names[testing]
colnames(pDataT) <- c("cellID", "cellType", "SubjectName")

allT <- list()
allP <- list()

for(i in 1:10){
  generator <- Generator(sce = test, phenoData = pDataT, Num.mixtures = 100, pool.size = 500, seed = i*100)
  allT[[i]] <- generator[["T"]]
  allP[[i]] <- generator[["P"]]
}

save(allT, allP, cellGEP, C, refProfiles.var, train_cellID,pDataC, cellnumber,file = "readyData_pbmc5.Rdata")
