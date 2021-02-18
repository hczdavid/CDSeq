#---------------#
# Run CDseq
#---------------#

library(CDSeq)
library(foreach)
library(doParallel)


load("readyData_pbmc5.Rdata")

# set mcmc = 700,  gene_subset_size = NULL, dilution_factor = 10
# set mcmc = 1000, gene_subset_size = NULL, dilution_factor = 10
# set mcmc = 2000, gene_subset_size = NULL, dilution_factor = 10
# set mcmc = 1000, gene_subset_size = 800, dilution_factor = 10, block_number = 10


registerDoParallel(detectCores()-1)
allCDSeq <- foreach(i = 1:10, .inorder = FALSE)%dopar%{
  
  result <- CDSeq(bulk_data =  allT[[i]], cell_type_number=5, mcmc_iterations = 2000, 
                gene_subset_size = NULL, block_number = 1,
                dilution_factor = 10, gene_length = NULL, 
                reference_gep = C)
  return(result)
}
stopImplicitCluster() 
save(allCDSeq, file = "allCDSeq5_noRR_2000_10.Rdata")





