#!/usr/bin/env Rscript
#
#
## DISTRIBUTIONS COMPARATION 
### -> Calculate the difference and the interference (if this difference is real or random)

################################ CONFIGURATION ############################

#### Algorithm and model that we are going to compare
model_name <- 'iPAE1146_PAO1'
#algorithm <- 'mtf'
algorithm <- 'lmoma'
#algorithm <- 'room'
outdir <- '../flux_comparison_results/'

#### Do you want to normalize the fluxes? YES (1) NO (0)
normalize_fluxes <- 0

#### Do you want to use the separated reactions? YES (1) NO (0)
sep_reactions <- 0

###########################################################################
############### YOU SHOULDN'T NEED TO TOUCH NOTHING BELOW #################

####################  SEPARATED REACTIONS #################################
if (sep_reactions == 1) {
  sep_dir <- "/Separate_reactions/"
  sep_name <- "SepReact_"
  row_names <- 0
  config_SR <- "Separated Reactions"
} else {
  sep_dir <- ""
  sep_name <- ""
  config_SR <- "NO Separated Reactions"
}
##########################################################################

#################### NORMALIZE ###########################################
if (normalize_fluxes == 1) {
  norm_name <- "Norm_"
  config_N <- "Normalized Fluxes"
} else {
  norm_name <- ""
  config_N <- "NO Normalized Fluxes"
}
##########################################################################

##### File of WT MTF
mtfdir <- paste0('../mtf', sep_dir)
name_file_mtf <- grep("_MTF.", list.files(mtfdir), value = TRUE)
file_wt_mtf <- paste(mtfdir, name_file_mtf, sep='/')

#####################################################################
library(tidyr)
library(readr)
library(BSDA)
library(BBmisc)

## Importation of MTF data of the WT
wt_mtf <- read.table(file_wt_mtf)
if (dim(wt_mtf)[2] == 3) {
  wt_mtf <- wt_mtf[,2:3]
}
colnames(wt_mtf) <- c('Reaction', 'Flux')


## Creates a folder, a list and a table to hold results
if (dir.exists(outdir) == FALSE) {
  dir.create(outdir)
}

## Creates a folder to save normalized fluxes
dir_norm <- paste0(outdir, "Normalized_Fluxes_", sep_name)
if (normalize_fluxes == 1 & dir.exists(dir_norm) == FALSE) {
  dir.create(dir_norm)
}
results <- list()
  
  ### KOLMOGOROV-SMIRNOV TEST 
KS_test_results <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(KS_test_results) = c("Gene", "D", "p-value")

  ### Z-TEST 
Z_test_results <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(Z_test_results) = c("Gene", "Z", "p-value")

  ### WILCOXON TEST 
W_test_results <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(W_test_results) = c("Gene", "W", "p-value")

## The loops reads each file, filters all values selecting only those that differ 
## from the wild type and normalizes both datasets for their comparison
mut_dir <- paste0(".", sep_dir)
pattern_mut <- paste0(algorithm, ".txt$")
mut_files <- list.files(mut_dir,pattern= pattern_mut)
no_file <- 0

for (i in mut_files ) {
  no_file <- no_file + 1
  path <- paste0("./", sep_dir, i)
  # remove from start until "_ko_" and then from "_ till the end
  gene <- gsub('_lb.*', '', gsub('.*_ko_', '', path))
  
  mut <- read.table(file = path)
  no_mtf_reactions <- dim(wt_mtf)[1]
  no_mut_fluxes <- dim(mut)[1]
  if (no_mtf_reactions != no_mut_fluxes) {
    cat('ERROR: no. of reactions in mutant', path, 'different from MTF\n')
    quit()
  }
  colnames(mut) <- c('Reaction', 'Flux')
  
  # We compare only the fluxes that are different
  fluxes_wt <- wt_mtf$Flux[wt_mtf$Flux != mut$Flux]
  fluxes_mut <- mut$Flux[mut$Flux != wt_mtf$Flux]
  
  ############################### NORMALIZE ##############################################
  if (normalize_fluxes == 1){
    # Normalize
    fluxes_wt <- normalize(fluxes_wt)
    fluxes_mut <- normalize(fluxes_mut)
    
    ## Data normalized of mutants is assigned to a variable, added to a list and written to a file 
    flx_tmp <- mut[wt_mtf$Flux != mut$Flux, ]
    flx_tmp$Flux_Normalized <- fluxes_mut
    dynamicVariableName <- paste0(gene, '_norm')
    assign(dynamicVariableName, flx_tmp)
    results[[dynamicVariableName]] <- flx_tmp
    filename <- paste0(dir_norm, dynamicVariableName,
                       '.tsv')
    write.table(flx_tmp, file = filename, row.names = FALSE, sep = "\t")
  } 
  ########################################################################################
  
  cat('\nComparing wild type vs', gene, 'knockout:         ', no_file, "/",length(mut_files))
  
####### KOLMOGOROV-SMIRNOV TEST 
  KS_test <- ks.test(fluxes_wt, fluxes_mut)
  
  ## Save the results of the function in a table
  KS_test_results[no_file,] <- c(gene, KS_test$statistic[[1]], KS_test$p.value[[1]])
  
  ## Print the results
  #print(KS_test)
  cat("\nKS test done")

####### Z-TEST 
  sigma_wt <- sd(fluxes_wt)
  sigma_mut <- sd(fluxes_mut)
  Z_test <- z.test(fluxes_wt, fluxes_mut, sigma.x = sigma_wt, sigma.y = sigma_mut)
  
  ## Save the results of the function in a table
  Z_test_results[no_file,] <- c(gene, Z_test$statistic[[1]], Z_test$p.value)
  
  ## Print the results
  #print(Z_test)
  cat("\nZ-test done")

####### WILCOXON TEST 
  W_test <- wilcox.test(fluxes_wt, fluxes_mut)
    
  ## Save the results of the function in a table
  W_test_results[no_file,] <- c(gene, W_test$statistic[[1]], W_test$p.value)
  
  ## Print the results
  #print(Z_test)
  cat("\nW-test done")

}

## Save the results in a document

#### KOLMOGOROV-SMIRNOV TEST
name_file_KS <- paste0("/KStest_", sep_name, norm_name)
file_KStest <- paste0(outdir, name_file_KS, algorithm, '.csv')
write.csv(KS_test_results, file_KStest)

#### Z-TEST
name_file_Z <- paste0("/Ztest_", sep_name, norm_name)
file_Ztest <- paste0(outdir, name_file_Z, algorithm, '.csv')
write.csv(Z_test_results, file_Ztest)

#### WILCOXON TEST
name_file_W <- paste0("/Wtest_", sep_name, norm_name)
file_Wtest <- paste0(outdir, name_file_W, algorithm, '.csv')
write.csv(W_test_results, file_Wtest)

#### ANALYSIS OF RESULTS - We show the info and save it into a file (Summary_Significants_mutants.txt)
sig_mut_Z <- sum(Z_test_results$`p-value` < 0.05)
sig_mut_KS <- sum(KS_test_results$`p-value` < 0.05)
sig_mut_W <- sum(W_test_results$`p-value` < 0.05)


info <- c("\n\nMODEL:", model_name, "-- ALGORITHM:", algorithm,
"\n\nThe number of mutations with a p-value < 0.05 are:
 Z-test ->", sig_mut_Z,
"\n KS-test ->", sig_mut_KS,
"\n W-test ->", sig_mut_W,
"\n\nCONFIGURATION:", config_SR, "--", config_N)

cat(info, file = "../Summary_Significants_mutants.txt", append = TRUE)
cat(info)
