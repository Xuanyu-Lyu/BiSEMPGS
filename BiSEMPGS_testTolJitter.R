# This is a script to find the best algorithm for the BiSEMPGS model. 
# Author: Xuanyu Lyu
# Date: 09/26/2024


# prepare the functions
source("BiSEMPGS_fitOpenMx.R")
library(crayon)

# the path to the dir with all txt files
# now we are using the 64k data for this test
data_path <- c("Data/Full_Model")
data_pattern <- c( "_48000.txt", "_32000.txt", "_64000.txt")
save_pattern <- c("_48000", "_32000", "_64000")
save_path <- "Analysis/FindBestTol_Full_Model"

# get all the names of the text files that end with 64000.txt
l_files <- list.files(data_path[1], pattern = data_pattern[3])

# feasibility tolerance vector and Optimality tolerance vectors
feaTol <- 10^seq(-3, -7, by = -1)
optTol <- 10^seq(-3, -7, by = -1)

feaTol <- 1e-4
optTol <- 1e-6


# run the simulation and save the list
for(i in 1:length(feaTol)){
    for(j in 1:length(optTol)){
        summary_list <- list()
        for (k in 1:10){
            fit <- fitBiSEMPGS_m2_tol(paste0(data_path[1], "/", l_files[k]), feaTol = feaTol[i], optTol = optTol[j], exhaustive = TRUE, extraTries = 10)
            summary_list[[l_files[k]]] <- fit
            cat(magenta("\nModel", l_files[k], "has been fitted with feaTol = ", feaTol[i], "and optTol = ", optTol[j], "\n"))
        }
        saveRDS(summary_list, paste0(save_path, "/m2_checka2_", "feaTol=", feaTol[i], "_optTol=", optTol[j], "_64000_summary_list.rds"))
        cat(green("\nSummary list has been saved with feaTol = ", feaTol[i], "and optTol = ", optTol[j], "\n"))
    }
}

# # try different values of distribution mean and variance
# save_path <- "Analysis/FindBestJitter_Full_Model"

# # vectors for jitter mean and variance

# v_jitterMean <- c(seq(0,3, by = .5))
# v_jitterVar <- c(seq(0,1, by = .2))

# # run the simulation and save the list

# for(i in 1:length(v_jitterMean)){
#     for(j in 1:length(v_jitterVar)){
#         summary_list <- list()
#         for (k in 1:100){
#             fit <- fitBiSEMPGS_m2_tol(paste0(data_path[1], "/", l_files[k]), jitterMean = v_jitterMean[i], jitterVar = v_jitterVar[j], extraTries = 5, exhaustive = TRUE)
#             summary_list[[l_files[k]]] <- fit
#             cat(magenta("\nModel", l_files[k], "has been fitted with jitterMean = ", v_jitterMean[i], "and jitterVar = ", v_jitterVar[j], "\n"))
#         }
#         saveRDS(summary_list, paste0(save_path, "/m2_", "jitterMean=", v_jitterMean[i], "_jitterVar=", v_jitterVar[j], "_64000_summary_list.rds"))
#         cat(green("\nSummary list has been saved with jitterMean = ", v_jitterMean[i], "and jitterVar = ", v_jitterVar[j], "\n"))
#     }
# }
