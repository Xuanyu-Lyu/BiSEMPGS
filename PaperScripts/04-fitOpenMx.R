source("PaperScripts/04-OpenMxFunctions.R")
library(crayon)
conditionNames <- c("Model_latent30", "Model_latent50", "Model_latent70", "Model_latent90")

folder_pattern <- c("nfam16000", "nfam32000", "nfam48000", "nfam64000", "nfam80000")
data_pattern <- c("_16000.txt", "_32000.txt", "_48000.txt",  "_64000.txt", "_80000.txt")
save_pattern <- c("_16000", "_32000", "_48000", "_64000", "_80000")

# for (i in 1:length(conditionNames)){
#     for (j in 1:length(data_pattern)){
#         data_path <- paste0("Data/Paper/", conditionNames[i], "/", folder_pattern[j])
#         l_files <- list.files(data_path, pattern = data_pattern[j])
#         # check if the data has been fitted
#         # if (file.exists(paste0("Analysis/Paper/", conditionNames[i], "/m2_paper", save_pattern[j], "_summary_list.rds"))){
#         #     cat("Summary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been fitted\n")
#         #     next
#         # }
#         summary_list <- list()
#         for (k in 1:length(l_files)){
#             fit <- fitBiSEMPGS_m2_tol(paste0(data_path, "/", l_files[k]), 
#                                       feaTol = 1e06, 
#                                       optTol = 1e08,
#                                       jitterMean = 0.5,
#                                       jitterVar = .1,
#                                       exhaustive = FALSE,
#                                       extraTries = 25)
#             summary_list[[l_files[k]]] <- fit
#             cat(magenta("\n", conditionNames[i], "\tn=", folder_pattern[j], "\tModel", l_files[k], "has been fitted\n"))
#         }
#         save_path <- paste0("Analysis/Paper/", conditionNames[i], "/m2_paper", save_pattern[j], "_summary_list.rds")
#         saveRDS(summary_list, save_path)
#         cat(green("\nSummary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been saved\n"))

#     }
# }


true_a <- list(c(0.438178, 0.464758),
            c(0.5656854,0.464758),
            c(0.669328,0.464758),
            c(0.7589466,0.464758))
for (i in 1:2){
    for (j in 1:1){
        data_path <- paste0("Data/Paper/", conditionNames[i], "/", folder_pattern[j])
        l_files <- list.files(data_path, pattern = data_pattern[j])
        # check if the data has been fitted
        # if (file.exists(paste0("Analysis/Paper/", conditionNames[i], "/m2_paper", save_pattern[j], "_summary_list_fixedA.rds"))){
        #     cat("Summary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been fitted\n")
        #     next
        # }
        summary_list <- list()
        for (k in 1:length(l_files)){
            fit <- fitBiSEMPGS_m2_tol(paste0(data_path, "/", l_files[k]), 
                                      #avalue = true_a[[i]],
                                      feaTol = 1e05, 
                                      optTol = 1e09,
                                      jitterMean = 0.5,
                                      jitterVar = .1,
                                      exhaustive = FALSE,
                                      extraTries = 25)
            summary_list[[l_files[k]]] <- fit
            cat(magenta("\n", conditionNames[i], "\tn=", folder_pattern[j], "\tModel", l_files[k], "has been fitted\n"))
        }
        save_path <- paste0("Analysis/Paper/test/", conditionNames[i], "/m2_paper", save_pattern[j], "_summary_list_fixedA.rds")
        saveRDS(summary_list, save_path)
        cat(green("\nSummary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been saved\n"))

    }
}