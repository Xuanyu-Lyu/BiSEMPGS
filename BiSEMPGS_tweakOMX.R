# a script to locally test openmx fit setup for the bivariate SEMPGS model
# author: Xuanyu Lyu
# date: 09/17/2024

source("BiSEMPGS_fitOpenMx.R")
conditionNames <- c("Full_Model", "MeasurePgs30", "MeasurePgs10", "MeasurePgsFully", 
"f11-decrease", "f12-decrease", "f11.12.21.22-decrease", "am11-decrease", "am12-decrease", "am11.12.21.22-decrease", "Full_Model_.5latent")
# the path to the dir with all txt files
data_path <- c(paste0("Data/", conditionNames[1]))
# the path to the dir to save the results
save_path <- paste0("Data/", conditionNames[1], "/Local_Analysis")

data_pattern <- c("_64000.txt")
save_pattern <- c( "_64000")
model_type <- "m2"
mxSetup <- "-.05lb_freeArg_VF1e-4"
n_models <- "All"
for (j in 1){
    # a list to save all the summary data
    summary_list <- list()

    # get all the names of the text files that end with 32000.txt
    l_files <- list.files(data_path[j], pattern = data_pattern[j])
    #l_files <- l_files[1:n_models]

    # fit the model for each data file
    for (i in 1:length(l_files)){
        # fit the model
        #fit <- fitBiSEMPGS_m2(paste0(data_path[j], "/", l_files[i]))
        
        # fit the model with fixed a
        fit <- fitBiSEMPGS_m2(paste0(data_path[j], "/", l_files[i]))
        # save the summary
        summary_list[[l_files[i]]] <- fit
        # save the fit
        cat("Model", l_files[i], "has been fitted\n")
    }
    if (!dir.exists(save_path)){
        dir.create(save_path)}

    # save the summary list
    saveRDS(summary_list, paste0(save_path, "/", model_type,mxSetup, save_pattern[j], "_nModel", n_models, "_summary_list.rds"))
}




# analyse the results
summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_forceJ_64000_nModelAll_summary_list.rds") # try large scale/ bad vf 

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-10lb_freeArg_forceJ_64000_nModelAll_summary_list.rds") # bad a

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_forceJ_VFfree_64000_nModelAll_summary_list.rds") # bad vf11 and ve11

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_forceJ_VFnolb_64000_nModelAll_summary_list.rds") # bad vf hc

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_forceJ_VF-3_noVYconst_64000_nModelAll_summary_list.rds")  #bad vf ve a

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_forceJ_VF-.05_noVFconst_64000_nModelAll_summary_list.rds")  # bad a vf

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_fixArg_VFnolb_64000_nModelAll_summary_list.rds")  # try large scale/ seems to be the best, if excluding outlier VFs

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_VF1e-4_64000_nModelAll_summary_list.rds") #  very bad a ve vf

# extract all the status code of openmx and put them into a vector
status_codes <- sapply(summary_list, function(x) x$statusCode)
summary(status_codes)
# extract all the estimates in the list and put each parameter as a column in a data frame
# Initialize an empty 78 column data frame
df <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
colnames(df) <- summary_list[[2]]$parameters$name
#colnames(df) 
# Loop over the elements in the summary_list
for(i in 1:length(summary_list)) {
    for(j in 1:nrow(summary_list[[i]]$parameters)){
        if (!is.null(summary_list[[i]]$parameters$Estimate[j])) {
            df[i,j] <- summary_list[[i]]$parameters$Estimate[j]
        } else {
            print(paste("NULL value at i =", i, "and j =", j))
        }
    }
}

# # initialize a df for the standard errors\
# df_se <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
# colnames(df_se) <- summary_list[[2]]$parameters$name
# colnames(df_se)
# # Loop over the elements in the summary_list
# for(i in 1:length(summary_list)) {
#     for(j in 1:nrow(summary_list[[i]]$parameters)){
#         if (!is.null(summary_list[[i]]$parameters$Std.Error[j])) {
#             df_se[i,j] <- summary_list[[i]]$parameters$Std.Error[j]
#         } else {
#             print(paste("NULL value at i =", i, "and j =", j))
#         }
#     }
# }
# # show how many rows of the dataframe has NA
# sum(is.na(df_se))/prod(dim(df_se))



df$status_codes <- status_codes
#aggregate(df$f11, by = list(df$status_codes), FUN = mean)
df <- df[-1,]
# get only the results with green status code
df <- df[df$status_codes %in% c("OK", "OK/green"),]
nrow(df)
# get only the results with positive f11 and f22 estimates
#df <- df[df$f11 > 0 & df$f22 > 0,]
#aggregate(df$f11, by = list(df$status_codes), FUN = mean)

#remove the

#remove all rows with any value smaller than -.045 exactly
#df <- df[!apply(df[,1:64] <= -0.048, 1, any), ]

#remove the outliers that are three sd away from the mean of VF11 VF12 and VF22
df <- df[abs(df$VF11 - mean(df$VF11)) < 3*sd(df$VF11),]
nrow(df)


library(psych)
describe(df)
