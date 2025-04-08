# fit the univariate model with simulated bivariate data to check the bias of the univariate estimates

# read the data to format it into the right format
# testData <- read.table("Data/Paper/UniModelBias/onlyAM/nfam32000/loop1.rds_32000.txt", header = TRUE)

# # subset the data to only include varibles whose names end with 1
# testData <- testData[, grepl("1$", names(testData))]
# head(testData)

# # arrange the order of columns to match the OpenMx script
# testData <- testData[, c("NTm1", "Tm1", "NTp1", "Tp1", "Ym1", "Yp1", "Yo1")]
# head(testData)

# # test the model fit function
# source("PaperScripts/08-uniModelBias_OpenMx.R")

# testFit <- fitUniSEMPGS(testData, max.cores = 2)

# param_names <- gsub("\\d+", "", testFit$param$matrix[1:13])
# param_names <- gsub("sigma", "VY", param_names)
# param_names <- gsub("E", "VE", param_names)

# df_estimate <-data.frame(matrix(NA, ncol = 13))
# colnames(df_estimate) <- param_names[1:13]
# head(df_estimate)

# df_estimate[1,] <- testFit$param$Estimate[1:13]

# head(df_estimate)


# OK now create a loop to run the model fit for all the data files

conditionNames <- c("onlyAM", "onlyVT", "bothAMVT")
source("PaperScripts/08-uniModelBias_OpenMx.R")
for(condition in 1:length(conditionNames)){
    data_path <- paste0("Data/Paper/UniModelBias/", conditionNames[condition], "/nfam32000")
    l_files <- list.files(data_path, pattern = ".txt")

    for (k in 1:length(l_files)){
        # load the data
        testData <- read.table(paste0(data_path, "/", l_files[k]), header = TRUE)
        # subset the data to only include varibles whose names end with 1
        testData <- testData[, grepl("1$", names(testData))]
        # arrange the order of columns to match the OpenMx script
        testData <- testData[, c("NTm1", "Tm1", "NTp1", "Tp1", "Ym1", "Yp1", "Yo1")]
        # fit the model
        testFit <- fitUniSEMPGS(testData, max.cores = 2)
        
        param_names <- gsub("\\d+", "", testFit$param$matrix[1:13])
        param_names <- gsub("sigma", "VY", param_names)
        param_names <- gsub("E", "VE", param_names)
        
        df_estimate <-data.frame(matrix(NA, ncol = 13))
        colnames(df_estimate) <- param_names[1:13]
        
        df_estimate[1,] <- testFit$param$Estimate[1:13]
        if(k == 1){
            df_final <- df_estimate
        } else {
            df_final <- rbind(df_final, df_estimate)
        }

    }
     # save the results
    write.table(df_final, file = paste0("Analysis/Paper/UniModelBias/", paste0(conditionNames[condition],"_uniSEMPGS_estimates.txt")), sep = "\t", row.names = FALSE, col.names = TRUE)

}
