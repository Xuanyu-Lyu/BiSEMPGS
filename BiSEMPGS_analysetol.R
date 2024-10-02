# a script to analyse the model fitting results using different tolerance values

# load the necessary libraries
library(psych)

# specify the folder
folder <- "Analysis/FindBestTol_Full_Model/smallTol"

# get all the files in the folder
l_files <- list.files(folder)

l_df <- list()
for (files in l_files){
    summary_list <- readRDS(paste0(folder, "/", files))
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
    df$status_codes <- status_codes
    #aggregate(df$f11, by = list(df$status_codes), FUN = mean)
    #df <- df[-1,]
    # get only the results with green status code
    df <- df[df$status_codes %in% c("OK", "OK/green"),]
    nrow(df)
    l_df[[files]] <- df
}

describe(l_df[[1]])

# get the summary of the dataframe and save the medians in a data frame
est_sum_df <- data.frame(matrix(ncol = length(l_df), nrow = 77))
for (i in 1:length(l_df)){
    est_sum_df[,i] <- c(describe(l_df[[i]])$mean,describe(l_df[[i]])[1,2])

}
rownames(est_sum_df) <- c(colnames(l_df[[1]]),"n")
name_l_files <- gsub("_64000_summary_list.rds", "", l_files)
colnames(est_sum_df) <- name_l_files


# check the graphics
names(l_df)
df <- l_df[["m2_checka_feaTol=1e-06_optTol=1e-07_64000_summary_list.rds"]]
