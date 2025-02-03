# This is a script to make all the sumplementary tables for the Multivariate SEM-PGS paper
conditionNames <- c("Model_r2_1", "Model_r2_2", "Model_r2_4", "Model_r2_8", "Model_r2_16")
sample_sizes <- c(4000, 8000, 16000, 32000, 48000, 64000)

# the function to get the df
getDf <- function(summary_list, fixed = FALSE) {
    status_codes <- sapply(summary_list, function(x) x$statusCode)
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
    df <- df[df$status_codes %in% c("OK", "OK/green"),]
    # exclude a that hit the lower bound
    if(!fixed){
        #df <- df[df$a11!=0.3 & df$a22!=0.3,]
        #df <- df[df$a11>0.4 & df$a22>0.4,]
        #df <- df[df$VY11>1 & df$VY22>1,]
       # df <- df[df$VE11>0 & df$VE12>0 & df$VE22>0,]
    }
  
    # exclude the outliers that is three standard deviations away from the mean, only applied to numerical variables
    #df <- df[apply(df[,1:(ncol(df)-15)], 2, function(x) all(abs(x - mean(x, na.rm = TRUE)) < 8*sd(x, na.rm = TRUE))),]

    # general lower bound variables
    varname_vector <- c()
    return(df)
}
# a function to get descriptive statistics for one parameter and return a vector
getDescriptive <- function(df, param, file_tv){
    med <- median(df[[param]], na.rm = TRUE)
    MAD <- mad(df[[param]], na.rm = TRUE)
    trueValue <- file_tv[file_tv$V1 == param,]$V2
    #cat(param, "\t", med, "\t", MAD, "\t", trueValue, "\n")
    # significance test if median is significantly different from the true value
    p_value <- wilcox.test(df[[param]], mu = trueValue, alternative = "two.sided")$p.value

    # compute the total variance and the variance from randomness
    var_total <- var(df[[param]], na.rm = TRUE)
    var_random <- sum((df[[param]] - trueValue)^2)/(length(df[[param]]) - 1)
    var_systematic <- var_random - var_total
    p_systematic <- var_systematic/var_random

    final_v <- c(med, MAD, trueValue, p_value, p_systematic)
    
    names(final_v) <- c("median", "MAD", "trueValue", "p_value", "proportion_systematic")
    return(final_v)

}


all_df <- list()
all_des <- list()
for (i in 1:length(conditionNames)){
    file_tv <- read.table(paste0("Data/Paper/Expected/",conditionNames[i],"_finalGen.txt"))
    df_list <- list()
    des_list <- list()
    for (j in 1:length(sample_sizes)){
        data_path <- paste0("Analysis/Paper/", conditionNames[i], "/m2_paper_version2_", sample_sizes[j],"_summary_list.rds")

        df_estimates <- getDf(readRDS(data_path))
        df_estimates <- df_estimates[,colnames(df_estimates) %in% file_tv$V1]
        #print(colnames(df_estimates))
        # get descriptive statistics for each parameter and save in a dataframe
        df_estimates_des <- data.frame(matrix(ncol = 6, nrow = ncol(df_estimates)))
        colnames(df_estimates_des) <- c("param","median", "MAD", "trueValue", "p_value", "proportion_systematic")
        for (k in 1:ncol(df_estimates)){
            df_estimates_des[k,1] <- colnames(df_estimates)[k]
            df_estimates_des[k,2:6] <- getDescriptive(df_estimates, colnames(df_estimates)[k], file_tv = file_tv)
        }

        df_list[[j]] <- df_estimates
        names(df_list)[j] <- paste0("n=", sample_sizes[j])
        des_list[[j]] <- df_estimates_des
        names(des_list)[j] <- paste0("n=", sample_sizes[j])
    }
    all_df[[i]] <- df_list
    names(all_df)[i] <- conditionNames[i]
    all_des[[i]] <- des_list
    names(all_des)[i] <- conditionNames[i]
    
}
