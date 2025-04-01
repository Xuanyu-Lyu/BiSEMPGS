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
# a functon to use permutation to get the p values
permutationTest <- function(X, null_median, n_iter=10000) {
  # Step 1: Compute the observed median
  m_obs <- median(X)
  # Step 2: Center the data so the median becomes the null median
  X_centered <- X - m_obs + null_median
  # Step 3: Generate the null distribution of the median using bootstrap resampling
  boot_medians <- replicate(n_iter, {
    sample_data <- sample(X_centered, size = length(X_centered), replace = TRUE)
    median(sample_data)
  })
  # Step 4: Calculate the two-tailed p-value
  p_lower <- mean(boot_medians <= m_obs)
  p_upper <- mean(boot_medians >= m_obs)
  p_value <- 2 * min(p_lower, p_upper)
  #p_value <- min(p_value, 1)  # Ensure p-value does not exceed 1
  # Return the results as a list
#   return(list(
#     observed_median = m_obs,
#     p_value = p_value
#     #,boot_medians = boot_medians
#   ))
    # Return the p value
    return(p_value)
}

# a function to get descriptive statistics for one parameter and return a vector
getDescriptive <- function(df, param, file_tv){
    med <- median(df[[param]], na.rm = TRUE)
    MAD <- mad(df[[param]], na.rm = TRUE)
    trueValue <- file_tv[file_tv$V1 == param,]$V2
    #cat(param, "\t", med, "\t", MAD, "\t", trueValue, "\n")
    # significance test if median is significantly different from the true value
    p_value <- permutationTest( df[[param]], null_median = trueValue, n_iter = 10000)
    #p_value <- wilcox.test(df[[param]], mu = trueValue, alternative = "two.sided")$p.value

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
        data_path <- paste0("Analysis/Paper/", conditionNames[i], "/m2_paper_version2_", sample_sizes[j],"_summary_list_fixedA.rds")

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

# create a seperate df for each condition
for (i in 1:length(conditionNames)){
    for (j in 1:length(sample_sizes)){
        all_des[[i]][[j]]$sample_sizes <- sample_sizes[j]
        # move the sample_sizes column to the first column
        all_des[[i]][[j]] <- all_des[[i]][[j]][,c(ncol(all_des[[i]][[j]]), 1:(ncol(all_des[[i]][[j]])-1))]
        if(j == 1){
            des <- all_des[[i]][[j]]
        } else {
           des <- rbind(des, all_des[[i]][[j]])
        }
    } 
    # order the rows by param and sample sizes
    des <- des[order(des$param, des$sample_sizes),]
    des <- des[, c("param","sample_sizes", "median", "MAD", "trueValue", "p_value", "proportion_systematic")]
    assign(paste0("des_", conditionNames[i]), des)

}

# get the latex code for each of the tables
library(kableExtra)
library(dplyr)

getLatex <- function(df){
 # Store the param column separately
  param_col <- df$param
  
  # Remove the param column from the data frame
  df <- df[, !names(df) %in% "param"]
  # Round numeric columns to four digits
  df <- df %>% mutate(across(where(is.numeric), round, 4))
  
  # Generate the LaTeX table without the param column
  kable(df, format = "latex", booktabs = TRUE, row.names = FALSE, longtable = TRUE) %>%
    kable_styling(latex_options = c("striped", "hold_position")) %>%
    pack_rows(index = table(param_col))
}

for (i in 1:length(conditionNames)){
    assign(paste0("latex_des_", conditionNames[i]), getLatex(get(paste0("des_", conditionNames[i]))))
}

# save the five latex code to a text file
all_latex <- c("r2pgs = .16", latex_des_Model_r2_16,
               "r2pgs = .08", latex_des_Model_r2_8,
               "r2pgs = .04", latex_des_Model_r2_4,
               "r2pgs = .02", latex_des_Model_r2_2,
               "r2pgs = .01", latex_des_Model_r2_1)
write(all_latex, "Analysis/Paper/Figure on manu/Appendix/latex_des_fixedA.txt")
