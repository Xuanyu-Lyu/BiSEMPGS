# fit the univariate model with simulated bivariate data to check the bias of the univariate estimates

test <- readRDS("Analysis/Paper/UniModelBias/biasbothAMVT_MVN_summary_list_500.rds")

test[[1]]$Minus2LogLikelihood

# OK now create a loop to run the model fit for all the data files

conditionNames <- c("onlyAM", "onlyVT", "bothAMVT", "bothAMVT_StrongCrossTrait")
source("PaperScripts/08-uniModelBias_OpenMx.R")
if(FALSE){
    for(condition in 1:length(conditionNames)){
        data_path <- paste0("Data/Paper/UniModelBias/", conditionNames[condition], "/nfam32000")
        l_files <- list.files(data_path, pattern = ".txt")

        for (k in 1:(length(l_files))){
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
            param_names[2] <- gsub("e", "VE", param_names[2])
            
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
}



## the code below is for the MVN data part.

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
        #df <- df[df$VE11>0 & df$VE12>0 & df$VE22>0,]
    }
    return(df)
}

for (condition in 1:4){
    summary_list_path <- paste0("Analysis/Paper/uniModelBias/bias", conditionNames[condition] ,"_MVN_summary_list_500.rds")
    summary_list <- readRDS(summary_list_path)
    assign(paste0("df_", conditionNames[condition]), getDf(summary_list))
    # change the column names
    # only read the colnames from a specific file
    v_colnames <- colnames(read.table("Analysis/Paper/UniModelBias/onlyAM_uniSEMPGS_estimates.txt", header = TRUE))
    v_colnames[2] <- gsub("e", "VE", v_colnames[2])
    #print(v_colnames)
    
    # Get the data frame
    df_temp <- get(paste0("df_", conditionNames[condition]))
    
    # Modify the column names
    colnames(df_temp)[1:13] <- v_colnames
    
    # Assign the modified data frame back to the original variable
    assign(paste0("df_", conditionNames[condition]), df_temp)
}

#df_checke <- read.table("Analysis/Paper/UniModelBias/onlyAM_uniSEMPGS_estimates_teste.txt", header = TRUE)


# # read the df final
df_final <- read.table("Analysis/Paper/UniModelBias/onlyAM_uniSEMPGS_estimates.txt", header = TRUE)
colnames(df_final)[2] <- gsub("e", "VE", colnames(df_final)[2])
conditionNames <- c("onlyAM", "onlyVT", "bothAMVT", "bothAMVT_StrongCrossTrait")

# Create a dataframe to store the true values from iterative math
df_true <- data.frame(matrix(NA, nrow = 13, ncol = 4))
colnames(df_true) <- conditionNames
rownames(df_true) <- colnames(df_final)

df_true[,1] <- c(.15, .6, 0.007591632, 0.1081802, 0.2, 0.7745967, 0.07893050, 0.19182940, 0.02864593, 1.9156826, 0.1868788, 0.7170956, 0.30150144)
df_true[,2] <- c(.15, .6, 0.007155983, 0.1071635, 0.2, 0.7745967, 0.07959381, 0.21451903, 0.02769168, 1.9739026, 0.1855578, 0.7191225, 0.3094718)
df_true[,3] <- c(.15, .6, 0.006934182, 0.1001404, 0.2, 0.7745967, 0.09258325, 0.1763591, 0.02634246, 2.0595697, 0.1898732, 0.7274934, 0.3490547)
df_true[,4] <- c(.0 , .6, 0.0004472696, -0.003577789, .2, 0.7745967, 0.07801939, -0.3643391, 0.0003908214, 1.8130360, 0.1397945, 0.5185387, 0.2732524)
# ready to plot the results

# chat's version of the permutation test
bootstrap_test_median <- function(data, m0, n_boot = 10000, seed = NULL) {
  # data:   Numeric vector of observations
  # m0:     Hypothesized median under H0
  # n_boot: Number of bootstrap iterations
  # seed:   Optional seed for reproducibility
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  data <- as.numeric(data)
  n <- length(data)
  
  # 1. Observed median and observed test statistic
  observed_median <- median(data, na.rm = TRUE)
  observed_stat <- observed_median - m0  # difference from hypothesized median
  observed_abs_stat <- abs(observed_stat)
  
  # 2. Center the data to reflect H0: shift so that its median = m0
  #    (remove the observed median, then add the hypothesized median)
  centered_data <- data - observed_median + m0
  
  # 3. Bootstrap replicates
  count_extreme <- 0
  for (i in seq_len(n_boot)) {
    # Sample with replacement from the centered data
    b_sample <- sample(centered_data, size = n, replace = TRUE)
    b_median <- median(b_sample, na.rm = TRUE)
    
    # Test statistic for the bootstrap sample
    b_stat <- b_median - m0
    if (abs(b_stat) >= observed_abs_stat) {
      count_extreme <- count_extreme + 1
    }
  }
  
  # 4. Two-sided p-value
  p_value <- count_extreme / n_boot
  
  return(p_value)
}

# a function to get descriptive statistics for one parameter and return a vector
getDescriptive <- function(df, param, file_tv){
    med <- median(df[[param]], na.rm = TRUE)
    MAD <- mad(df[[param]], na.rm = TRUE)
    trueValue <- file_tv[names(file_tv) == param]
    #cat(param, "\t", med, "\t", MAD, "\t", trueValue, "\n")
    # significance test if median is significantly different from the true value
    p_value <- bootstrap_test_median( df[[param]], m0 = trueValue, n_boot = 10000)
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
# descriptive_df <- data.frame(matrix(ncol = 5, nrow = length(v_true_values)))
# colnames(descriptive_df) <- c("median", "MAD", "trueValue", "p_value", "proportion_systematic")

# # fill the data frame with the descriptive statistics
# for (i in 1:length(v_true_values)) {
#     param <- names(v_true_values)[i]
#     # get the true value from the vector
#     trueValue <- v_true_values[i]
#     # get the descriptive statistics for the parameter
#     descriptive_df[i,] <- getDescriptive(df, param, v_true_values)
# }
# rownames(descriptive_df) <- names(v_true_values)

# a function to generate decriptive statistics for all parameters
makeDescriptiveDf <- function(df_estimates, v_true_values){
    descriptive_df <- data.frame(matrix(ncol = 5, nrow = length(v_true_values)))
    colnames(descriptive_df) <- c("median", "MAD", "trueValue", "p_value", "proportion_systematic")
    
    # fill the data frame with the descriptive statistics
    for (i in 1:length(v_true_values)) {
        param <- names(v_true_values)[i]
        # get the true value from the vector
        trueValue <- v_true_values[i]
        # get the descriptive statistics for the parameter
        descriptive_df[i,] <- getDescriptive(df=df_estimates, param, v_true_values)
    }
    rownames(descriptive_df) <- names(v_true_values)
    return(descriptive_df)
}

# make plots for each parameter
library(ggplot2)
library(patchwork)

# create a function to plot the histogram
plot_histogram <- function(df, param, trueValue, median, MAD, p_value) {
      # Calculate the maximum density value
    #     max_density <- max(density(df[[param]], na.rm = TRUE)$y)
  
    #   # Set the y-axis upper limit with a 5% expansion
    #     upper_limit <- max_density * 1.05
    p <- ggplot(df, aes_string(x = param)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "#1f77b4", alpha = 0.9) +
        geom_vline(xintercept = trueValue, color = "#d62728", linetype = "dashed", size = 1) +
        geom_vline(xintercept = median, color = "#b5cf6b", linetype = "dashed", size = 1) +
        #geom_vline(xintercept = median + MAD, color = "orange", linetype = "dashed", size = 1) +
       # geom_vline(xintercept = median - MAD, color = "orange", linetype = "dashed", size = 1) +
        scale_x_continuous(labels = function(x) sprintf("%.2f", x)) + 
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(title = paste(param),
             x = NULL,
             y = NULL) + 
             theme_minimal() +
        theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) +
            annotate("text", x = Inf, y = Inf, label = paste("p = ", round(p_value, 3)), color = "black", 
                    hjust = 1.1, vjust = 1.5)
    return(p)
}

# create a function to combine plots for specific parameters
combine_plots <- function(df, params, descriptive_df, ncol=4) {
    plots <- list()
    descriptive_df <- descriptive_df[]
    for (param in params) {
        trueValue <- descriptive_df[param,"trueValue"]
        median <- descriptive_df[param,"median"]
        MAD <- descriptive_df[param,"MAD"]
        p_value <- descriptive_df[param,"p_value"]
        plots[[param]] <- plot_histogram(df, param, trueValue, median, MAD, p_value)
    }
    combined_plot <- wrap_plots(plots, ncol = ncol)
}

# get plots for the first condition
df1 <- df_onlyAM[,1:13]
v_true_values1 <- df_true[,1]
names(v_true_values1) <- colnames(df1)
descriptive_df1 <- makeDescriptiveDf(df1, v_true_values1)
# plot them
p1 <- combine_plots(df1, names(v_true_values1)[-2], descriptive_df1, ncol = 3)
p1
ggsave(paste0("Analysis/Paper/", "UniBias_onlyAM_figure.png") , p1, width = 10, height = 12, type = "cairo-png", dpi = 600)


# get plots for the second condition
df2 <- df_onlyVT[,1:13]
#colnames(df2)[2] <- gsub("e", "VE", colnames(df2)[2])
v_true_values2 <- df_true[,2]
names(v_true_values2) <- colnames(df2)
descriptive_df2 <- makeDescriptiveDf(df2, v_true_values2)
# plot them
p2 <- combine_plots(df2, names(v_true_values2)[-2], descriptive_df2, ncol = 3)
p2
ggsave(paste0("Analysis/Paper/", "UniBias_onlyVT_figure.png") , p2, width = 10, height = 12, type = "cairo-png", dpi = 600)


# get plots for the third condition
df3 <- df_bothAMVT[,1:13]
#colnames(df3)[2] <- gsub("e", "VE", colnames(df3)[2])
v_true_values3 <- df_true[,3]
names(v_true_values3) <- colnames(df3)
descriptive_df3 <- makeDescriptiveDf(df3, v_true_values3)
# plot them
p3 <- combine_plots(df3, names(v_true_values3)[-2], descriptive_df3, ncol = 3)
p3
ggsave(paste0("Analysis/Paper/", "UniBias_bothAMVT_figure.png") , p3, width = 10, height = 12, type = "cairo-png", dpi = 600)

# get plots for the fourth condition
df4 <- df_bothAMVT_StrongCrossTrait[,1:13]
#colnames(df4)[2] <- gsub("e", "VE", colnames(df4)[2])
v_true_values4 <- df_true[,4]
names(v_true_values4) <- colnames(df4)
descriptive_df4 <- makeDescriptiveDf(df4, v_true_values4)
# plot them
p4 <- combine_plots(df4, names(v_true_values4)[-2], descriptive_df4, ncol = 3)
p4
ggsave(paste0("Analysis/Paper/", "UniBias_bothAMVT_StrongCrossTrait_figure.png") , p4, width = 10, height = 12, type = "cairo-png", dpi = 600)