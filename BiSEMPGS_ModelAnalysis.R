# This script is designed for looking into the output estimates and try to find a better way of fitting the BiSEMPGS model

# read the summary list
summary_list <- readRDS("Summary/Full_Model/m2_16000_summary_list.rds")
summary_list[[3]]$parameters
# extract all the estimates in the list and put each parameter as a column in a data frame
# Initialize an empty 78 column data frame
df <- data.frame(matrix(ncol = nrow(summary_list[[1]]$parameters), nrow = length(summary_list)))
colnames(df) <- summary_list$loop1.rds_16000.txt$parameters$name

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

plot(df$VY11, ylim = c(0,5))
abline(h = 1.7292875, col = "red", lwd = 2)
summary(df$VY11)


plot(df$f11, ylim = c(0,1))
abline(h = 0.15, col = "red", lwd = 2)
# Now df is a data frame where each column is the estimates from each element in the summary_list