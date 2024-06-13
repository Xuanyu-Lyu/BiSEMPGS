# This script is designed for looking into the output estimates and try to find a better way of fitting the BiSEMPGS model

# read the summary list
# no constraints on gc and hc, 16000 samples, lb = -.05
summary_list <- readRDS("Analysis/Full_Model/m2_16000_summary_list.rds")
# constraints on everything, 16000 samples, lb = -.05
summary_list <- readRDS("Analysis/Full_Model/m2_allConst_16000_summary_list.rds")
# no constraints on gc and hc, 16000 samples, lb = -.55
summary_list <- readRDS("Analysis/Full_Model/m2_.55lb_16000_summary_list.rds")

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

sum(df$VY11 <1.9, na.rm = TRUE)

library(ggplot2)

# a plot for three VY estimates
true_values <- c(VY11 = 1.7292875, VY12 = 0.3693813,  VY22 = 1.1455869)

df_long <- tidyr::pivot_longer(df, c("VY11", "VY12",  "VY22"), names_to = "Variable", values_to = "Value")
# Add an index variable
df_long$Index <- 1:nrow(df_long)

ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four f estimates
true_values <- c(f11 = 0.15, f12 = 0.1, f21 = 0.05, f22 = 0.1)

df_long <- tidyr::pivot_longer(df, c("f11", "f12",  "f21", "f22"), names_to = "Variable", values_to = "Value")
# Add an index variable
df_long$Index <- 1:nrow(df_long)

ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

plot(df$VY11, ylim = c(0,5))
abline(h = 1.7292875, col = "red", lwd = 2)
summary(df$VY11)


plot(df$f22, ylim = c(0,1))
abline(h = 0.10, col = "red", lwd = 2)
# Now df is a data frame where each column is the estimates from each element in the summary_list