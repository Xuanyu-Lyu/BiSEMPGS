# This script is designed for looking into the output estimates and try to find a better way of fitting the BiSEMPGS model

# read the summary list
# no constraints on gc and hc, 16000 samples, lb = -.05
summary_list <- readRDS("Analysis/Full_Model/m2_16000_summary_list.rds")
# constraints on everything, 16000 samples, lb = -.05
summary_list <- readRDS("Analysis/Full_Model/m2_allConst_16000_summary_list.rds")
# no constraints on gc and hc, 16000 samples, lb = -.55
summary_list <- readRDS("Analysis/Full_Model/m2_.55lb_16000_summary_list.rds")
# no constraints on gc and hc, 48000 samples, lb = -.55, smaller tolerance, adjust f starting values
summary_list <- readRDS("Analysis/Full_Model/m2_.55lb_smallerTol_adjustf_48000_summary_list.rds")
# no constraints on gc and hc, 48000 samples, lb = -.05, smaller tolerance, new setup for constraints
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_48000_summary_list.rds")

# extract all the status code of openmx and put them into a vector
status_codes <- sapply(summary_list, function(x) x$statusCode)
summary(status_codes)
# extract all the estimates in the list and put each parameter as a column in a data frame
# Initialize an empty 78 column data frame
df <- data.frame(matrix(ncol = nrow(summary_list[[1]]$parameters), nrow = length(summary_list)))
colnames(df) <- summary_list[[1]]$parameters$name
colnames(df) 
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
aggregate(df$f11, by = list(status_codes), FUN = mean)

library(ggplot2)

# a plot for three VY estimates
true_values <- c(VY11 = 1.7292875, VY12 = 0.3693813,  VY22 = 1.1455869)
df_long <- tidyr::pivot_longer(df, c("VY11", "VY12",  "VY22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four f estimates
true_values <- c(f11 = 0.15, f12 = 0.1, f21 = 0.05, f22 = 0.1)
df_long <- tidyr::pivot_longer(df, c("f11", "f12",  "f21", "f22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for two delta estimates
true_values <- c(delta11 = sqrt(.49*.5), delta22 = sqrt(.16*.3))
df_long <- tidyr::pivot_longer(df, c("delta11", "delta22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for two a estimates
true_values <- c(a11 = sqrt(.49*.5), a22 = sqrt(.16*.7))
df_long <- tidyr::pivot_longer(df, c("a11", "a22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four v estimates
true_values <- c(v11 = 0.20214807, v12 = 0.09995314, v21 = 0.08287606, v22 = 0.07056632)
df_long <- tidyr::pivot_longer(df, c("v11", "v12",  "v21", "v22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four w estimates
true_values <- c(w11 = 0.20148136, w12 = 0.08849580, w21 = 0.08263773, w22 = 0.05576877)
df_long <- tidyr::pivot_longer(df, c("w11", "w12",  "w21", "w22"), names_to = "Variable", values_to = "Value")
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