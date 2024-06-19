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
# no constraints on gc and hc, 16000 samples, lb = -.05, smaller tolerance, adjust f starting values
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_16000_summary_list.rds")
# no constraints on gc and hc, 32000 samples, lb = -.05, smaller tolerance, adjust f starting values
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_32000_summary_list.rds")
# no constraints on gc and hc, 48000 samples, lb = -.05, smaller tolerance, new setup for constraints
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_48000_summary_list.rds")
# no constraints on gc and hc, 64000 samples, lb = -.05, smaller tolerance, new setup for constraints
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_64000_summary_list.rds")

# no constraints on gc and hc, 16000 samples, lb = -.05, smaller tolerance, new setup, sqrt wvconstraint
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_16000_summary_list (1).rds")
# no constraints on gc and hc, 16000 samples, lb = -.05, smaller tolerance, new setup, no wvconstraint
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_freewv_16000_summary_list.rds")



# extract all the status code of openmx and put them into a vector
status_codes <- sapply(summary_list, function(x) x$statusCode)
summary(status_codes)
# extract all the estimates in the list and put each parameter as a column in a data frame
# Initialize an empty 78 column data frame
df <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
colnames(df) <- summary_list[[2]]$parameters$name
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
df <- df[-1,]
library(ggplot2)
library(tidyr)
# a plot for three VY estimates
true_values <- c(VY11 = 1.7292875, VY12 = 0.3693813,  VY22 = 1.1455869)
df_long <- tidyr::pivot_longer(df, c("VY11", "VY12", "VY22"), names_to = "Variable", values_to = "Value")
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

# a plot for four Omega estimates
true_values <- c(Omega11 = 0.42719846, Omega12 = 0.1124501, Omega21 = 0.07169126, Omega22 = 0.1497856)
df_long <- tidyr::pivot_longer(df, c("Omega11", "Omega12",  "Omega21", "Omega22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four Gamma estimates
true_values <- c(Gamma11 = 0.42874770, Gamma12 = 0.09775105, Gamma21 = 0.07171512, Gamma22 = 0.21521151)
df_long <- tidyr::pivot_longer(df, c("Gamma11", "Gamma12",  "Gamma21", "Gamma22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for three gc estimates
true_values <- c(gc11 = 0.04243337, gc12 = 0.009898546, gc22 = 0.006565526)
df_long <- tidyr::pivot_longer(df, c("gc11", "gc12",  "gc22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  ylim(0,.05)
  theme_minimal()

# a plot for three hc estimates
true_values <- c(hc11 = 0.04322826, hc12 = 0.01270896, hc22 = 0.01348196)
df_long <- tidyr::pivot_longer(df, c("hc11", "hc12",  "hc22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  ylim(0,.05)+
  theme_minimal()


# plot(df$VY11, ylim = c(0,5))
# abline(h = 1.7292875, col = "red", lwd = 2)
# summary(df$VY11)


# plot(df$f22, ylim = c(0,1))
# abline(h = 0.10, col = "red", lwd = 2)
# # Now df is a data frame where each column is the estimates from each element in the summary_list




