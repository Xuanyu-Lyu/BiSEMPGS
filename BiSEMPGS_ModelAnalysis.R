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

# no constraints on gc and hc, 48k samples, lb = -.05, smaller tolerance, new setup, fixed a
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_fixedA_48000_summary_list.rds")
# no constraints on gc and hc, 32k samples, lb = -.05, smaller tolerance, new setup, fixed a
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_fixedA_32000_summary_list.rds")


# no constraints on gc and hc, 32k samples, lb = -.05, smaller tolerance, new setup, fixed a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_fixedArg_32000_summary_list.rds")

# no constraints on gc and hc, 64k samples, lb = -.05, smaller tolerance, new setup, fixed a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_fixedArg_64000_summary_list.rds")

# no constraints on gc and hc, 64k samples, lb = -.05, smaller tolerance, new setup, fixed a, rg, closer h
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_smallerTol_newSetup_fixedArg_closerh_64000_summary_list.rds")
#summary_list[[2]]$parameters$Std.Error[1]

# no constraints on gc and hc, 64k samples, lb = -.001 for f and hc, smaller tolerance, new setup, fixed a, rg, closer h
summary_list <- readRDS("Analysis/Full_Model/m2_.001lbfh_sTol_newSetup_fixedArg_closerh_64000_summary_list.rds")

# constraint on hc but not gc, 48k samples, lb = -0.5, fixed a, rg, 100 models
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_fixedArg_closerh_hcCon_48000_nModel100_summary_list.rds")

# constraint on hc but not gc, 32k samples, lb = -0.5, fixed a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_fixedArg_closerh_hcCon_32000_nModelAll_summary_list.rds")

# constraint on hc but not gc, 48k samples, lb = -0.5, fixed a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_fixedArg_closerh_hcCon_48000_nModelAll_summary_list.rds")

# constraint on hc but not gc, 64k samples, lb = -0.5, fixed a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_fixedArg_closerh_hcCon_64000_nModelAll_summary_list.rds")

# constraint on hc but not gc, 32k samples, lb = -0.5, free a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_32000_nModelAll_summary_list.rds")

# constraint on hc but not gc, 48k samples, lb = -0.5, free a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_48000_nModelAll_summary_list.rds")

# constraint on hc but not gc, 64k samples, lb = -0.5, free a, rg
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_64000_nModelAll_summary_list.rds")

# constraint on hc but not gc, 32k samples, lb = -0.5, free a, rg, constraint on j
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_forceJ_32000_nModelAll_summary_list.rds") # does not work

# constraint on hc but not gc, 48k samples, lb = -0.5, free a, rg, constraint on j
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_forceJ_48000_nModelAll_summary_list.rds") # does not work

# constraint on hc but not gc, 64k samples, lb = -0.5, free a, rg, constraint on j
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_forceJ_64000_nModelAll_summary_list.rds") 
# does not work well, a is underestimated. When using stringent filtering, the estimates are not that bad.

# constraint on hc but not gc, 64k samples, lb = -0.5, free a, rg, constraint on j, using matt's constraints
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_forceJ_tweakConst_48000_nModelAll_summary_list.rds") # does not work, a and delta are too off.

# constraint on hc but not gc, 64k samples, lb = -0.5, free a, rg, constraint on j, using my constraints + Omega Constraints
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_closerh_hcCon_forceJ_tweakConst2_48000_nModel100_summary_list.rds") # does not work, a and delta are too off.

# constraint on hc but not gc, 48k samples, lb = -10, free a, rg, constraint on j
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_forceJ_-10lb_48000_nModel100_summary_list.rds") 

# constraint on hc but not gc, 48k samples, lb = 1e4, free a, rg, constraint on j
summary_list <- readRDS("Analysis/Full_Model/m2_.05lb_freeArg_forceJ_1e4lb_48000_nModel100_summary_list.rds") 



# new setup for mxTryHard function, 48k samples, lb = -.05, free a, rg, constraint on j, VF lb=.-1
summary_list <- readRDS("Analysis/Full_Model/m2_-.05lb_freeArg_VF-.1_forceJ_48000_nModelAll_summary_list.rds") 

# 25 model demonstration, new setup for mxTryHard function, 64k samples, lb = -.05, free a, rg, constraint on j, VF lb=.-.05
summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_VF-1_64000_nModelAll_summary_list.rds")

summary_list <- readRDS("Data/Full_Model/Local_Analysis/m2-.05lb_freeArg_VF-1_a.2_64000_nModelAll_summary_list.rds") #  bad a 

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

#remove the outliers that are three sd away from the mean of VY11 VY12 and VY22
df <- df[abs(df$VY11 - mean(df$VY11)) < 3*sd(df$VY11),]

# remove lower than 0 a
#df <- df[df$VF11 > 0 & df$VF22 > 0,]
nrow(df)


library(psych)
describe(df)
library(ggplot2)
library(tidyr)
#save_path <- "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Analysis/Full_Model/DistributionFigures"
# a plot for three VY estimates
true_values <- c(VY11 = 1.7478366, VY12 = 0.3401594,  VY22 = 1.1359723)
df_long <- tidyr::pivot_longer(df, c("VY11", "VY12", "VY22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
# Calculate means for VY variables
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
g1 <- ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") 
  #theme_minimal()
g1
#ggsave(paste0(save_path,"/VY_fixarg.png"), g1, width = 8, height = 8, type = "cairo-png", dpi = 400)

# a plot for three VF estimates
true_values <- c(VF11 = 0.17274308, VF12 = 0.08830186, VF22 = 0.05280104)
df_long <- tidyr::pivot_longer(df, c("VF11", "VF12", "VF22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
g1 <- ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") 

g1

# a plot for three VE estimates
true_values <- c(VE11 = 0.51000000, VE12 = 0.06545227, VE22 = 0.84000000)
df_long <- tidyr::pivot_longer(df, c("VE11", "VE12", "VE22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
g1 <- ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free")
g1

# a plot for four f estimates
describe(df[,c("f11","f12","f21","f22")], trim = 0)
true_values <- c(f11 = 0.15, f12 = 0.1, f21 = 0.05, f22 = 0.1)
df_long <- tidyr::pivot_longer(df, c("f11", "f12",  "f21", "f22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
g1 = ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") 
g1
#ggsave(paste0(save_path,"/f_fixarg.png"), g1, width = 8, height = 8, type = "cairo-png", dpi = 400)



# a plot for two delta estimates
true_values <- c(delta11 = sqrt(.49*.5), delta22 = sqrt(.16*.3))
df_long <- tidyr::pivot_longer(df, c("delta11", "delta22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for two a estimates
true_values <- c(a11 = sqrt(.49*.5), a22 = sqrt(.16*.7))
df_long <- tidyr::pivot_longer(df, c("a11", "a22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four v estimates
true_values <- c(v11 = 0.20765721, v12 = 0.10233074, v21 = 0.08560434, v22 = 0.07192265)
df_long <- tidyr::pivot_longer(df, c("v11", "v12",  "v21", "v22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
g1 = ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") 
g1
#ggsave(paste0(save_path,"/v_fixarg.png"), g1, width = 8, height = 8, type = "cairo-png", dpi = 400)


# a plot for four w estimates
true_values <- c(w11 = 0.20374099, w12 = 0.07627324, w21 = 0.08258827, w22 = 0.05076394)
df_long <- tidyr::pivot_longer(df, c("w11", "w12",  "w21", "w22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four Omega estimates
true_values <- c(Omega11 = 0.4337504, Omega12 = 0.08469132, Omega21 = 0.0634763, Omega22 = 0.1440672)
df_long <- tidyr::pivot_longer(df, c("Omega11", "Omega12",  "Omega21", "Omega22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four Gamma estimates
true_values <- c(Gamma11 = 0.4364954, Gamma12 = 0.09868633, Gamma21 = 0.0733564, Gamma22 = 0.2164179)
df_long <- tidyr::pivot_longer(df, c("Gamma11", "Gamma12",  "Gamma21", "Gamma22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for three gc estimates
true_values <- c(gc11 = 0.04243337, gc12 = 0.009898546, gc22 = 0.006565526)
df_long <- tidyr::pivot_longer(df, c("gc11", "gc12",  "gc22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  ylim(-.05,.05)+
  theme_minimal()

# a plot for three hc estimates
true_values <- c(hc11 = 0.04322826, hc12 = 0.01270896, hc22 = 0.01348196)
df_long <- tidyr::pivot_longer(df, c("hc11", "hc12",  "hc22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
g1 <- ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  ylim(-.05,.25)
g1
#ggsave(paste0(save_path,"/hc_fixarg.png"), g1, width = 8, height = 8, type = "cairo-png", dpi = 400)


# a plot for four mu estimates
true_values <- c(mu11 = 0.2226081, mu12 =0.02948167, mu21 = -0.04587792, mu22 = 0.2490391)
df_long <- tidyr::pivot_longer(df, c("mu11", "mu12",  "mu21", "mu22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# a plot for four ht estimates
true_values <- c(ht11 = 0.04322826, ht12 =0.01599563, ht21 = 0.009422301, ht22 = 0.01348196)
df_long <- tidyr::pivot_longer(df, c("ht11", "ht12",  "ht21", "ht22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()


# a plot for four gt estimates
true_values <- c(gt11 = 0.04243337, gt12 =0.01205056, gt21 = 0.007746527, gt22 = 0.006565526)
df_long <- tidyr::pivot_longer(df, c("gt11", "gt12",  "gt21", "gt22"), names_to = "Variable", values_to = "Value")
df_long$Index <- 1:nrow(df_long)
means <- aggregate(Value ~ Variable, data = df_long, FUN = median)
ggplot(df_long, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values[Variable]), color = "red") +
  geom_hline(data = means, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# plot(df$VY11, ylim = c(0,5))
# abline(h = 1.7292875, col = "red", lwd = 2)
# summary(df$VY11)


# plot(df$f22, ylim = c(0,1))
# abline(h = 0.10, col = "red", lwd = 2)
# # Now df is a data frame where each column is the estimates from each element in the summary_list










