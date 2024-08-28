# this script is designed for making plots using the .5 latent simulated data
library(psych)
library(ggplot2)
# read a summary list
# no constraints on gc and hc, 16000 samples, lb = -.05, smaller tolerance, new setup, no wvconstraint
summary_list1 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_16000_summary_list.rds")
summary_list2 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_32000_summary_list.rds")
summary_list3 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_48000_summary_list.rds")
summary_list4 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_64000_summary_list.rds")

#summary_list1 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_adelta_16000_summary_list.rds")
#summary_list2 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_adelta_32000_summary_list.rds")
#summary_list3 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_adelta_48000_summary_list.rds")
#summary_list4 <- readRDS("Analysis/Full_Model_.5latent/m2_.05lb_smallerTol_newSetup_adelta_64000_summary_list.rds")

save_path <- "Analysis/Full_Model_.5latent"
# a function to extract the estimates from the summary list
list2df <- function(summary_list){
    # extract all the status code of openmx and put them into a vector
    status_codes <- sapply(summary_list, function(x) x$statusCode)
    #summary(status_codes)
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
    return(df)
}

# extract the estimates from the summary list
df1 <- list2df(summary_list1)
df2 <- list2df(summary_list2)
df3 <- list2df(summary_list3)
df4 <- list2df(summary_list4)
df1$SampleSize <- "16k"
df2$SampleSize <- "32k"
df3$SampleSize <- "48k"
df4$SampleSize <- "64k"

# exclude the rows with a status code is not "OK" or "OK/green"
df1 <- df1[df1$status_codes %in% c("OK", "OK/green"),]
df2 <- df2[df2$status_codes %in% c("OK", "OK/green"),]
df3 <- df3[df3$status_codes %in% c("OK", "OK/green"),]
df4 <- df4[df4$status_codes %in% c("OK", "OK/green"),]

# get the SE of the estimates with four digits in four sample sizes
as.data.frame(describe(df1, trim = 0 ))[,c(1,4,13)]
as.data.frame(describe(df2, trim = 0 ))[,c(1,4,13)]
as.data.frame(describe(df3, trim = 0 ))[,c(1,4,13)]
as.data.frame(describe(df4, trim = 0 ))[,c(1,4,13)]
# pick the estimates of interest and paste them into a new dataframe 
df_plot <- rbind(df1, df2, df3,df4)
# exclude the rows with a status code is not "OK" or "OK/green"
df_plot <- df_plot[df_plot$status_codes %in% c("OK", "OK/green"),]

# only include rows that all estimates are positive
df_plot <- df_plot[rowSums(df_plot[,c(1:23,28:64)] < 0) == 0,]


# test the accuracy of mu
Vy_matrix <- matrix(c(mean(df_plot$VY11), mean(df_plot$VY12), mean(df_plot$VY12), mean(df_plot$VY22)), nrow = 2, byrow = T)
mu_matrix <- matrix(c(mean(df_plot$mu11), mean(df_plot$mu12), mean(df_plot$mu21), mean(df_plot$mu22)), nrow = 2, byrow = T)

Vy_matrix <- matrix(c(median(df_plot$VY11), median(df_plot$VY12), median(df_plot$VY12), median(df_plot$VY22)), nrow = 2, byrow = T)
mu_matrix <- matrix(c(median(df_plot$mu11), median(df_plot$mu12), median(df_plot$mu21), median(df_plot$mu22)), nrow = 2, byrow = T)

mu_matrix <- matrix(c(0.2207511, 0.01191484, -0.06126167, 0.2467678), nrow = 2, byrow = T)
matecov <- Vy_matrix%*%mu_matrix%*%t(Vy_matrix)
# put the Vy and Matecov into a 4 by 4 matrix
Vy_4by4 <- matrix(rbind(cbind(Vy_matrix, t(matecov)),
                        cbind(matecov, Vy_matrix)), nrow = 4, byrow = T)
Vy_4by4 |> cov2cor()
head(df_plot)

# make the plot for f11

library(ggplot2)

true_value <- 0.15
g1<-ggplot(df_plot, aes(x = SampleSize, y = f11)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.125, .175)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/f11.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# make the plot for hc11
true_value <- 0.01
g1<-ggplot(df_plot, aes(x = SampleSize, y = hc11)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.02, .06)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
# plot for f12

true_value <- 0.10
g1<-ggplot(df_plot, aes(x = SampleSize, y = f12)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.075, .125)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/f12.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for f22

true_value <- 0.10
g1<-ggplot(df_plot, aes(x = SampleSize, y = f22)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.075, .125)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/f22.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)


# plot for delta11
true_value <- sqrt(.486*.5)
g1<-ggplot(df_plot, aes(x = SampleSize, y = delta11)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.44, .54)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/delta11.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for delta22
true_value <- sqrt(.159*.5)
g1<-ggplot(df_plot, aes(x = SampleSize, y = delta22)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.24, .34)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/delta22.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for a11
true_value <- sqrt(.486*.5)
g1<-ggplot(df_plot, aes(x = SampleSize, y = a11)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.44, .54)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/a11.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for a22
true_value <- sqrt(.159*.5)
g1<-ggplot(df_plot, aes(x = SampleSize, y = a22)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.24, .34)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/a22.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for mu11
true_value <- 0.2207511
g1<-ggplot(df_plot, aes(x = SampleSize, y = mu11)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.2, .25)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/mu11.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for mu12
true_value <- 0.02661484
g1<-ggplot(df_plot, aes(x = SampleSize, y = mu12)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(0.01, .04)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/mu12.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)

# plot for k12
true_value <- .05
g1<-ggplot(df_plot, aes(x = SampleSize, y = k12)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  coord_cartesian(ylim = c(.035, .065)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
g1
ggsave(paste0(save_path,"/k12.png"), g1, width = 4, height = 6, type = "cairo-png", dpi = 600)




# Comparison with regression estimation of f11 and f12

# read the data
df_reg <- read.table("Data/Full_Model_.5latent/loop1.rds_64000.txt", sep = "\t", header = T)
colnames(df_reg)
df_reg$PGSp1 <- df_reg$NTp1 + df_reg$Tp1
df_reg$PGSo1 <- df_reg$Tp1 + df_reg$Tm1
df_reg$PGSp2 <- df_reg$NTp2 + df_reg$Tp2
df_reg$PGSo2 <- df_reg$Tp2 + df_reg$Tm2
# model for same-trait vertical transmission
m1 <- lm(Yo1 ~ Yp1 + Yo2 + PGSp1 + PGSo1, data = df_reg)
summary(m1)
m1$coefficients[2]

# model for cross-trait vertical transmission
m2 <- lm(Yo1 ~ Yp2 + Yp1+ PGSp1 + PGSo1 + PGSp2 + PGSo2, data = df_reg)
summary(m2)

# a loop to read all the data files, get the regression estimates and put them into a data frame
df_reg <- data.frame("f11" = as.numeric(rep(NA,125)), "f12" = as.numeric(rep(NA,125)))
# get all the file names in the datafolder
file_names <- list.files("Data/Full_Model_.5latent", pattern = "_64000", full.names = T)
for (loop in 1:125){
    df_reg_loop <- read.table(file_names[loop], sep = "\t", header = T)
    df_reg_loop$PGSp1 <- df_reg_loop$NTp1 + df_reg_loop$Tp1
    df_reg_loop$PGSo1 <- df_reg_loop$Tp1 + df_reg_loop$Tm1
    df_reg_loop$PGSp2 <- df_reg_loop$NTp2 + df_reg_loop$Tp2
    df_reg_loop$PGSo2 <- df_reg_loop$Tp2 + df_reg_loop$Tm2
    m1 <- lm(Yo1 ~ Yp1 +  PGSp1 + PGSo1, data = df_reg_loop)
    m2 <- lm(Yo1 ~ Yp2 + Yp1+ PGSp1 + PGSo1 + PGSp2 + PGSo2, data = df_reg_loop)
    df_reg[loop,1] <- m1$coefficients[2]
    df_reg[loop,2] <- m2$coefficients[2]
}

df_reg$type <- "Regression"
df_sem <- df_plot[,c("f11", "f12")]
df_sem$type <- "SEM"
df_compare <- rbind(df_reg, df_sem)
head(df_compare)

# make a plot for f11 comparison
true_value = 0.15
g1<-ggplot(df_compare, aes(x = type, y = f11)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  coord_cartesian(ylim = c(0.1, .5)) +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1))
g1
ggsave(paste0(save_path,"/f11_compare.png"), g1, width = 6, height = 6, type = "cairo-png", dpi = 600)

# make a plot for f12 comparison
true_value = 0.1
g1<-ggplot(df_compare, aes(x = type, y = f12)) +
  geom_boxplot(width = .4, fill = "#00a6ffd1", color = "#1313a3", size = 1, staplewidth = .2 ,outlier.shape = 5) +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "#2a6fef") +
  coord_cartesian(ylim = c(0.05, .15)) +
  geom_hline(aes(yintercept = true_value), color = "gold", size = 1.25, linetype = "24") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1))
g1
ggsave(paste0(save_path,"/f12_compare.png"), g1, width = 6, height = 6, type = "cairo-png", dpi = 400)

mean(df_plot$VY11)
