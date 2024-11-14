# the R script to create plots that estimates by sample size

# load the necessary libraries
library(ggplot2)
library(patchwork)
summary_list1 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_16000_summary_list.rds")
summary_list2 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_32000_summary_list.rds")
summary_list3 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_48000_summary_list.rds")
summary_list4 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_64000_summary_list.rds")
summary_list5 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_80000_summary_list.rds")

summary_list1_fixedA <- readRDS("Analysis/Paper/Model_latent30/m2_paper_16000_summary_list_fixedA.rds")
summary_list2_fixedA <- readRDS("Analysis/Paper/Model_latent30/m2_paper_32000_summary_list_fixedA.rds")
summary_list3_fixedA <- readRDS("Analysis/Paper/Model_latent30/m2_paper_48000_summary_list_fixedA.rds")
summary_list4_fixedA <- readRDS("Analysis/Paper/Model_latent30/m2_paper_64000_summary_list_fixedA.rds")
summary_list5_fixedA <- readRDS("Analysis/Paper/Model_latent30/m2_paper_80000_summary_list_fixedA.rds")

getDf <- function(summary_list) {
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
    df <- df[df$a11!=0.2 & df$a22!=0.2,]

    # general lower bound variables
    varname_vector <- c()
    return(df)
}

getSe <- function(summary_list){
    status_codes <- sapply(summary_list, function(x) x$statusCode)
    df_se <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
    colnames(df_se) <- summary_list[[2]]$parameters$name
    colnames(df_se)
    # Loop over the elements in the summary_list
    for(i in 1:length(summary_list)) {
        for(j in 1:nrow(summary_list[[i]]$parameters)){
            if (!is.null(summary_list[[i]]$parameters$Std.Error[j])) {
                df_se[i,j] <- summary_list[[i]]$parameters$Std.Error[j]
            } else {
                print(paste("NULL value at i =", i, "and j =", j))
            }
        }
    }
    df_se$status_codes <- status_codes
    df_se <- df_se[df_se$status_codes %in% c("OK", "OK/green"),]
    #df_se <- df_se[df_se$a11!=0.3 & df_se$a22!=0.3,]
    return(df_se)
}


df1 <- getDf(summary_list1)
df2 <- getDf(summary_list2)
df3 <- getDf(summary_list3)
df4 <- getDf(summary_list4)
df5 <- getDf(summary_list5)

df_se1 <- getSe(summary_list1)
df_se2 <- getSe(summary_list2)
df_se3 <- getSe(summary_list3)
df_se4 <- getSe(summary_list4)
df_se5 <- getSe(summary_list5)

df1_fixedA <- getDf(summary_list1_fixedA)
df2_fixedA <- getDf(summary_list2_fixedA)
df3_fixedA <- getDf(summary_list3_fixedA)
df4_fixedA <- getDf(summary_list4_fixedA)
df5_fixedA <- getDf(summary_list5_fixedA)

df_se1_fixedA <- getSe(summary_list1_fixedA)
df_se2_fixedA <- getSe(summary_list2_fixedA)
df_se3_fixedA <- getSe(summary_list3_fixedA)
df_se4_fixedA <- getSe(summary_list4_fixedA)
df_se5_fixedA <- getSe(summary_list5_fixedA)

# plot the estimates of f11 as a function of the sample size
getDfPlot <- function(param){
    df_plot <- rbind(cbind(df1[[param]], df_se1[[param]], rep("16k", nrow(df1))),
                     cbind(df2[[param]], df_se2[[param]], rep("32k", nrow(df2))),
                     cbind(df3[[param]], df_se3[[param]], rep("48k", nrow(df3))),
                     cbind(df4[[param]], df_se4[[param]], rep("64k", nrow(df4))),
                     cbind(df5[[param]], df_se5[[param]], rep("80k", nrow(df5))))
    df_plot <- as.data.frame(df_plot)
    df_plot[,3] <- as.factor(df_plot[,3])
    df_plot[,1:2] <- apply(df_plot[,1:2], 2, as.numeric)
    colnames(df_plot) <- c(param,paste0("se_", param), "sample_size")
    df_plot$sample_size <- as.factor(df_plot$sample_size)
    return(df_plot)
}

getDfPlot_fixedA <- function(param){
    df_plot <- rbind(cbind(df1_fixedA[[param]], df_se1_fixedA[[param]], rep("16k", nrow(df1_fixedA))),
                     cbind(df2_fixedA[[param]], df_se2_fixedA[[param]], rep("32k", nrow(df2_fixedA))),
                     cbind(df3_fixedA[[param]], df_se3_fixedA[[param]], rep("48k", nrow(df3_fixedA))),
                     cbind(df4_fixedA[[param]], df_se4_fixedA[[param]], rep("64k", nrow(df4_fixedA))),
                     cbind(df5_fixedA[[param]], df_se5_fixedA[[param]], rep("80k", nrow(df5_fixedA))))
    df_plot <- as.data.frame(df_plot)
    df_plot[,3] <- as.factor(df_plot[,3])
    df_plot[,1:2] <- apply(df_plot[,1:2], 2, as.numeric)
    colnames(df_plot) <- c(param,paste0("se_", param), "sample_size")
    df_plot$sample_size <- as.factor(df_plot$sample_size)
    return(df_plot)
}

getDfSumm <- function(df_plot, func = "median"){
    if (func == "median") {
        df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = function(x) median(x, na.rm = TRUE))
        colnames(df_summ) <- c("sample_size", "center")

    } else if (func == "mean") {
        df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = function(x) mean(x, na.rm = TRUE))
        colnames(df_summ) <- c("sample_size", "center")

    }
    #df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = mean)
    df_summ$se <- aggregate(df_plot[,2], by = list(df_plot$sample_size), FUN = function(x) mean(x, na.rm = TRUE))[,2]
    return(df_summ)
}

# Define a function to create a prettier plot
create_pretty_plot <- function(df_summ, param_name, color1 = "blue", file_tv = "Data/Paper/Expected/Model_latent30_finalGen.txt") {
    true_value_df <- read.table(file_tv, header = FALSE)
    
    true_value <- true_value_df[true_value_df$V1 == param_name, 2]
    print(true_value)
    se_y <- max(df_summ$se, na.rm=TRUE)*4
    lim_y <- c(mean(df_summ$center - se_y), mean(df_summ$center + se_y))
    ggplot(df_summ, aes(x = sample_size, y = center)) +
    geom_point(size = 3, color = color1) +
    geom_errorbar(aes(ymin = center - se, ymax = center + se), width = 0.2, color = color1) +
    geom_line(aes(group = 1), color = color1, size = 1) +
    geom_hline(aes(yintercept = true_value), color = "#b2182b", size = 1.25, linetype = "24") +
    coord_cartesian(ylim = lim_y) +
    labs(title = paste("Estimates of", param_name),
         x = "Sample size",
         y = paste(param_name)) +
      theme_minimal() +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black", size = 1),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_text(size = 10),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 10))+
    scale_color_manual(values = c(color1))
}

#library(ggplot2)

# Modified create_pretty_plot function
create_pretty_plot_se <- function(df_summ_fixed, df_summ_nonfixed, param_name, 
                               color_fixed = "#1f77b4", color_nonfixed = "#ad494a", 
                               file_tv = "Data/Paper/Expected/Model_latent30_finalGen.txt") {
  
  # Add a group identifier to each dataframe
  df_summ_fixed$A_status <- "Fixed a"
  df_summ_nonfixed$A_status <- "Estimated a"
  
  # Combine the dataframes
  df_combined <- dplyr::bind_rows(df_summ_fixed, df_summ_nonfixed)
  
  # Read the true value from the file
  true_value_df <- read.table(file_tv, header = FALSE)
  true_value <- true_value_df[true_value_df$V1 == param_name, 2]
  
  # Calculate y limits based on SE
  se_y <- max(df_combined$se, na.rm = TRUE) * 4
  lim_y <- c(mean(df_combined$center - se_y), mean(df_combined$center + se_y))
  
  # Define color palette
  my_palette <- c("Fixed a" = color_fixed, "Estimated a" = color_nonfixed)
  
  # Create the plot
  p <- ggplot(df_combined, aes(x = sample_size, y = se, color = A_status)) +
    geom_point(size = 3) +
    #geom_errorbar(aes(ymin = center - se, ymax = center + se), width = 0.2) +
    geom_line(aes(group = A_status), size = 1) +
    #geom_hline(yintercept = true_value, linetype = "dashed", color = "#b2182b", size = 1.25) +
    #coord_cartesian(ylim = lim_y) +
    labs(title = paste("SE of", param_name),
         x = "Sample Size",
         y = paste(param_name, "SE"),
         color = "a Status") +
    scale_color_manual(values = my_palette) +
    guides(color = "none")+
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  return(p)
}

# df_plot <- getDfPlot("f11")
# df_summ <- getDfSumm(df_plot)
# df_plot_fixedA <- getDfPlot_fixedA("f11")
# df_summ_fixedA <- getDfSumm(df_plot_fixedA)
# create_pretty_plot_se(df_summ_fixedA, df_summ, "f11")

# Get the data for the plot
#df_plot <- getDfPlot("a11")
#df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = mean)
#df_summ <- getDfSumm(df_plot)
# Create the plot
#pretty_plot <- create_pretty_plot(df_summ, "delta11")

# Display the plot
print(pretty_plot)

# a color palette for the plots
#my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")
my_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#6b6ecf", "#b5cf6b", "#9c9ede", "#e7969c", "#cedb9c",
  "#e7ba52", "#9edae5", "#dbdb8d", "#ad494a", "#393b79"
)
# Define a function to create plots for a selection of parameters and combine them
create_combined_plot <- function(params, ncol = 2) {
  plots <- list()
  
  for (param in params) {
    df_plot <- getDfPlot(param)
    df_summ <- getDfSumm(df_plot)
    plot <- create_pretty_plot(df_summ, param, color1 = sample(my_palette, 1))
    plots[[param]] <- plot
  }
  
  combined_plot <- wrap_plots(plots, ncol = ncol, widths = rep(1, length(params)))
  return(combined_plot)
}

# Define the parameters you want to plot
params <- c("f11", "mu11")
params <- c("f11", "mu11", "delta11", "w11","VY11","gc11")


# Create the combined plot
combined_plot <- create_combined_plot(params, ncol = 3)

# Display the combined plot
print(combined_plot)
ggsave("Analysis/Paper/p1.png", combined_plot, width = 10, height = 6, type = "cairo-png", dpi = 600)

params2 <- c("f12", "mu12", "v12" ,"w12","VY12","gc12")

combined_plot2 <- create_combined_plot(params2, ncol = 3)
print(combined_plot2)
ggsave("Analysis/Paper/p2.png", combined_plot2, width = 10, height = 6, type = "cairo-png", dpi = 600)


create_combined_plot_se <- function(params, ncol = 2) {
  plots <- list()
  
  for (param in params) {
    df_plot <- getDfPlot(param)
    df_summ <- getDfSumm(df_plot)
    df_plot_fixedA <- getDfPlot_fixedA(param)
    df_summ_fixedA <- getDfSumm(df_plot_fixedA)
    plot <- create_pretty_plot_se(df_summ_fixedA, df_summ, param)
    plots[[param]] <- plot
  }
  
  combined_plot <- wrap_plots(plots, ncol = ncol, widths = rep(1, length(params)))
  return(combined_plot)
}

combined_plot_se <- create_combined_plot_se(params, ncol = 3)
print(combined_plot_se)
ggsave("Analysis/Paper/p5.png", combined_plot_se, width = 10, height = 6, type = "cairo-png", dpi = 600)

combined_plot_se2 <- create_combined_plot_se(params2, ncol = 3)
print(combined_plot_se2)
ggsave("Analysis/Paper/p6.png", combined_plot_se2, width = 10, height = 6, type = "cairo-png", dpi = 600)
