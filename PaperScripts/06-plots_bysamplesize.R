# the R script to create plots that estimates by sample size

# load the necessary libraries
library(ggplot2)
library(patchwork)
summary_list1 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_16000_summary_list.rds")
summary_list2 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_32000_summary_list.rds")
summary_list3 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_48000_summary_list.rds")
summary_list4 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_64000_summary_list.rds")
summary_list5 <- readRDS("Analysis/Paper/Model_latent30/m2_paper_80000_summary_list.rds")

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
    #df <- df[df$a11!=0.3 & df$a22!=0.3,]
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

# plot the estimates of f11 as a function of the sample size
getDfPlot <- function(param){
    df_plot <- rbind(cbind(df1[[param]], df_se1[[param]], rep(16000, nrow(df1))),
                     cbind(df2[[param]], df_se2[[param]], rep(32000, nrow(df2))),
                     cbind(df3[[param]], df_se3[[param]], rep(48000, nrow(df3))),
                     cbind(df4[[param]], df_se4[[param]], rep(64000, nrow(df4))),
                     cbind(df5[[param]], df_se5[[param]], rep(80000, nrow(df5))))
    df_plot <- as.data.frame(df_plot)
    colnames(df_plot) <- c(param,paste0("se_", param), "sample_size")
    df_plot$sample_size <- as.factor(df_plot$sample_size)
    return(df_plot)
}

getDfSumm <- function(df_plot, func = "mean"){
    if (func == "median") {
        df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = median)
        colnames(df_summ) <- c("sample_size", "median")

    } else if (func == "mean") {
        df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = mean)
        colnames(df_summ) <- c("sample_size", "mean")

    }
    #df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = mean)
    df_summ$se <- aggregate(df_plot[,2], by = list(df_plot$sample_size), FUN = mean)[,2]
    return(df_summ)
}

# Define a function to create a prettier plot
create_pretty_plot <- function(df_summ, param_name, color1 = "blue") {
  ggplot(df_summ, aes(x = sample_size, y = mean)) +
    geom_point(size = 3, color = color1) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = color1) +
    geom_line(aes(group = 1), color = color1, size = 1) +
    labs(title = paste("Estimates of", param_name),
         x = "Sample size",
         y = paste(param_name)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    scale_color_manual(values = c(color1))
}

# Get the data for the plot
df_plot <- getDfPlot("delta11")
df_summ <- getDfSumm(df_plot)
# Create the plot
pretty_plot <- create_pretty_plot(df_summ, "delta11")

# Display the plot
print(pretty_plot)

# a color palette for the plots
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")
# Define a function to create plots for a selection of parameters and combine them
create_combined_plot <- function(params, ncol = 2) {
  plots <- list()
  
  for (param in params) {
    df_plot <- getDfPlot(param)
    df_summ <- getDfSumm(df_plot)
    plot <- create_pretty_plot(df_summ, param, color1 = sample(my_palette, 1))
    plots[[param]] <- plot
  }
  
  combined_plot <- wrap_plots(plots, ncol = ncol)
  return(combined_plot)
}

# Define the parameters you want to plot
params <- c("delta11", "f11", "a11", "mu11")

# Create the combined plot
combined_plot <- create_combined_plot(params)

# Display the combined plot
print(combined_plot)

params2 <- c("w12","v12","VY12","gc12","hc12","ic12")
combined_plot2 <- create_combined_plot(params2)
print(combined_plot2)
