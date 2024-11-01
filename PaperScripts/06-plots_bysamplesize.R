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
    se_y <- max(df_summ$se, na.rm=TRUE)*4
    lim_y <- c(mean(df_summ$mean - se_y), mean(df_summ$mean + se_y))
  ggplot(df_summ, aes(x = sample_size, y = mean)) +
    geom_point(size = 3, color = color1) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = color1) +
    geom_line(aes(group = 1), color = color1, size = 1) +
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

# Get the data for the plot
df_plot <- getDfPlot("a11")
#df_summ <- aggregate(df_plot[,1], by = list(df_plot$sample_size), FUN = mean)
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
  
  combined_plot <- wrap_plots(plots, ncol = ncol, widths = rep(1, length(params)))
  return(combined_plot)
}

# Define the parameters you want to plot
params <- c("f11", "mu11", "a11")
params <- c("f11", "mu11", "a11", "w11","VY11","gc11")


# Create the combined plot
combined_plot <- create_combined_plot(params, ncol = 3)

# Display the combined plot
print(combined_plot)
ggsave("Analysis/Paper/p1.png", combined_plot, width = 10, height = 6, type = "cairo-png", dpi = 600)

params2 <- c("f12", "mu12", "k12", "w12","VY12","gc12")

combined_plot2 <- create_combined_plot(params2, ncol = 3)
print(combined_plot2)
ggsave("Analysis/Paper/p2.png", combined_plot2, width = 10, height = 6, type = "cairo-png", dpi = 600)

