# Sensitivity Analysis for Fixed A Parameter - Analysis and Plotting
# Author: Generated for sensitivity analysis
# Date: August 5, 2025
# This script analyzes the results from the sensitivity analysis with fixed A parameter

# Load necessary libraries
library(ggplot2)
library(patchwork)
library(dplyr)

# Define the condition and sample size
condition_name <- "r2pgs04"
sample_size_name <- "32k"

# Define the A levels that were tested
a_level_names <- c("level1", "level2", "level3", "level4", "level5", 
                   "level6", "level7", "level8", "level9", "level10")

# Calculate the true A values for reference
true_a11 <- sqrt(0.64 * (0.60/0.64))  # sqrt(0.6) ≈ 0.7746
true_a22 <- sqrt(0.36 * 0.8)          # sqrt(0.288) ≈ 0.5367

# Calculate the actual A values for each level
a_values_list <- list(
    #level0 = c(true_a11 - 0.150, true_a22),   # -0.150 from true a11
    level1 = c(true_a11 - 0.125, true_a22),   # -0.125 from true a11
    level2 = c(true_a11 - 0.100, true_a22),   # -0.100 from true a11
    level3 = c(true_a11 - 0.075, true_a22),   # -0.075 from true a11
    level4 = c(true_a11 - 0.050, true_a22),   # -0.050 from true a11
    level5 = c(true_a11 - 0.025, true_a22),   # -0.025 from true a11
    level6 = c(true_a11, true_a22),           # True value for a11
    level7 = c(true_a11 + 0.025, true_a22),   # +0.025 from true a11
    level8 = c(true_a11 + 0.050, true_a22),   # +0.050 from true a11
    level9 = c(true_a11 + 0.075, true_a22),   # +0.075 from true a11
    level10 = c(true_a11 + 0.100, true_a22),   # +0.100 from true a11
    level11 = c(true_a11 + 0.125, true_a22)   # +0.125 from true a11
)

cat("Loading sensitivity analysis results...\n")

# Load the results for each A level
summary_lists <- list()
for (level_name in a_level_names) {
    file_path <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name, 
                       "/", level_name, "_fixedA_summary_list_100.rds")
    if (file.exists(file_path)) {
        summary_lists[[level_name]] <- readRDS(file_path)
        cat("Loaded", level_name, "\n")
    } else {
        cat("Warning: File not found:", file_path, "\n")
    }
}

# Function to extract parameter estimates from summary list
getDf <- function(summary_list) {
    status_codes <- sapply(summary_list, function(x) x$statusCode)
    df <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
    colnames(df) <- summary_list[[2]]$parameters$name
    
    # Loop over the elements in the summary_list
    for(i in 1:length(summary_list)) {
        for(j in 1:nrow(summary_list[[i]]$parameters)){
            df[i,j] <- summary_list[[i]]$parameters$Estimate[j]
        }
    }
    df$status_codes <- status_codes
    df <- df[df$status_codes %in% c("OK", "OK/green"),]
    return(df)
}

# Extract data frames for each A level
df_list <- list()
for (level_name in a_level_names) {
    if (level_name %in% names(summary_lists)) {
        df_list[[level_name]] <- getDf(summary_lists[[level_name]])
        cat("Extracted data for", level_name, "- rows:", nrow(df_list[[level_name]]), "\n")
    }
}

# Function to create plot data for a parameter across A levels
getDfPlot_FixedA <- function(param) {
    df_plot <- data.frame()
    
    for (level_name in a_level_names) {
        if (level_name %in% names(df_list)) {
            df_temp <- df_list[[level_name]]
            if (param %in% colnames(df_temp) && nrow(df_temp) > 0) {
                a11_value <- a_values_list[[level_name]][1]  # Get the a11 value for this level
                temp_data <- data.frame(
                    param_value = df_temp[[param]],
                    a11_fixed = rep(a11_value, nrow(df_temp)),
                    level = rep(level_name, nrow(df_temp))
                )
                df_plot <- rbind(df_plot, temp_data)
            }
        }
    }
    
    colnames(df_plot)[1] <- param
    return(df_plot)
}

# Function to get summary statistics
getDfSumm_FixedA <- function(df_plot, param) {
    df_summ <- aggregate(df_plot[[param]], by = list(df_plot$a11_fixed), FUN = median, na.rm = TRUE)
    colnames(df_summ) <- c("a11_fixed", "median")
    df_summ$MAD <- aggregate(df_plot[[param]], by = list(df_plot$a11_fixed), FUN = function(x) mad(x, na.rm = TRUE))[,2]
    return(df_summ)
}

# Function to extract true values from the level where A is at its true value (level6)
extractTrueValues <- function() {
    # level6 corresponds to fixing A at its true value
    true_level <- "level6"
    
    #if (true_level %in% names(df_list)) {
    if (FALSE) {
        df_true <- df_list[[true_level]]
        
        # Calculate median estimates for all parameters when A is at true value
        true_values <- list()
        param_names <- c("f11", "mu11", "w11", "v11", "delta11", "gc11",
                        "f12", "mu12", "w12", "v12", "VY12", "gc12")
        
        for (param in param_names) {
            if (param %in% colnames(df_true)) {
                true_values[[param]] <- median(df_true[[param]], na.rm = TRUE)
                cat("True value for", param, ":", round(true_values[[param]], 4), "\n")
            } else {
                cat("Warning: Parameter", param, "not found in data\n")
            }
        }
        
        return(true_values)
    } else {
        cat("Warning: level6 (true A value) not found in data. Using fallback values.\n")
        # Fallback to calculated values
        vg1 <- 0.64
        vg2 <- 0.36
        prop.h2.latent1 <- 0.60/0.64
        prop.h2.latent2 <- 0.8
        
        true_values <- list(
            f11 = 0.15,
            mu11 = 0.19695153,
            w11 = 0.09375280,
            v11 = 0.3473978,
            delta11 = sqrt(vg1*(1-prop.h2.latent1)),
            gc11 = 0.007098673,
            f21 = 0.05,
            mu21 = -0.06759629,
            w21 = 0.04014979,
            v21 = 0.1431288,
            VY12 = 0.4864288,
            gc12 = 0.002961545
        )
        return(true_values)
    }
}

# Extract true values from the level where A is fixed at its true value
true_values <- extractTrueValues()

# Function to create a violin plot for a parameter showing distribution at each A level
create_sensitivity_violin_plot <- function(param_name, color1 = "#1f77b4") {
    df_plot <- getDfPlot_FixedA(param_name)
    
    if (nrow(df_plot) == 0) {
        cat("Warning: No data for parameter", param_name, "\n")
        return(NULL)
    }
    
    # Get true value
    true_value <- true_values[[param_name]]
    if (is.null(true_value)) {
        # Fallback: use median when A is at true value
        df_true_level <- df_plot[abs(df_plot$a11_fixed - true_a11) < 0.001, ]
        if (nrow(df_true_level) > 0) {
            true_value <- median(df_true_level[[param_name]], na.rm = TRUE)
        } else {
            true_value <- median(df_plot[[param_name]], na.rm = TRUE)
        }
    }
    
    # Create the violin plot
    p <- ggplot(df_plot, aes(x = factor(round(a11_fixed, 3)), y = .data[[param_name]])) +
        geom_violin(fill = color1, alpha = 0.7, color = color1) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
        #geom_hline(yintercept = true_value, color = "#d62728", size = 1.25, linetype = "dashed") +
        geom_vline(xintercept = which(abs(sort(unique(df_plot$a11_fixed)) - true_a11) == min(abs(sort(unique(df_plot$a11_fixed)) - true_a11))), 
                   color = "#de1d1df0", size = 1.3, linetype = "dotted", alpha = 0.8) +
        labs(title = paste(param_name),
             x = expression(paste("Fixed ", a[11], " value")),
             y = paste(param_name)) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1)
        )
    
    return(p)
}

# Color palette
my_palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)

# Function to create combined violin plots
create_combined_sensitivity_violin_plot <- function(params, ncol = 3) {
    plots <- list()
    
    for (i in seq_along(params)) {
        param <- params[i]
        color <- my_palette[((i-1) %% length(my_palette)) + 1]
        plot <- create_sensitivity_violin_plot(param, color1 = color)
        if (!is.null(plot)) {
            plots[[param]] <- plot
        }
    }
    
    if (length(plots) > 0) {
        combined_plot <- wrap_plots(plots, ncol = ncol)
        return(combined_plot)
    } else {
        return(NULL)
    }
}

# Define the parameters for the two figures
params1 <- c("f11", "mu11", "w11", "v11", "delta11", "gc11")
params2 <- c("f21", "mu21",  "w21", "v21", "VY12","gc12")

cat("Creating first combined violin plot (11 parameters)...\n")
combined_plot1 <- create_combined_sensitivity_violin_plot(params1, ncol = 3)

if (!is.null(combined_plot1)) {
    print(combined_plot1)
    output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
    ggsave(paste0(output_dir, "/sensitivity_violin_plot_11_params.png"), 
           combined_plot1, width = 16, height = 8, type = "cairo-png", dpi = 600)
    cat("Saved first combined violin plot\n")
}

cat("Creating second combined violin plot (12 parameters)...\n")
combined_plot2 <- create_combined_sensitivity_violin_plot(params2, ncol = 3)

if (!is.null(combined_plot2)) {
    print(combined_plot2)
    output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
    ggsave(paste0(output_dir, "/sensitivity_violin_plot_12_params.png"), 
           combined_plot2, width = 16, height = 8, type = "cairo-png", dpi = 600)
    cat("Saved second combined violin plot\n")
}

cat("\n=== Sensitivity analysis plotting completed ===\n")
cat("Plots saved in:", paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name, "/"), "\n")
cat("Red dashed line: True parameter value\n")
cat("Green dotted line: True a11 value\n")