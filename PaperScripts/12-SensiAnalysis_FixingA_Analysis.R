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
        param_names <- c("f11", "f12", "f21", "f22", "mu11", "mu12", "mu21", "mu22", 
                        "w11", "w21", "v11", "v21", "delta11", "gc11", "gc12",
                        "VY11", "VY12", "VY22")
        
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
            f12 = 0.1,
            f21 = 0.05,
            f22 = 0.1,
            mu11 = 0.19695153,
            mu12 = 0.1,
            mu21 = -0.06759629,
            mu22 = 0.3,
            w11 = 0.09375280,
            w21 = 0.04014979,
            v11 = 0.3473978,
            v21 = 0.1431288,
            delta11 = sqrt(vg1*(1-prop.h2.latent1)),
            gc11 = 0.007098673,
            gc12 = 0.002961545,
            VY11 = 2.0,
            VY12 = 0.4864288,
            VY22 = 1.5
        )
        return(true_values)
    }
}

# Extract true values from the level where A is fixed at its true value
true_values <- extractTrueValues()

# Function to create a violin plot for a parameter showing distribution at each A level
create_sensitivity_violin_plot <- function(param_name, color1 = "#1f77b4") {
    df_plot <- getDfPlot_FixedA(param_name)
    
    # Add VF matrix elements as columns if we have the necessary parameters
    # Get all necessary parameters for VF calculation
    df_f11 <- getDfPlot_FixedA("f11")
    df_f12 <- getDfPlot_FixedA("f12")
    df_f21 <- getDfPlot_FixedA("f21") 
    df_f22 <- getDfPlot_FixedA("f22")
    df_mu11 <- getDfPlot_FixedA("mu11")
    df_mu12 <- getDfPlot_FixedA("mu12")
    df_mu21 <- getDfPlot_FixedA("mu21")
    df_mu22 <- getDfPlot_FixedA("mu22")
    df_VY11 <- getDfPlot_FixedA("VY11")
    df_VY12 <- getDfPlot_FixedA("VY12")
    df_VY22 <- getDfPlot_FixedA("VY22")
    
    # Check if all necessary parameters exist and have the same structure
    if (nrow(df_f11) > 0 && nrow(df_f12) > 0 && nrow(df_f21) > 0 && nrow(df_f22) > 0 &&
        nrow(df_mu11) > 0 && nrow(df_mu12) > 0 && nrow(df_mu21) > 0 && nrow(df_mu22) > 0 &&
        nrow(df_VY11) > 0 && nrow(df_VY12) > 0 && nrow(df_VY22) > 0) {
        
        # Merge all parameters by a11_fixed and level to ensure proper alignment
        all_params <- merge(df_f11, df_f12, by = c("a11_fixed", "level"), suffixes = c("", ".f12"))
        all_params <- merge(all_params, df_f21, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_f22, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_mu11, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_mu12, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_mu21, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_mu22, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_VY11, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_VY12, by = c("a11_fixed", "level"))
        all_params <- merge(all_params, df_VY22, by = c("a11_fixed", "level"))
        
        # Calculate VF matrix elements for each row
        # VF = 2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f)
        VF11_vals <- VF12_vals <- VF22_vals <- numeric(nrow(all_params))
        
        for (i in 1:nrow(all_params)) {
            # Construct matrices for this row
            f_mat <- matrix(c(all_params$f11[i], all_params$f12[i], 
                             all_params$f21[i], all_params$f22[i]), nrow = 2, byrow = TRUE)
            
            mu_mat <- matrix(c(all_params$mu11[i], all_params$mu12[i],
                              all_params$mu21[i], all_params$mu22[i]), nrow = 2, byrow = TRUE)
            
            VY_mat <- matrix(c(all_params$VY11[i], all_params$VY12[i],
                              all_params$VY12[i], all_params$VY22[i]), nrow = 2, byrow = TRUE)
            
            # Calculate VF matrix: VF = 2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f)
            term1 <- 2 * f_mat %*% VY_mat %*% t(f_mat)
            term2 <- f_mat %*% VY_mat %*% mu_mat %*% VY_mat %*% t(f_mat)
            term3 <- f_mat %*% VY_mat %*% t(mu_mat) %*% VY_mat %*% t(f_mat)
            
            VF_mat <- term1 + term2 + term3
            
            # Store the matrix elements
            VF11_vals[i] <- VF_mat[1, 1]
            VF12_vals[i] <- VF_mat[1, 2]
            VF22_vals[i] <- VF_mat[2, 2]
        }
        
        # Add VF elements to the original df_plot if it matches the structure
        if (nrow(df_plot) == nrow(all_params)) {
            df_plot$VF11 <- VF11_vals
            df_plot$VF12 <- VF12_vals  
            df_plot$VF22 <- VF22_vals
        }
    } 

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

# Function to create VF data for plotting (memory optimized)
create_VF_data <- function() {
    # Get all necessary parameters for VF calculation
    cat("Loading parameter data for VF calculation...\n")
    
    # Check if required parameters exist first
    required_params <- c("f11", "f12", "f21", "f22", "mu11", "mu12", "mu21", "mu22", "VY11", "VY12", "VY22")
    param_check <- sapply(required_params, function(p) {
        df_temp <- getDfPlot_FixedA(p)
        return(nrow(df_temp) > 0)
    })
    
    if (!all(param_check)) {
        missing_params <- required_params[!param_check]
        cat("Warning: Missing required parameters for VF analysis:", paste(missing_params, collapse = ", "), "\n")
        return(NULL)
    }
    
    # Get base structure from f11 (smallest dataset to work with)
    df_base <- getDfPlot_FixedA("f11")
    n_total <- nrow(df_base)
    cat("Processing", n_total, "observations for VF calculation\n")
    
    # Sample data if too large to prevent memory issues
    if (n_total > 5000) {
        cat("Large dataset detected. Sampling 5000 observations to prevent memory issues...\n")
        sample_idx <- sample(n_total, 5000)
        df_base <- df_base[sample_idx, ]
        n_total <- 5000
    }
    
    # Get all parameters but only for the sampled/filtered observations
    param_data <- list()
    for (param in required_params) {
        df_temp <- getDfPlot_FixedA(param)
        if (n_total < nrow(df_temp)) {
            # Match the sampling if we reduced the dataset
            df_temp <- df_temp[sample_idx, ]
        }
        param_data[[param]] <- df_temp
    }
    
    # Pre-allocate result vectors
    VF11_vals <- numeric(n_total)
    VF12_vals <- numeric(n_total)
    VF22_vals <- numeric(n_total)
    
    # Process in smaller chunks to save memory
    chunk_size <- 500
    n_chunks <- ceiling(n_total / chunk_size)
    
    cat("Processing VF calculations in", n_chunks, "chunks...\n")
    
    for (chunk in 1:n_chunks) {
        start_idx <- (chunk - 1) * chunk_size + 1
        end_idx <- min(chunk * chunk_size, n_total)
        chunk_indices <- start_idx:end_idx
        
        cat("Processing chunk", chunk, "of", n_chunks, "(rows", start_idx, "to", end_idx, ")\n")
        
        # Calculate VF for this chunk
        for (i in chunk_indices) {
            # Get parameter values for this observation
            f11 <- param_data$f11[[1]][i]
            f12 <- param_data$f12[[1]][i]  
            f21 <- param_data$f21[[1]][i]
            f22 <- param_data$f22[[1]][i]
            
            mu11 <- param_data$mu11[[1]][i]
            mu12 <- param_data$mu12[[1]][i]
            mu21 <- param_data$mu21[[1]][i]
            mu22 <- param_data$mu22[[1]][i]
            
            VY11 <- param_data$VY11[[1]][i]
            VY12 <- param_data$VY12[[1]][i]
            VY22 <- param_data$VY22[[1]][i]
            
            # Construct matrices (more memory efficient)
            # VF = 2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f)
            
            # Calculate VF elements directly without creating full matrices
            # For a 2x2 matrix, we can compute elements individually
            
            # Term 1: 2 * f %*% VY %*% t(f)
            # f %*% VY = [f11*VY11 + f12*VY12, f11*VY12 + f12*VY22; f21*VY11 + f22*VY12, f21*VY12 + f22*VY22]
            fVY11 <- f11*VY11 + f12*VY12
            fVY12 <- f11*VY12 + f12*VY22
            fVY21 <- f21*VY11 + f22*VY12  
            fVY22 <- f21*VY12 + f22*VY22
            
            # (f %*% VY) %*% t(f) 
            term1_11 <- 2 * (fVY11*f11 + fVY12*f21)
            term1_12 <- 2 * (fVY11*f12 + fVY12*f22)
            term1_22 <- 2 * (fVY21*f12 + fVY22*f22)
            
            # Term 2: f %*% VY %*% mu %*% VY %*% t(f)
            # VY %*% mu = [VY11*mu11 + VY12*mu21, VY11*mu12 + VY12*mu22; VY12*mu11 + VY22*mu21, VY12*mu12 + VY22*mu22]
            VYmu11 <- VY11*mu11 + VY12*mu21
            VYmu12 <- VY11*mu12 + VY12*mu22
            VYmu21 <- VY12*mu11 + VY22*mu21
            VYmu22 <- VY12*mu12 + VY22*mu22
            
            # (VY %*% mu) %*% VY = above result multiplied by VY
            VYmuVY11 <- VYmu11*VY11 + VYmu12*VY12
            VYmuVY12 <- VYmu11*VY12 + VYmu12*VY22
            VYmuVY21 <- VYmu21*VY11 + VYmu22*VY12
            VYmuVY22 <- VYmu21*VY12 + VYmu22*VY22
            
            # f %*% (VY %*% mu %*% VY)
            fVYmuVY11 <- f11*VYmuVY11 + f12*VYmuVY21
            fVYmuVY12 <- f11*VYmuVY12 + f12*VYmuVY22
            fVYmuVY21 <- f21*VYmuVY11 + f22*VYmuVY21
            fVYmuVY22 <- f21*VYmuVY12 + f22*VYmuVY22
            
            # (f %*% VY %*% mu %*% VY) %*% t(f)
            term2_11 <- fVYmuVY11*f11 + fVYmuVY12*f21
            term2_12 <- fVYmuVY11*f12 + fVYmuVY12*f22
            term2_22 <- fVYmuVY21*f12 + fVYmuVY22*f22
            
            # Term 3: f %*% VY %*% t(mu) %*% VY %*% t(f) (similar calculation with mu transposed)
            # VY %*% t(mu) where t(mu) has mu11, mu21 in first row and mu12, mu22 in second row
            VYmut11 <- VY11*mu11 + VY12*mu12
            VYmut12 <- VY11*mu21 + VY12*mu22
            VYmut21 <- VY12*mu11 + VY22*mu12
            VYmut22 <- VY12*mu21 + VY22*mu22
            
            # Continue similar calculations...
            VYmutVY11 <- VYmut11*VY11 + VYmut12*VY12
            VYmutVY12 <- VYmut11*VY12 + VYmut12*VY22
            VYmutVY21 <- VYmut21*VY11 + VYmut22*VY12
            VYmutVY22 <- VYmut21*VY12 + VYmut22*VY22
            
            fVYmutVY11 <- f11*VYmutVY11 + f12*VYmutVY21
            fVYmutVY12 <- f11*VYmutVY12 + f12*VYmutVY22
            fVYmutVY21 <- f21*VYmutVY11 + f22*VYmutVY21
            fVYmutVY22 <- f21*VYmutVY12 + f22*VYmutVY22
            
            term3_11 <- fVYmutVY11*f11 + fVYmutVY12*f21
            term3_12 <- fVYmutVY11*f12 + fVYmutVY12*f22
            term3_22 <- fVYmutVY21*f12 + fVYmutVY22*f22
            
            # Final VF matrix elements
            VF11_vals[i] <- term1_11 + term2_11 + term3_11
            VF12_vals[i] <- term1_12 + term2_12 + term3_12
            VF22_vals[i] <- term1_22 + term2_22 + term3_22
        }
        
        # Force garbage collection after each chunk
        gc()
    }
    
    # Create result data frames
    df_VF11 <- data.frame(
        VF11 = VF11_vals,
        a11_fixed = df_base$a11_fixed,
        level = df_base$level
    )
    
    df_VF12 <- data.frame(
        VF12 = VF12_vals,
        a11_fixed = df_base$a11_fixed,
        level = df_base$level
    )
    
    df_VF22 <- data.frame(
        VF22 = VF22_vals,
        a11_fixed = df_base$a11_fixed,
        level = df_base$level
    )
    
    cat("VF calculation completed successfully\n")
    return(list(VF11 = df_VF11, VF12 = df_VF12, VF22 = df_VF22))
}

# Function to create VF-specific violin plots
create_VF_violin_plots <- function() {
    VF_data <- create_VF_data()
    
    if (is.null(VF_data)) {
        cat("Cannot create VF plots: missing data\n")
        return(NULL)
    }
    
    plots <- list()
    colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")
    
    for (i in 1:2) {
        param_name <- names(VF_data)[i]
        df_plot <- VF_data[[param_name]]
        color <- colors[i]
        
        p <- ggplot(df_plot, aes(x = factor(round(a11_fixed, 3)), y = .data[[param_name]])) +
            geom_violin(fill = color, alpha = 0.7, color = color) +
            geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
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
        
        plots[[param_name]] <- p
    }
    
    combined_plot <- wrap_plots(plots, ncol = 2)
    return(combined_plot)
}

# Function to create combined violin plots for regular parameters
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

# Define the parameters for the figures
params1 <- c("f11", "f12", "f21", "f22", "mu11", "mu12")
params2 <- c("mu21", "mu22", "VY11", "VY12", "VY22", "w11")
params3 <- c("w21", "v11", "v21", "delta11", "gc11", "gc12")

# Define VF parameters as a separate group for analysis
params_VF <- c("VF11", "VF12")

cat("Creating first combined violin plot (f and mu parameters)...\n")
#combined_plot1 <- create_combined_sensitivity_violin_plot(params1, ncol = 3)

# if (!is.null(combined_plot1)) {
#     print(combined_plot1)
#     output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
#     ggsave(paste0(output_dir, "/sensitivity_violin_plot_f_mu_params.png"), 
#            combined_plot1, width = 16, height = 8, type = "cairo-png", dpi = 600)
#     cat("Saved first combined violin plot\n")
# }

# cat("Creating second combined violin plot (VY and other parameters)...\n")
# combined_plot2 <- create_combined_sensitivity_violin_plot(params2, ncol = 3)

# if (!is.null(combined_plot2)) {
#     print(combined_plot2)
#     output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
#     ggsave(paste0(output_dir, "/sensitivity_violin_plot_VY_other_params.png"), 
#            combined_plot2, width = 16, height = 8, type = "cairo-png", dpi = 600)
#     cat("Saved second combined violin plot\n")
# }

# cat("Creating third combined violin plot (remaining parameters)...\n")
# combined_plot3 <- create_combined_sensitivity_violin_plot(params3, ncol = 3)

# if (!is.null(combined_plot3)) {
#     print(combined_plot3)
#     output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
#     ggsave(paste0(output_dir, "/sensitivity_violin_plot_remaining_params.png"), 
#            combined_plot3, width = 16, height = 8, type = "cairo-png", dpi = 600)
#     cat("Saved third combined violin plot\n")
# }

cat("Creating VF matrix violin plot...\n")
combined_plot_VF <- create_VF_violin_plots()

if (!is.null(combined_plot_VF)) {
    print(combined_plot_VF)
    output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    ggsave(paste0(output_dir, "/sensitivity_violin_plot_VF_matrix.png"), 
           combined_plot_VF, width = 16, height = 6, type = "cairo-png", dpi = 600)
    cat("Saved VF matrix violin plot\n")
}

# Additional analysis: Create a summary of VF elements as a function of fixed a11
cat("Creating VF summary analysis...\n")

# Function to create VF summary data across all a11 levels using the computed VF data
create_VF_summary_from_data <- function() {
    VF_data <- create_VF_data()
    
    if (is.null(VF_data)) {
        cat("Warning: No VF data available for summary analysis\n")
        return(NULL)
    }
    
    # Get unique a11_fixed values from any of the VF data frames
    a11_levels <- sort(unique(VF_data$VF11$a11_fixed))
    
    # Calculate median and MAD VF values for each a11 level
    VF_summary <- data.frame(
        a11_fixed = a11_levels,
        VF11_median = numeric(length(a11_levels)),
        VF12_median = numeric(length(a11_levels)), 
        VF22_median = numeric(length(a11_levels)),
        VF11_mad = numeric(length(a11_levels)),
        VF12_mad = numeric(length(a11_levels)),
        VF22_mad = numeric(length(a11_levels))
    )
    
    for (i in seq_along(a11_levels)) {
        a11_val <- a11_levels[i]
        
        # Calculate summary statistics for each VF element
        VF_summary$VF11_median[i] <- median(VF_data$VF11$VF11[VF_data$VF11$a11_fixed == a11_val], na.rm = TRUE)
        VF_summary$VF12_median[i] <- median(VF_data$VF12$VF12[VF_data$VF12$a11_fixed == a11_val], na.rm = TRUE)
        VF_summary$VF22_median[i] <- median(VF_data$VF22$VF22[VF_data$VF22$a11_fixed == a11_val], na.rm = TRUE)
        
        VF_summary$VF11_mad[i] <- mad(VF_data$VF11$VF11[VF_data$VF11$a11_fixed == a11_val], na.rm = TRUE)
        VF_summary$VF12_mad[i] <- mad(VF_data$VF12$VF12[VF_data$VF12$a11_fixed == a11_val], na.rm = TRUE)
        VF_summary$VF22_mad[i] <- mad(VF_data$VF22$VF22[VF_data$VF22$a11_fixed == a11_val], na.rm = TRUE)
    }
    
    return(VF_summary)
}

# Create VF summary using the improved function
VF_summary <- create_VF_summary_from_data()

if (!is.null(VF_summary)) {
    # Save VF summary
    output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    write.csv(VF_summary, paste0(output_dir, "/VF_summary_by_a11.csv"), row.names = FALSE)
    cat("Saved VF summary to VF_summary_by_a11.csv\n")
    
    # Print summary
    cat("VF Matrix Elements Summary:\n")
    print(VF_summary)
}

cat("\n=== Sensitivity analysis plotting completed ===\n")
cat("Plots saved in:", paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name, "/"), "\n")
cat("Generated plots:\n")
cat("- f and mu parameters: sensitivity_violin_plot_f_mu_params.png\n")
cat("- VY and other parameters: sensitivity_violin_plot_VY_other_params.png\n") 
cat("- Remaining parameters: sensitivity_violin_plot_remaining_params.png\n")
cat("- VF matrix elements: sensitivity_violin_plot_VF_matrix.png\n")
cat("Red dotted line: True a11 value\n")
cat("VF matrix calculated from: VF = 2*f*VY*t(f) + f*VY*mu*VY*t(f) + f*VY*t(mu)*VY*t(f)\n")