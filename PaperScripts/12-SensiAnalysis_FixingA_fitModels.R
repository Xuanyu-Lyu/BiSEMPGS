# Sensitivity analysis for BiSEMPGS model with fixed A parameter at different levels
# Author: Generated for sensitivity analysis
# Date: August 4, 2025
# This script fits the BiSEMPGS model with parameter A fixed at different levels
# Using already simulated MVN data for r2pgs04 condition with 32k sample size

# Load required functions
source("PaperScripts/04-OpenMxFunctions.R")

# Define the condition and sample size (only r2pgs04 with 32k as requested)
condition_name <- "r2pgs04"
sample_size_name <- "32k"

# Calculate the true A values for r2pgs04 condition (condition 3)
# vg1 = 0.64, vg2 = 0.36, prop.h2.latent1 = 0.60/0.64, prop.h2.latent2 = 0.8
true_a11 <- sqrt(0.64 * (0.60/0.64))  # sqrt(0.6) ≈ 0.7746
true_a22 <- sqrt(0.36 * 0.8)          # sqrt(0.288) ≈ 0.5367

# Define different fixed A values for sensitivity analysis (10 levels around true values)
# Only varying a11 (trait 1), keeping a22 (trait 2) fixed at true value
# Using steps of 0.025 to create 10 levels from -0.05 to +0.05 around true a11
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

a_value_names <- names(a_values_list)

cat("Starting sensitivity analysis for fixed A parameter\n")
cat("Condition:", condition_name, "\n")
cat("Sample size:", sample_size_name, "\n")
cat("True A values: a11 =", round(true_a11, 4), ", a22 =", round(true_a22, 4), "\n")
cat("A value levels to test:", paste(a_value_names, collapse = ", "), "\n")
cat("Detailed A values:\n")
for(i in seq_along(a_values_list)) {
    cat("  ", names(a_values_list)[i], ": a11 =", round(a_values_list[[i]][1], 4), 
        ", a22 =", round(a_values_list[[i]][2], 4), "\n")
}

# Check if data directory exists
data_dir <- paste0("Data/Paper/test_v12_ss/", condition_name, "_", sample_size_name)
if (!dir.exists(data_dir)) {
    stop("Data directory does not exist: ", data_dir)
}

# List available data files
data_files <- list.files(data_dir, pattern = "samples_iter_empi_.*\\.tsv", full.names = FALSE)
if (length(data_files) == 0) {
    stop("No data files found in: ", data_dir)
}

cat("Found", length(data_files), "data files\n")

# Extract iteration numbers from filenames and select first 100
iter_numbers <- as.numeric(gsub("samples_iter_empi_(\\d+)\\.tsv", "\\1", data_files))
iter_numbers <- sort(iter_numbers)
iter_numbers_to_use <- iter_numbers[1:min(100, length(iter_numbers))]

cat("Using", length(iter_numbers_to_use), "iterations for analysis\n")

# Loop through each A value level
for (a_idx in seq_along(a_values_list)) {
    a_level_name <- a_value_names[a_idx]
    a_values <- a_values_list[[a_idx]]
    
    cat("\n=== Processing A level:", a_level_name, "===\n")
    cat("A values: a11 =", a_values[1], ", a22 =", a_values[2], "\n")
    
    # Initialize summary list for this A level
    summary_list_fixedA <- list()
    
    # Loop through selected iterations
    for (i in seq_along(iter_numbers_to_use)) {
        iter_num <- iter_numbers_to_use[i]
        
        # Construct data file path
        data_file_path <- paste0(data_dir, "/samples_iter_empi_", iter_num, ".tsv")
        
        if (!file.exists(data_file_path)) {
            cat("Warning: Data file not found:", data_file_path, "\n")
            next
        }
        
        # Fit the model with fixed A values
        tryCatch({
            fit_fixedA <- fitBiSEMPGS_m2_tol_fixH2(
                data_path = data_file_path,
                avalue = a_values,
                feaTol = 1e-6, 
                optTol = 1e-9,
                jitterMean = 0.5,
                jitterVar = 0.1,
                exhaustive = FALSE,
                extraTries = 5
            )
            
            summary_list_fixedA[[paste0("iter_", iter_num)]] <- fit_fixedA
            cat("Completed iteration", i, "/", length(iter_numbers_to_use), "(iter_num:", iter_num, ")\n")
            
        }, error = function(e) {
            cat("Error in iteration", iter_num, ":", e$message, "\n")
            summary_list_fixedA[[paste0("iter_", iter_num)]] <- list(error = e$message)
        })
    }
    
    # Save results for this A level
    output_dir <- paste0("Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    save_path <- paste0(output_dir, "/", a_level_name, "_fixedA_summary_list_100.rds")
    saveRDS(summary_list_fixedA, save_path)
    
    cat("Saved results for A level", a_level_name, "to:", save_path, "\n")
    cat("Total models fitted:", length(summary_list_fixedA), "\n")
}

cat("\n=== Sensitivity analysis completed ===\n")
cat("Results saved in: Analysis/Paper/SensiAnalysis_FixedA/", condition_name, "_", sample_size_name, "/\n")