exp_latent30 <- readRDS("Data/Paper/Expected/Model_latent30_finalGen.rds")
exp_latent50 <- readRDS("Data/Paper/Expected/Model_latent50_finalGen.rds")
exp_latent70 <- readRDS("Data/Paper/Expected/Model_latent70_finalGen.rds")
exp_latent90 <- readRDS("Data/Paper/Expected/Model_latent90_finalGen.rds")

calculate_means <- function(simulations) {
  # Get the names of the parameters
  param_names <- names(simulations[[1]])
  
  # Initialize a list to store the means
  mean_params <- list()
  
  # Loop through each parameter and calculate the element-wise mean
  for (param in param_names) {
    # Extract the parameter matrices from all simulations
    param_matrices <- lapply(simulations, function(sim) sim[[param]])
    
    # Calculate the element-wise mean
    mean_matrix <- Reduce("+", param_matrices) / length(param_matrices)
    
    # Store the mean matrix in the list
    mean_params[[param]] <- mean_matrix
  }
  
  return(mean_params)
}



break_into_named_vector <- function(mean_params) {
  named_vector <- c()
  
  for (param in names(mean_params)) {
    matrix <- mean_params[[param]]
    for (i in 1:nrow(matrix)) {
      for (j in 1:ncol(matrix)) {
        new_name <- paste0(param, i, j)
        named_vector[new_name] <- matrix[i, j]
      }
    }
  }
  
  return(named_vector)
}

# Calculate the means
means_latent30 <- calculate_means(exp_latent30)
means_latent50 <- calculate_means(exp_latent50)
means_latent70 <- calculate_means(exp_latent70)
means_latent90 <- calculate_means(exp_latent90)

v_means_latent30 <- break_into_named_vector(means_latent30)
v_means_latent50 <- break_into_named_vector(means_latent50)
v_means_latent70 <- break_into_named_vector(means_latent70)
v_means_latent90 <- break_into_named_vector(means_latent90)

v_means_latent30 <- c(v_means_latent30, 
                      "a11" = 0.438178,   "a22" = 0.464758, "delta11" = 0.669328, "delta22" = 0.3794733, 
                      "f11" = 0.15, "f12" = 0.1, "f21" = 0.05, "f22" = 0.1, "k12" = .05)
v_means_latent50 <- c(v_means_latent50, "a11" = 0.5656854,   "a22" = 0.464758, "delta11" = 0.5656854, "delta22" = 0.3794733, 
                      "f11" = 0.15, "f12" = 0.1, "f21" = 0.05, "f22" = 0.1, "k12" = .05)
v_means_latent70 <- c(v_means_latent70, "a11" = 0.669328,   "a22" = 0.464758, "delta11" = 0.438178, "delta22" = 0.3794733, 
                      "f11" = 0.15, "f12" = 0.1, "f21" = 0.05, "f22" = 0.1, "k12" = .05)
v_means_latent90 <- c(v_means_latent90, "a11" = 0.7589466,   "a22" = 0.464758, "delta11" = 0.2529822, "delta22" = 0.3794733, 
                      "f11" = 0.15, "f12" = 0.1, "f21" = 0.05, "f22" = 0.1, "k12" = .05)

# save the named list into a txt with two columns
write.table(as.data.frame(v_means_latent30), file = "Data/Paper/Expected/Model_latent30_finalGen.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(v_means_latent50), file = "Data/Paper/Expected/Model_latent50_finalGen.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(v_means_latent70), file = "Data/Paper/Expected/Model_latent70_finalGen.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(v_means_latent90), file = "Data/Paper/Expected/Model_latent90_finalGen.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)