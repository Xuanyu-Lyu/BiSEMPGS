# this script is designed to fit the BiSEMPGS model with the txt data files.
# OpenMx summary will be saved  in a list and then saved as a .rds file in the analysis subfolder

# the path to the data
source("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/BiSEMPGS_fitOpenMx.R")
conditionNames <- c("Full_Model", "MeasurePgs30", "MeasurePgs10", "MeasurePgsFully", 
"f11-decrease", "f12-decrease", "f11.12.21.22-decrease", "am11-decrease", "am12-decrease", "am11.12.21.22-decrease", "Full_Model_.5latent")
# the path to the dir with all txt files
data_path <- c( paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/", conditionNames[11], "/nfam16000"),
                paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/", conditionNames[11], "/nfam48000"),
                paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/", conditionNames[11], "/nfam32000"),
                paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/", conditionNames[11], "/nfam64000"))
# the path to the dir to save the results
save_path <- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Analysis/", conditionNames[11])

data_pattern <- c("_16000.txt", "_48000.txt", "_32000.txt", "_64000.txt")
save_pattern <- c("_16000","_48000", "_32000", "_64000")
model_type <- "m2"
mxSetup <- "_.05lb_smallerTol_newSetup_adelta"
for (j in 1:4){
    # a list to save all the summary data
    summary_list <- list()

    # get all the names of the text files that end with 32000.txt
    l_files <- list.files(data_path[j], pattern = data_pattern[j])

    # fit the model for each data file
    for (i in 1:length(l_files)){
        # fit the model
        fit <- fitBiSEMPGS_m2(paste0(data_path[j], "/", l_files[i]))
        # save the summary
        summary_list[[l_files[i]]] <- fit
        # save the fit
        cat("Model", l_files[i], "has been fitted\n")
    }
    if (!dir.exists(save_path)){
        dir.create(save_path)}

    # save the summary list
    saveRDS(summary_list, paste0(save_path, "/", model_type,mxSetup, save_pattern[j], "_summary_list.rds"))
}
