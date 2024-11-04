# This script is designed for merging the 16k txt files into larger samples

# the path to the txt files
# conditionNames <- c("Full_Model", "MeasurePgs30", "MeasurePgs10", "MeasurePgsFully", 
# "f11-decrease", "f12-decrease", "f11.12.21.22-decrease", "am11-decrease", "am12-decrease", "am11.12.21.22-decrease", "Full_Model_test")

conditionNames <- c("Model_latent30", "Model_latent50", "Model_latent70", "Model_latent90")




# a function to merge the txt files
merge_txt <- function(data_path, save_path, merge_count = 2){
    l_files <- list.files(data_path, pattern = "16000.txt")
    finalFileCount <- floor(length(l_files) / merge_count)
    for (i in 1:finalFileCount){
        # check if the data has been created
        if (file.exists(paste0(save_path, "/", sub("16000.txt$", "", l_files[(i-1)*merge_count + 1]), merge_count*16000, ".txt"))){
            cat("Data file", paste0(save_path, "/", sub("16000.txt$", "", l_files[(i-1)*merge_count + 1]), merge_count*16000, ".txt"), "has been created\n")
            next
        }
        # read the first file
        data_df <- read.table(paste0(data_path, "/", l_files[(i-1)*merge_count + 1]), header = TRUE)
        for (j in 2:merge_count){
            data_df <- rbind(data_df, read.table(paste0(data_path, "/", l_files[(i-1)*merge_count + j]), header = TRUE))
        }
        # write the data as a txt file
        if (!dir.exists(save_path)){
		    dir.create(save_path)
	    }
        write.table(data_df, file = paste0(save_path, "/", sub("16000.txt$", "", l_files[(i-1)*merge_count + 1]), merge_count*16000, ".txt"), sep = "\t", row.names = FALSE)
        cat("Data file", paste0(save_path, "/", sub("16000.txt$", "", l_files[(i-1)*merge_count + 1]), merge_count*16000, ".txt"), "has been saved\n")
    }
} 

# # run the function
# for (i in 1:length(data_path)){
#     data_path <- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/", conditionNames[i], "/nfam16000")
#     save_path <- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/", conditionNames[i], "/nfam32000")
#     save_path2<- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/", conditionNames[i], "/nfam48000")
#     save_path3<- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/", conditionNames[i], "/nfam64000")
#     save_path4<- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/", conditionNames[i], "/nfam80000")
#     merge_txt(data_path = data_path, save_path = save_path, 2)
#     merge_txt(data_path = data_path, save_path = save_path2, 3)
#     merge_txt(data_path = data_path, save_path = save_path3, 4)
#     merge_txt(data_path = data_path, save_path = save_path4, 5)
# }

# run the function
for (i in 1:length(conditionNames)){
    data_path <- paste0("Data/Paper/", conditionNames[i], "/nfam16000")
    save_path <- paste0("Data/Paper/", conditionNames[i], "/nfam32000")
    save_path2<- paste0("Data/Paper/", conditionNames[i], "/nfam48000")
    save_path3<- paste0("Data/Paper/", conditionNames[i], "/nfam64000")
    save_path4<- paste0("Data/Paper/", conditionNames[i], "/nfam80000")
    merge_txt(data_path = data_path, save_path = save_path, 2)
    merge_txt(data_path = data_path, save_path = save_path2, 3)
    merge_txt(data_path = data_path, save_path = save_path3, 4)
    merge_txt(data_path = data_path, save_path = save_path4, 5)
}

