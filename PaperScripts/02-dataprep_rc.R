# a script to convert the .rds simulated data into .txt files for the BiSEMPGS model
# 12/10/24 update: dont need to run 03-createMoreTxt.R as we simulated the data with the target_n 

# the path to the data
# conditionNames <- c("Full_Model", "MeasurePgs30", "MeasurePgs10", "MeasurePgsFully", 
# "f11-decrease", "f12-decrease", "f11.12.21.22-decrease", "am11-decrease", "am12-decrease", "am11.12.21.22-decrease", "Full_Model_test")
#conditionNames <- c("Model_latent30", "Model_latent50", "Model_latent70", "Model_latent90")
conditionNames <- c("Model_r2_16", "Model_r2_8", "Model_r2_4", "Model_r2_2", "Model_r2_1")

# a data prep function

data_prep <- function(data_path, save_path, target_n){
    # get the number of .rds files in the folder
    l_files <- list.files(data_path, pattern = ".rds")
    # read all the .rds file one by one and save them as .txt files
    for (i in 1:length(l_files)){
        
        # check if the data has been converted to txt
        if (file.exists(paste0(save_path, "/", l_files[i],"_",target_n, ".txt"))){
            cat("Data file", l_files[i], "has been converted to txt\n")
            next
        }
        # load the data
        data <- readRDS(paste0(data_path, "/", l_files[i]))

        data_df <- data$PHEN

        data_df <- data_df[,c("ID", "Father.ID", "Mother.ID",
                              "Y1P","Y2P","Y1M","Y2M","Y1","Y2",
                              "TPO1","TPO2","NTPO1","NTPO2",
                              "TMO1","TMO2","NTMO1","NTMO2")] |> as.data.frame()
        
        # remove the rows with the same Father.ID
        data_df <- data_df[!duplicated(data_df$Father.ID),]
        print(nrow(data_df))
        # remove the columns with ID
        data_df <- data_df[,-(1:3)]
        
        # change the column names to what we defined in the OpenMx script
        colnames(data_df) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")

        # sample the data to the target_n
        data_df <- data_df[sample(1:nrow(data_df), target_n),]
        # write the data as a txt file
        if (!dir.exists(save_path)){
		    dir.create(save_path)
	    }
        write.table(data_df, file = paste0(save_path, "/", l_files[i],"_",target_n, ".txt"), sep = "\t", row.names = FALSE)
        cat("Data file", l_files[i], "has been saved as a txt file\n")
        #return(data_df)
    }
}

# run the function
for (i in 1:5){
    data_path <- paste0("/scratch/alpine/xuly4739/BiSEMPGS/Data/", conditionNames[i])
    v_sample <- c(.4e4,.8e4 ,1.6e4, 3.2e4, 4.8e4, 6.4e4)
    for(j in 1:length(v_sample)){
        save_path <- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/", conditionNames[i], "/nfam", v_sample[j])
        data_prep(data_path, save_path, v_sample[j])
         }
    }


# data_path <- paste0("/scratch/alpine/xuly4739/BiSEMPGS/Data/", conditionNames[11])
# save_path <- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/", conditionNames[11], "/nfam16000")
# data_prep(data_path, save_path, 1.6e4)

