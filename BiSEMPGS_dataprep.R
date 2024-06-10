# this is a script to prepare the data for the BiSEMPGS model from the simulated data

# the path to the data
data_path <- "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/testdata/loop1.rds"
save_path <- "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/testdata"

# load the data
data <- readRDS(data_path)
data_df <- data$PHEN
colnames(data_df)
data_df <- data_df[,c("ID", "Father.ID", "Mother.ID",
                      "Y1P","Y2P","Y1M","Y2M","Y1","Y2",
                      "TPO1","TPO2","NTPO1","NTPO2",
                      "TMO1","TMO2","NTMO1","NTMO2")] |> as.data.frame()
nrow(data_df)

summary(table(data_df$Father.ID) |> as.factor())
# if more than one rows has the same Father.ID, remove until only one row is left


# remove the rows with the same Father.ID
data_df <- data_df[!duplicated(data_df$Father.ID),]
nrow(data_df)
# remove the columns with ID
data_df <- data_df[,-(1:3)]

# change the column names to what we defined in the OpenMx script
colnames(data_df) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
# write the data as a txt file
write.table(data_df, file = paste0(save_path, "/loop1.txt"), sep = "\t", row.names = FALSE)
