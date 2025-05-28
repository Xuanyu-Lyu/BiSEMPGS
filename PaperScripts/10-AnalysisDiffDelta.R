# read the file

testData <- read.table("Data/Paper/differentDelta/OnlyUni/loop1.rds_12000.txt", header = TRUE)

# only keep trait 1 columns
# subset the data to only include varibles whose names end with 1
testData <- testData[, grepl("1$", names(testData))]
# arrange the order of columns to match the OpenMx script
testData <- testData[, c("NTm1", "Tm1", "NTp1", "Tp1", "Ym1", "Yp1", "Yo1")]
cov(testData)

source("PaperScripts/08-uniModelBias_OpenMx.R")

test_summary <- fitUniSEMPGS_diffDelta(testData)

test_summary
