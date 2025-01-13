# test the openmx VE22

test_df <- read.table("Data/Paper/Model_r2_8/nfam48000/loop2.rds_48000.txt", header = TRUE)
source("PaperScripts/04-OpenMxFunctions.R")
model1 <- fitBiSEMPGS_m2_tol("Data/Paper/Model_r2_8/nfam48000/loop2.rds_48000.txt")

summary(model1)
mxEval(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra + VE, model1)
mxEval(VE, model1)

mxEval(VY, model1)
mxEval(VY_Algebra, model1)
