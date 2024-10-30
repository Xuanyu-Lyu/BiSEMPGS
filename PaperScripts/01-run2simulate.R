# this script is used to simulate the data for the BiSEMPGS model to test the math

# load the Gene Evolve simulation functions
source("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Simulate.Multivariate.NOLD.AM.FUNCTIONS_May2021-mck7.R")
#library(expm)
#setwd("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS")

# fetch the job array information
args <- commandArgs(trailingOnly = TRUE)
array_idx <- as.numeric(args[1])

# saving directory
save_dir_summary <- "/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS"
save_dir_data <- "/scratch/alpine/xuly4739/BiSEMPGS"

# # create a list of starting parameters for different conditions
# conditionNames <- c("Full_Model", "MeasurePgs30", "MeasurePgs10", "MeasurePgsFully", 
# "f11-decrease", "f12-decrease", "f11.12.21.22-decrease", "am11-decrease", "am12-decrease", "am11.12.21.22-decrease", "Full_Model_test")
# # starting parameter list 1 is for the full model with different sample sizes
# startingParamList1 <- list(vg1 = rep(.49,11),
# 						   vg2 = rep(.16,11),
# 						   am11 = c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0,0.4,0,.4),
# 						   am12 = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0,0,.2),
# 						   am21 = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0,.1),
# 						   am22 = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0,.3),
# 						   f11 = c(0.15,0.15,0.15,0.15,0.075,0.15,0.075,0.15,0.15,0.15,.15),
# 						   f12 = c(0.1,0.1,0.1,0.1,0.1,0.05,0.05,0.1,0.1,0.1,.1),
# 						   f21 = c(0.05,0.05,0.05,0.05,0.05,0.05,0.025,0.05,0.05,0.05,.05),
# 						   f22 = c(0.1,0.1,0.1,0.1,0.1,0.1,0.05,0.1,0.1,0.1,.1),
# 						   Nfam = c(5e4, rep(5e4,10)),
# 						   rg = rep(.1,11),
# 						   re = rep(.1,11),
# 						   prop.h2.latent1 = c(0.5,0.7,0.9,0,rep(0.5,6),.3),
# 						   prop.h2.latent2 = c(0.7,0.7,0.9,0,rep(0.7,6),.7))
# # condition 11 a11 = 0.3834058 a22 = 0.334664
# sqrt(.16*.7)
# new starting parameter list for the paper
startingParamList1 <- list(vg1 = rep(.64,4),
						   vg2 = rep(.36,4),
						   am11 = rep(0.4,4),
						   am12 = rep(0.2,4),
						   am21 = rep(0.1,4),
						   am22 = rep(0.3,4),
						   f11 = rep(0.15,4),
						   f12 = rep(0.1,4),
						   f21 = rep(0.05,4),
						   f22 = rep(0.1,4),
						   Nfam = rep(5e4, 4),
						   rg = rep(.1,4),
						   re = rep(.1,4),
						   prop.h2.latent1 = c(.3,.5,.7,.9),
						   prop.h2.latent2 = c(.6,.6,.6,.6))

conditionNames <- c("Model_latent30", "Model_latent50", "Model_latent70", "Model_latent90")


for (condition in 1:4){
	# WILDCARD parameters
	pop.size <-  startingParamList1["Nfam"][[1]][[condition]] #maybe something like 2e4, or 20000, when running for real
	num.cvs <- 200 #maybe 25
	seed <- 62*condition
	num.gen <-  20 #8 should be sufficient to get to AM equilibrium
	avoid.inb <- TRUE #avoid inbreeding?
	save.covariances <- TRUE #save the covariance matrices?
	save.history <- TRUE #save the data for each generation?
	min.maf <- .1 #.2 for real
	max.maf <- .5 #.5

	#USER INPUT VARIABLES
	#VG
	vg1 <- startingParamList1["vg1"][[1]][[condition]] #trait 1 vg
	vg2 <- startingParamList1["vg2"][[1]][[condition]] #trait 2 vg
	rg <- startingParamList1["rg"][[1]][[condition]] #genetic CORRELATION @t0 bw trait 1 and trait 2 for both obs. PGS and latent PGS (assumed to be the same). NOTE: this is NOT the full rg at t0. It is the rg bw PGSs, and the rg bw LGSs. The full rg may be a bit different (Simpson's paradox)
	(k2.matrix <- matrix(c(1,rg,rg,1),nrow=2,byrow=T)) #k2 matrix is 2 * k matrix - i.e., genotypic (instead of haplotypic) var/covar at t0
	prop.h2.latent1 <- startingParamList1["prop.h2.latent1"][[1]][[condition]] #  trait 1, e.g., height
	prop.h2.latent2 <- startingParamList1["prop.h2.latent2"][[1]][[condition]] # trait 2, e.g., IQ

	#AM - these are NOT the mu copaths. They are the CORRELATIONS between male & female traits
	am11 <-  startingParamList1["am11"][[1]][[condition]] #height.m-height.f am across 2 its
	am12 <-  startingParamList1["am12"][[1]][[condition]] #height.m-iq.f; e.g., tall males choose smart females
	am21 <-  startingParamList1["am21"][[1]][[condition]] #iq.m-height.f
	am22 <-  startingParamList1["am22"][[1]][[condition]] #iq.m - iq.f

	#VT
	f11 <- startingParamList1["f11"][[1]][[condition]] # regression of offspring trait 1 F on parental trait 1
	f12 <- startingParamList1["f12"][[1]][[condition]] # regression of offspring trait 1 F on parental trait 2
	f21 <- startingParamList1["f21"][[1]][[condition]] # regression of offspring trait 2 F on parental trait 1
	f22 <- startingParamList1["f22"][[1]][[condition]] # regression of offspring trait 2 F on parental trait 2

	#E
	re <- startingParamList1["re"][[1]][[condition]] #environmental CORRELATION between trait 1 & trait 2

	# IMPLIED variables
	#This section just converts the above inputs into matrices in our algebra

	#Create lower delta matrix (d stands for delta)
	(vg.obs1 <- vg1*(1-prop.h2.latent1))
	(vg.obs2 <- vg2*(1-prop.h2.latent2))
	(covg.obs <- rg*sqrt(vg.obs1)*sqrt(vg.obs2))
	d11 <- sqrt(vg.obs1)
	#d21 <- covg.obs/d11
	d21 <- 0 #new way
	d22 <- sqrt(vg.obs2 - d21^2)
	(delta.mat <- matrix(c(d11,0,d21,d22),byrow=T,nrow=2))
	(vg.obs <- delta.mat %*% k2.matrix %*% t(delta.mat)) #Check
	matrix(c(vg.obs1,covg.obs,covg.obs,vg.obs2),nrow=2) #Check

	#Create lower a matrix
	(vg.lat1 <- vg1*(prop.h2.latent1))
	(vg.lat2 <- vg2*(prop.h2.latent2))
	(covg.lat <- rg*sqrt(vg.lat1)*sqrt(vg.lat2))
	a11 <- sqrt(vg.lat1)
	#a21 <- covg.lat/a11
	a21 <- 0 #new way
	a22 <- sqrt(vg.lat2 - a21^2)
	(a.mat <- matrix(c(a11,0,a21,a22),byrow=T,nrow=2))
	(vg.lat <- a.mat %*% k2.matrix %*% t(a.mat)) #Check
	matrix(c(vg.lat1,covg.lat,covg.lat,vg.lat2),nrow=2) #Check

	#Find the total genetic covariance
	#(covg <- vg.obs[1,2] + vg.lat[1,2])
	(covg.mat <- vg.obs + vg.lat)
	a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) #check

	#Create E lower matrix
	ve1 <- 1 - vg1
	ve2 <- 1 - vg2
	(cove <- re*sqrt(ve1*ve2))
	(cove.mat <- matrix(c(ve1,cove,cove,ve2),nrow=2,byrow=T))

	#Create var-covar(Y) @t0. NOTE: we define the base population as the population before any AM or VT has occurred.
	(COVY <- covg.mat+cove.mat)
	a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) + cove.mat  #Check; these two should have 1 on the diagonals and the off-diags should be exactly the same if we've done everything above correctly.
	#Note that our actual model estimates k and j separately, and thus our model does not make the assumption made in this simulation that the observed and latent rg at t0 are the same (or does it? do we actually estimate j12, or do we set it equal to k12?)


	#Note: we will add influences of AM and VT AFTER time 0
	#AM - make sure to change am.list if you want AM to change across generations
	(mate.cor.mat <- matrix(c(am11,am12,am21,am22),nrow=2,byrow=T))
	am.list <- rep(list(mate.cor.mat),num.gen+1)

	#VF
	(f.mat <- matrix(c(f11,f12,f21,f22),nrow=2,byrow=T))
	(covf.mat <- 2 * (f.mat %*% COVY %*% t(f.mat)))

	maf.vector <- runif(num.cvs,min.maf,max.maf)  #Can change the distribution of MAFs here
	gentp.var <- maf.vector*(1-maf.vector)*2

	alphas.pre <- mvrnorm(num.cvs,c(0,0),matrix(c(1,rg,rg,1),nrow=2),empirical=TRUE)
	#CAREFUL: when you specify empirical=TRUE, sum(alphas.pre^2) is NOT num.cvs. It is num.cvs-1
	colSums(alphas.pre^2) #should be num.cvs-1. Compensate for that below using num.cvs-1
	alphas <- alphas.pre * cbind(sqrt(1/((num.cvs-1)*gentp.var)),sqrt(1/((num.cvs-1)*gentp.var)))

	cv.info <- data.frame(maf=maf.vector,alpha1=alphas[,1],alpha2=alphas[,2]) #we'll use this for both the observed and latent
	ObjectsKeep <- as.character(ls())
	# write a loop to run the simulation 100 times and save all the summary data in a list
	#l.summaryLast <- list()
	#l.all <- list()
	for (i in 1:10){
	loop_index <- (array_idx-1)*10 + i 
	# if the file is already there, skip
	if (file.exists(paste0(save_dir_data,"/Data/",conditionNames[condition],"/loop",loop_index,".rds"))){
		cat(conditionNames[condition],"/Simulation",loop_index,"already exists\n")
		next
	}
	AM.DATA <- AM.SIMULATE(CV.INFO=cv.info, NUM.GENERATIONS=num.gen, POP.SIZE=pop.size, AVOID.INB=avoid.inb, SAVE.EACH.GEN=save.history, SAVE.COVS=save.covariances, SEED=seed*loop_index, 
						cove.mat=cove.mat, fmat=f.mat, amat=a.mat, dmat=delta.mat, cor.list=am.list, covy=COVY, k2.matrix=k2.matrix)
	SUMMARY.last <- AM.DATA$SUMMARY.RES[[num.gen]]
	#l.summaryLast[[i]] <- SUMMARY.last
	#l.all[[i]] <- AM.DATA
	# test if a folder exist, if not, create one
	if (!dir.exists(paste0(save_dir_summary, "/Summary/",conditionNames[condition]))){
		dir.create(paste0(save_dir_summary, "/Summary/",conditionNames[condition]))
	}		
	if (!dir.exists(paste0(save_dir_data,"/Data/",conditionNames[condition]))){
		dir.create(paste0(save_dir_data, "/Data/",conditionNames[condition]))
	}
	# save the data
	saveRDS(SUMMARY.last, file=paste0(save_dir_summary,"/Summary/",conditionNames[condition],"/loop",loop_index,".rds"))
	saveRDS(AM.DATA, file=paste0(save_dir_data,"/Data/",conditionNames[condition],"/loop",loop_index,".rds"))
	cat(conditionNames[condition],"/Simulation",loop_index,"done\n")
	rm(list = setdiff(ls(), c(ObjectsKeep, "i", "ObjectsKeep")))
}


}


