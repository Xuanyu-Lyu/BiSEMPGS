# this script is used to simulate the data for the BiSEMPGS model to test the math

source("Simulate.Multivariate.NOLD.AM.FUNCTIONS_May2021-mck7.R")
library(expm)
#setwd("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS")

# create a list of starting parameters for different conditions
conditionNames <- c("f_trans_full", "am_trans_full", "delta_trans_full", "f_half_full", "am_half_full", "delta_half_full", "normal_full")
startingParamList <- list(vg1 = c(.36,.36,.09,.36,.36,.18,.36),
						  vg2 = c(.09,.09,.36,.09,.09,.045,.09),
						  am11 = c(.4,.3,.4,.4,.2,.4,.4),
						  am12 = c(.2,.2,.2,.2,.1,.2,.2),
						  am21 = c(.2,.2,.2,.2,.1,.2,.2),
						  am22 = c(.3,.4,.3,.3,.15,.3,.3),
						  f11 = c(.1,.15,.15,.15/2,.15,.15,.15),
						  f12 = c(.1,.1,.1,.05,.1,.1,.1),
						  f21 = c(.05,.05,.05,.025,.05,.05,.05),
						  f22 = c(.15,.1,.1,.05,.1,.1,.1))

for (condition in 7){
	# WILDCARD parameters
	pop.size <- 20000 #maybe something like 2e4, or 20000, when running for real
	num.cvs <- 200 #maybe 25
	seed <- 62*condition
	num.gen <-  15 #8 should be sufficient to get to AM equilibrium
	avoid.inb <- TRUE #avoid inbreeding?
	save.covariances <- TRUE #save the covariance matrices?
	save.history <- TRUE #save the data for each generation?
	min.maf <- .1 #.2 for real
	max.maf <- .5 #.5

	#USER INPUT VARIABLES
	#VG
	vg1 <- startingParamList["vg1"][[1]][[condition]] #trait 1 vg
	vg2 <- startingParamList["vg2"][[1]][[condition]] #trait 2 vg
	rg <- 0 #genetic CORRELATION @t0 bw trait 1 and trait 2 for both obs. PGS and latent PGS (assumed to be the same). NOTE: this is NOT the full rg at t0. It is the rg bw PGSs, and the rg bw LGSs. The full rg may be a bit different (Simpson's paradox)
	(k2.matrix <- matrix(c(1,rg,rg,1),nrow=2,byrow=T)) #k2 matrix is 2 * k matrix - i.e., genotypic (instead of haplotypic) var/covar at t0
	prop.h2.latent1 <- 0 #  trait 1, e.g., height
	prop.h2.latent2 <- 0 # trait 2, e.g., IQ

	#AM - these are NOT the mu copaths. They are the CORRELATIONS between male & female traits
	am11 <-  startingParamList["am11"][[1]][[condition]] #height.m-height.f am across 2 its
	am12 <-  startingParamList["am12"][[1]][[condition]] #height.m-iq.f; e.g., tall males choose smart females
	am21 <-  startingParamList["am21"][[1]][[condition]] #iq.m-height.f
	am22 <-  startingParamList["am22"][[1]][[condition]] #iq.m - iq.f

	#VT
	f11 <- startingParamList["f11"][[1]][[condition]] # regression of offspring trait 1 F on parental trait 1
	f12 <- startingParamList["f12"][[1]][[condition]] # regression of offspring trait 1 F on parental trait 2
	f21 <- startingParamList["f21"][[1]][[condition]] # regression of offspring trait 2 F on parental trait 1
	f22 <- startingParamList["f22"][[1]][[condition]] # regression of offspring trait 2 F on parental trait 2

	#E
	re <- 0 #environmental CORRELATION between trait 1 & trait 2

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
	for (i in 101:300){
		AM.DATA <- AM.SIMULATE(CV.INFO=cv.info, NUM.GENERATIONS=15, POP.SIZE=pop.size, AVOID.INB=avoid.inb, SAVE.EACH.GEN=save.history, SAVE.COVS=save.covariances, SEED=seed*i, 
							cove.mat=cove.mat, fmat=f.mat, amat=a.mat, dmat=delta.mat, cor.list=am.list, covy=COVY, k2.matrix=k2.matrix)
		SUMMARY.last <- AM.DATA$SUMMARY.RES[[15]]
		#l.summaryLast[[i]] <- SUMMARY.last
		#l.all[[i]] <- AM.DATA
		# test if a folder exist, if not, create one
		if (!dir.exists(paste0("Summary/",conditionNames[condition]))){
			dir.create(paste0("Summary/",conditionNames[condition]))
		}		
		if (!dir.exists(paste0("Data/",conditionNames[condition]))){
			dir.create(paste0("Data/",conditionNames[condition]))
		}
		saveRDS(SUMMARY.last, file=paste0("Summary/",conditionNames[condition],"/loop",i,".rds"))
		saveRDS(AM.DATA, file=paste0("Data/",conditionNames[condition],"/loop",i,".rds"))
		cat(conditionNames[condition],"Simulation",i,"done\n")
		rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep")))
	}
	# save the summary data into a rds file
	
}


