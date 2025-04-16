# This script is for checking the bias of univariate model when the population effects are bivariate using MVN simulation
conditionNames <- c("onlyAM", "onlyVT", "bothAMVT", "bothAMVT_StrongCrossTrait")
source("PaperScripts/04-OpenMxFunctions.R")
startingParamList1 <- list(vg1 = rep(.64,4),
						   vg2 = rep(.36,4),
						   am11 = c(rep(0.4,3),.05),
						   am12 = c(.2,0,.2,.4),
						   am21 = c(.2,0,.2,.4), 
						   am22 = c(rep(0.3,3),.1),
						   f11 = c(rep(0.15,3),0),
						   f12 = c(0,.1,.1,.3),
						   f21 = c(0,.05,.05,.25), 
						   f22 = rep(0.1,4),
						   Nfam = rep(8e4, 4),
						   rg = rep(0,4),
						   re = rep(0,4),
						   prop.h2.latent1 = rep(.60/.64,4),
						   prop.h2.latent2 = rep(.8,4))

for (condition in 1:4){
			# WILDCARD parameters
		pop.size <-  startingParamList1["Nfam"][[1]][[condition]] #maybe something like 2e4, or 20000, when running for real

		#USER INPUT VARIABLES
		#VG
		vg1 <- startingParamList1["vg1"][[1]][[condition]] #trait 1 vg
		vg2 <- startingParamList1["vg2"][[1]][[condition]] #trait 2 vg
		rg <- startingParamList1["rg"][[1]][[condition]] #genetic CORRELATION @t0 bw trait 1 and trait 2 for both obs. PGS and latent PGS (assumed to be the same). NOTE: this is NOT the full rg at t0. It is the rg bw PGSs, and the rg bw LGSs. The full rg may be a bit different (Simpson's paradox)
		k2.matrix <- matrix(c(1,rg,rg,1),nrow=2,byrow=T) #k2 matrix is 2 * k matrix - i.e., genotypic (instead of haplotypic) var/covar at t0
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
		vg.obs1 <- vg1*(1-prop.h2.latent1)
		vg.obs2 <- vg2*(1-prop.h2.latent2)
		covg.obs <- rg*sqrt(vg.obs1)*sqrt(vg.obs2)
		d11 <- sqrt(vg.obs1)
		#d21 <- covg.obs/d11
		d21 <- 0 #new way
		d22 <- sqrt(vg.obs2 - d21^2)
		delta.mat <- matrix(c(d11,0,d21,d22),byrow=T,nrow=2)
		vg.obs <- delta.mat %*% k2.matrix %*% t(delta.mat) #Check


		#Create lower a matrix
		vg.lat1 <- vg1*(prop.h2.latent1)
		vg.lat2 <- vg2*(prop.h2.latent2)
		covg.lat <- rg*sqrt(vg.lat1)*sqrt(vg.lat2)
		a11 <- sqrt(vg.lat1)
		#a21 <- covg.lat/a11
		a21 <- 0 #new way
		a22 <- sqrt(vg.lat2 - a21^2)
		a.mat <- matrix(c(a11,0,a21,a22),byrow=T,nrow=2)
		vg.lat <- a.mat %*% k2.matrix %*% t(a.mat) #Check
		matrix(c(vg.lat1,covg.lat,covg.lat,vg.lat2),nrow=2) #Check

		#Find the total genetic covariance
		covg.mat <- vg.obs + vg.lat
		

		#Create E lower matrix
		ve1 <- 1 - vg1
		ve2 <- 1 - vg2
		cove <- re*sqrt(ve1*ve2)
		cove.mat <- matrix(c(ve1,cove,cove,ve2),nrow=2,byrow=T)

		#Create var-covar(Y) @t0. NOTE: we define the base population as the population before any AM or VT has occurred.
		COVY <- covg.mat+cove.mat

		#Note: we will add influences of AM and VT AFTER time 0
		#AM - make sure to change am.list if you want AM to change across generations
		mate.cor.mat <- matrix(c(am11,am12,am21,am22),nrow=2,byrow=T)

		#VF
		f.mat <- matrix(c(f11,f12,f21,f22),nrow=2,byrow=T)
		covf.mat <- 2 * (f.mat %*% COVY %*% t(f.mat))
		gens = 15
		a.t0 <- a.mat
		delta.t0 <- delta.mat
		j.t0 <- k2.matrix*.5
		k.t0 <- k2.matrix*.5
		f.t0 <- f.mat
		rmate.t0 <- mate.cor.mat
		covE.t0 <- cove.mat
		#rmate.t0 <- obs.cor

		#Define and initialize at t0 parameters that evolve
		exp.gc <- exp.hc <- exp.ic <- exp.gt <- exp.ht <- exp.itlo <- exp.itol <- exp.w <- exp.v <- exp.VY <- exp.VY.diagonal <- exp.corVY<- exp.mu <- exp.VF <- exp.Omega <- exp.Gamma <- exp.thetaT <- exp.thetaNT <- exp.thetaLT <- exp.thetaLNT <- mate.cov <- mate.cor <- exp.VGO <- exp.VGL <- exp.COVLO <- vector(mode="list",length=gens-1)
		#initialize parameters that begin as 0's in base population (before VT and AM)
		exp.VGO[[1]] <- exp.VGL[[1]] <- exp.COVLO[[1]] <- exp.gc[[1]] <- exp.hc[[1]] <- exp.ic[[1]] <- exp.w[[1]] <- exp.v[[1]] <- matrix(0,nrow=2,ncol=2)
		# define the true 4X4 covariance matrix for Vy at t0
		exp.covY.full <- exp.corY.full <- exp.covY.full.diagonal <- vector(mode="list",length=gens-1)
		#initialize parameters that are non-zero in base population (think of these as the first parents to pass on VT and engage in AM)
		exp.VY[[1]] <- COVY
		exp.VY.diagonal[[1]] <- diag(diag(exp.VY[[1]]))
		exp.VF[[1]] <- covf.mat
		exp.mu[[1]] <- solve(COVY) %*% rmate.t0 %*% solve(t(COVY))
		mate.cov[[1]] <- COVY %*% exp.mu[[1]] %*% COVY 
		exp.covY.full[[1]] <- rbind(cbind(exp.VY[[1]], t(mate.cov[[1]])),cbind(mate.cov[[1]], exp.VY[[1]]))
		exp.corY.full[[1]] <- cov2cor(exp.covY.full[[1]])
		mate.cor[[1]] <- mate.cov[[1]]/COVY
		mate.cor[[1]][1,2] <- mate.cov[[1]][1,2]/sqrt(COVY[1,1]*COVY[2,2])
		mate.cor[[1]][2,1] <- mate.cov[[1]][2,1]/sqrt(COVY[1,1]*COVY[2,2])

	for (it in 2:gens){
			(exp.Omega[[it]] <- 2*delta.t0%*%exp.gc[[it-1]] + delta.t0 %*% k.t0 + .5*exp.w[[it-1]] + 2*a.t0%*%exp.ic[[it-1]])
			#cat("Omega",it,"\n",exp.Omega[[it]],"\n")
			(exp.Gamma[[it]] <- 2*a.t0%*%exp.hc[[it-1]] + 2*delta.t0%*%t(exp.ic[[it-1]]) + a.t0%*%j.t0 + .5*exp.v[[it-1]])
			#cat("Gamma",it,"\n",exp.Gamma[[it]],"\n")

			(exp.VY[[it]] <- 2 * delta.t0 %*% t(exp.Omega[[it]]) + 2 * a.t0 %*% t(exp.Gamma[[it]]) + exp.w[[it-1]] %*% t(delta.t0) + exp.v[[it-1]] %*% t(a.t0) + exp.VF[[it-1]] + covE.t0 )
			#cat("VY","\n")
			#print(exp.VY[[it]])

			# two potential ways: fix mu or fix mate.cor
			exp.VY.diagonal[[it]] <- exp.VY[[it]]
			exp.VY.diagonal[[it]][1,2] <- 0
			exp.VY.diagonal[[it]][2,1] <- 0

			exp.covY.full.diagonal[[it]] <- cbind(rbind(exp.VY.diagonal[[it]],matrix(0,2,2)),rbind(matrix(0,2,2),exp.VY.diagonal[[it]]))

			exp.corVY[[it]] <- cov2cor(exp.VY[[it]])
			exp.corY.full[[it]] <- exp.corY.full[[1]]
			exp.corY.full[[it]][1,2] <- exp.corVY[[it]][1,2]
			exp.corY.full[[it]][2,1] <- exp.corVY[[it]][2,1]
			exp.corY.full[[it]][3,4] <- exp.corVY[[it]][1,2]
			exp.corY.full[[it]][4,3] <- exp.corVY[[it]][2,1]
			#exp.covY.full[[it]] <- cor2cov(exp.corY.full[[it]],exp.VY[[it]][1,1],exp.VY[[it]][2,2],exp.VY[[it]][1,1],exp.VY[[it]][2,2]) 
			exp.covY.full[[it]] <- sqrt(exp.covY.full.diagonal[[it]]) %*% exp.corY.full[[it]] %*% sqrt(exp.covY.full.diagonal[[it]])
			mate.cov[[it]] <- exp.covY.full[[it]][3:4,1:2] ## this is the mate.cov fixing mate.cor
			mate.cor[[it]] <- exp.corY.full[[1]][3:4,1:2]
			(exp.mu[[it]] <- solve(exp.VY[[it]]) %*% mate.cov[[it]] %*% solve(t(exp.VY[[it]]))) 
			#cat("mu",it,"\n",exp.mu[[it]],"\n")
			(exp.gt[[it]] <- t(exp.Omega[[it]])%*%exp.mu[[it]]%*%exp.Omega[[it]])
			#cat("gt",it,"\n",exp.gt[[it]],"\n")
			(exp.gc[[it]] <- .5*(exp.gt[[it]] + t(exp.gt[[it]])))
			#cat("gc",it,"\n",exp.gc[[it]],"\n")
			(exp.ht[[it]] <- t(exp.Gamma[[it]])%*%exp.mu[[it]]%*%exp.Gamma[[it]])
			#cat("ht",it,"\n",exp.ht[[it]],"\n")
			(exp.hc[[it]] <- .5*(exp.ht[[it]] + t(exp.ht[[it]])))  
			#cat("hc",it,"\n",exp.hc[[it]],"\n")
			(exp.w[[it]] <- 2*f.t0%*%exp.Omega[[it]] + f.t0%*%exp.VY[[it]]%*%exp.mu[[it]]%*%exp.Omega[[it]] + f.t0%*%exp.VY[[it]]%*%t(exp.mu[[it]])%*%exp.Omega[[it]])
			#cat("w", "\n")
			#print(exp.w[[it]])
			(exp.v[[it]] <- 2*f.t0%*%exp.Gamma[[it]] + f.t0%*%exp.VY[[it]]%*%exp.mu[[it]]%*%exp.Gamma[[it]] + f.t0%*%exp.VY[[it]]%*%t(exp.mu[[it]])%*%exp.Gamma[[it]])  
			#cat("v", "\n")
			#print(exp.v[[it]])
			(exp.VF[[it]] <- 2*f.t0%*%exp.VY[[it]]%*%t(f.t0) + f.t0%*%exp.VY[[it]]%*%exp.mu[[it]]%*%exp.VY[[it]]%*%t(f.t0) +  f.t0%*%exp.VY[[it]]%*%t(exp.mu[[it]])%*%exp.VY[[it]]%*%t(f.t0))
			#cat("VF","\n")
			#print(exp.VF[[it]])
			(exp.itlo[[it]] <- t(exp.Gamma[[it]])%*%exp.mu[[it]]%*%exp.Omega[[it]])
			
			(exp.itol[[it]] <- t(exp.Omega[[it]])%*%exp.mu[[it]]%*%exp.Gamma[[it]])

			(exp.ic[[it]] <- .25*(exp.itol[[it]] + t(exp.itol[[it]]) + exp.itlo[[it]] + t(exp.itlo[[it]])))
			#cat("ic",it,"\n",exp.ic[[it]],"\n")
			(exp.VGO[[it]] <- 2*delta.t0%*%k.t0%*%t(delta.t0) + 4*delta.t0%*%exp.gc[[it]]%*%t(delta.t0))
			(exp.VGL[[it]] <- 2*a.t0%*%j.t0%*%t(a.t0) + 4*a.t0%*%exp.hc[[it]]%*%t(a.t0))
			(exp.COVLO[[it]] <- 4*delta.t0%*%t(exp.ic[[it]])%*%t(a.t0) + 4*a.t0%*%exp.ic[[it]]%*%t(delta.t0))}


	Omega=exp.Omega[[gens]]
	#theta=exp.thetaNT[[gens]],
	Gamma=exp.Gamma[[gens]]
	gc=exp.gc[[gens]]
	hc=exp.hc[[gens]]
	itlo=exp.itlo[[gens]]
	itol=exp.itol[[gens]]
	ic=exp.ic[[gens]]
	VF=exp.VF[[gens]]
	w=exp.w[[gens]]
	v=exp.v[[gens]]
	VY=exp.VY[[gens]]
	mate.cov=mate.cov[[gens]]
	mu=exp.mu[[gens]]
	VGO=exp.VGO[[gens]]
	VGL=exp.VGL[[gens]]
	COVLO=exp.COVLO[[gens]]
	delta=delta.t0
	a=a.t0
	k=k2.matrix/2
	j=k2.matrix/2
	f=f.t0
	gt=exp.gt[[gens]]
	ht=exp.ht[[gens]]


	thetaNT <- 2 * delta %*% gc + 2 * a %*% ic + .5 * w
	thetaT  <- 1 * delta %*% k + thetaNT
	Yp_PGSm <- VY %*% mu %*% Omega
	Ym_PGSp <- VY %*% t(mu) %*% Omega
	Yp_Ym   <- VY %*% mu %*% VY
	Ym_Yp   <- VY %*% t(mu) %*% VY
	Yo_Yp   <- delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY
	Yo_Ym   <- delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY


	CMatrix <- rbind(
		#     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
		cbind(VY,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
		cbind(Ym_Yp,     VY,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
		cbind(Yo_Yp,     Yo_Ym,     VY,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
		cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
		cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
		cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
		cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc))  #NTm1 NTm2

	colnames(CMatrix) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
	rownames(CMatrix) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
	Mean <- rep(0, nrow(CMatrix))
	names(Mean) <- colnames(CMatrix)
	library(MASS)
	num_samples <- 6.4e4
	for (iter in 1:500){
	    set.seed(123+iter) 
	    samples <- mvrnorm(n = num_samples, mu = rep(0, nrow(CMatrix)), Sigma = CMatrix, empirical = FALSE)
	    # Save each sample to a file as a tsv
	    sample_file <- paste0("Data/Paper/uniModelBias/MVN/", conditionNames[condition],"/samples_iter_empi_", iter, ".tsv")
		# check if the file already exists
		if (file.exists(sample_file)) {
			# if it exists, skip saving
			next
		}
	    write.table(samples, file = sample_file, sep = "\t", row.names = FALSE, col.names = colnames(CMatrix), quote = FALSE)

	}
	summary_list <- list()
	for (i in 1:500){
	    DataPath <- paste0("Data/Paper/uniModelBias/MVN/", conditionNames[condition],"/samples_iter_empi_", i, ".tsv")
		testData <- read.table(DataPath, header = TRUE, sep = "\t")
		# subset the data to only include variables whose names end with 1
        testData <- testData[, grepl("1$", colnames(testData))]
        # arrange the order of columns to match the OpenMx script
        testData <- testData[, c("NTm1", "Tm1", "NTp1", "Tp1", "Ym1", "Yp1", "Yo1")]
        # fit the model
        fit <- fitUniSEMPGS(testData, max.cores = 2)
	    summary_list[[paste0("iter_", i)]] <- fit
	    cat("\nMVN samples iteration", i, "has been fitted\n")
	}

	# save the list of summaries
	save_path <- paste0("Analysis/Paper/uniModelBias/bias", conditionNames[condition] ,"_MVN_summary_list_500.rds")
	saveRDS(summary_list, save_path)
}