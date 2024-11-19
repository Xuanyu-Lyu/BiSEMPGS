
save_dir_data <- "/scratch/alpine/xuly4739/BiSEMPGS/Data"
save_path_expected <- "/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Data/Paper/Expected"
if (!dir.exists(save_path_expected)){
dir.create(save_path_expected)
}   

# x <- readRDS("Data/testdata/loop1.rds")

# names(x)
# x[[1]][[length(x[[1]])]]
# x[[4]] 
# x$COVARIANCES
# cor(as.data.frame(x[[4]])[,c("Y1P", "Y2P", "Y1M", "Y2M")])[1:2,3:4]

getObsExp <- function(simulatedRDS, delta, a, f, k){
    gens = 2
    delta.t0 = delta
    a.t0 = a
    f.t0 = f
    k.t0 = k
    j.t0 = k.t0
    target <- simulatedRDS[[1]][[length(simulatedRDS[[1]])]]
    #Define and initialize at t0 parameters that evolve
    exp.gc <- exp.hc <- exp.ic <- exp.gt <- exp.ht <- exp.itlo <- exp.itol <- exp.w <- exp.v <- exp.VY <- exp.VY.diagonal <- exp.corVY<- exp.mu <- exp.VF <- exp.Omega <- exp.Gamma <- exp.thetaT <- exp.thetaNT <- exp.thetaLT <- exp.thetaLNT <- mate.cov <- mate.cor <- exp.VGO <- exp.VGL <- exp.COVLO <- vector(mode="list",length=gens-1)
    #initialize parameters that begin as 0's in base population (before VT and AM)
    exp.VGO[[1]] <- exp.VGL[[1]] <- exp.COVLO[[1]] <- exp.gc[[1]] <- exp.hc[[1]] <- exp.ic[[1]] <- exp.w[[1]] <- exp.v[[1]] <- matrix(0,nrow=2,ncol=2)
    # define the true 4X4 covariance matrix for Vy at t0
    exp.covY.full <- exp.corY.full <- exp.covY.full.diagonal <- vector(mode="list",length=gens-1)
    #initialize parameters that are non-zero in base population (think of these as the first parents to pass on VT and engage in AM)
    (exp.VY[[1]] <- target$VP)
    (exp.VY.diagonal[[1]] <- diag(diag(exp.VY[[1]])))
    (exp.VF[[1]] <- target$covF)
    mate.cor[[1]] <- cor(as.data.frame(simulatedRDS[[4]]))[c("Y1P", "Y2P"), c("Y1M", "Y2M")]
    rmate.t0 <- mate.cor[[1]]
    print(rmate.t0)
    covE.t0 <- target$covE
    exp.w[[1]] <- target$w
    exp.v[[1]] <- target$q
    exp.Omega[[1]] <- target$omega
    exp.Gamma[[1]] <- target$gamma
     mate.cov[[1]] <- cov(as.data.frame(simulatedRDS[[4]]))[c("Y1P", "Y2P"), c("Y1M", "Y2M")]
    (exp.mu[[1]] <- solve(exp.VY[[1]]) %*% mate.cov[[1]] %*% solve(t(exp.VY[[1]])))

     (exp.covY.full[[1]] <- rbind(cbind(exp.VY[[1]], t(mate.cov[[1]])),cbind(mate.cov[[1]], exp.VY[[1]])))
     (exp.corY.full[[1]] <- cov2cor(exp.covY.full[[1]]))
    # mate.cor[[1]] <- mate.cov[[1]]/COVY
    # mate.cor[[1]][1,2] <- mate.cov[[1]][1,2]/sqrt(COVY[1,1]*COVY[2,2])
    # mate.cor[[1]][2,1] <- mate.cov[[1]][2,1]/sqrt(COVY[1,1]*COVY[2,2])
    # mate.cor[[1]]


    for (it in 2:gens){
            (exp.gt[[it]] <- t(exp.Omega[[it-1]])%*%exp.mu[[it-1]]%*%exp.Omega[[it-1]])
            cat("gt",it,"\n",exp.gt[[it]],"\n")
            (exp.gc[[it]] <- .5*(exp.gt[[it]] + t(exp.gt[[it]])))
            cat("gc",it,"\n",exp.gc[[it]],"\n")
            (exp.ht[[it]] <- t(exp.Gamma[[it-1]])%*%exp.mu[[it-1]]%*%exp.Gamma[[it-1]])
            cat("ht",it,"\n",exp.ht[[it]],"\n")
            (exp.hc[[it]] <- .5*(exp.ht[[it]] + t(exp.ht[[it]])))  
            cat("hc",it,"\n",exp.hc[[it]],"\n")


            (exp.itlo[[it]] <- t(exp.Gamma[[it-1]])%*%exp.mu[[it-1]]%*%exp.Omega[[it-1]])
            
            (exp.itol[[it]] <- t(exp.Omega[[it-1]])%*%exp.mu[[it-1]]%*%exp.Gamma[[it-1]])

            (exp.ic[[it]] <- .5*(exp.itol[[it]] + t(exp.itlo[[it]])))
            cat("ic",it,"\n",exp.ic[[it]],"\n")
            
            (exp.w[[it]] <- 2*f.t0%*%exp.Omega[[it-1]] + f.t0%*%exp.VY[[it-1]]%*%exp.mu[[it-1]]%*%exp.Omega[[it-1]] + f.t0%*%exp.VY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.Omega[[it-1]])
            cat("w", "\n")
            print(exp.w[[it]])
            (exp.v[[it]] <- 2*f.t0%*%exp.Gamma[[it-1]] + f.t0%*%exp.VY[[it-1]]%*%exp.mu[[it-1]]%*%exp.Gamma[[it-1]] + f.t0%*%exp.VY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.Gamma[[it-1]])  
            cat("v", "\n")
            print(exp.v[[it]])

            (exp.Omega[[it]] <- 2*delta.t0%*%exp.gc[[it]] + delta.t0 %*% k.t0 + .5*exp.w[[it]] + 2*a.t0%*%exp.ic[[it]])
            cat("Omega",it,"\n",exp.Omega[[it]],"\n")
            (exp.Gamma[[it]] <- 2*a.t0%*%exp.hc[[it-1]] + 2*delta.t0%*%t(exp.ic[[it]]) + a.t0%*%j.t0 + .5*exp.v[[it]])
            cat("Gamma",it,"\n",exp.Gamma[[it]],"\n")
            (exp.VF[[it]] <- 2*f.t0%*%exp.VY[[it-1]]%*%t(f.t0) + f.t0%*%exp.VY[[it-1]]%*%exp.mu[[it-1]]%*%exp.VY[[it-1]]%*%t(f.t0) +  f.t0%*%exp.VY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.VY[[it-1]]%*%t(f.t0))
            cat("VF","\n")
            print(exp.VF[[it]])

            (exp.VY[[it]] <- 2 * delta.t0 %*% t(exp.Omega[[it]]) + 2 * a.t0 %*% t(exp.Gamma[[it]]) + exp.w[[it]] %*% t(delta.t0) + exp.v[[it]] %*% t(a.t0) + exp.VF[[it]] + covE.t0 )
            cat("VY","\n")
            print(exp.VY[[it]])

            # exp.VY.diagonal[[it]] <- exp.VY[[it]]
            # exp.VY.diagonal[[it]][1,2] <- 0
            # exp.VY.diagonal[[it]][2,1] <- 0

            # exp.covY.full.diagonal[[it]] <- cbind(rbind(exp.VY.diagonal[[it]],matrix(0,2,2)),rbind(matrix(0,2,2),exp.VY.diagonal[[it]]))

            # exp.corVY[[it]] <- cov2cor(exp.VY[[it]])
            # exp.corY.full[[it]] <- exp.corY.full[[1]]
            # exp.corY.full[[it]][1,2] <- exp.corVY[[it]][1,2]
            # exp.corY.full[[it]][2,1] <- exp.corVY[[it]][2,1]
            # exp.corY.full[[it]][3,4] <- exp.corVY[[it]][1,2]
            # exp.corY.full[[it]][4,3] <- exp.corVY[[it]][2,1]
            # #exp.covY.full[[it]] <- cor2cov(exp.corY.full[[it]],exp.VY[[it]][1,1],exp.VY[[it]][2,2],exp.VY[[it]][1,1],exp.VY[[it]][2,2]) 
            # exp.covY.full[[it]] <- sqrt(exp.covY.full.diagonal[[it]]) %*% exp.corY.full[[it]] %*% sqrt(exp.covY.full.diagonal[[it]])
            # mate.cov[[it]] <- exp.covY.full[[it]][3:4,1:2] ## this is the mate.cov fixing mate.cor
            # mate.cor[[it]] <- exp.corY.full[[1]][3:4,1:2]
            (exp.mu[[it]] <- solve(exp.VY[[it]]) %*% mate.cov[[it-1]] %*% solve(t(exp.VY[[it]]))) 
            cat("mu",it,"\n",exp.mu[[it]],"\n")
            

            
            (exp.VGO[[it]] <- 2*delta.t0%*%k.t0%*%t(delta.t0) + 4*delta.t0%*%exp.gc[[it]]%*%t(delta.t0))
            (exp.VGL[[it]] <- 2*a.t0%*%j.t0%*%t(a.t0) + 4*a.t0%*%exp.hc[[it]]%*%t(a.t0))
            (exp.COVLO[[it]] <- 4*delta.t0%*%t(exp.ic[[it]])%*%t(a.t0) + 4*a.t0%*%exp.ic[[it]]%*%t(delta.t0))

    }
        EXP <- list(
            Omega=exp.Omega[[gens-1]],
            #theta=exp.thetaNT[[gens]],
            Gamma=exp.Gamma[[gens-1]],
            gc=exp.gc[[gens]],
            hc=exp.hc[[gens]],
            itlo=exp.itlo[[gens]],
            itol=exp.itol[[gens]],
            ic=exp.ic[[gens]],
            VF=exp.VF[[gens-1]],
            w=exp.w[[gens-1]],
            v=exp.v[[gens-1]],
            VY=exp.VY[[gens-1]],
            mate.cov=mate.cov[[gens-1]],
            mu=exp.mu[[gens-1]],
            VGO=exp.VGO[[gens]],
            VGL=exp.VGL[[gens]],
            COVLO=exp.COVLO[[gens]]
            )
        return(EXP)

            # (cov.fin.gen <- exp.VY[[gens]]%*%exp.mu[[gens]]%*%exp.VY[[gens]])
            # cor.fin.gen <- cov.fin.gen/exp.VY[[gens]]
            # cor.fin.gen[1,2] <- cov.fin.gen[1,2]/(sqrt(exp.VY[[gens]][1,1])*sqrt(exp.VY[[gens]][2,2]))
            # cor.fin.gen[2,1] <- cov.fin.gen[2,1]/(sqrt(exp.VY[[gens]][1,1])*sqrt(exp.VY[[gens]][2,2]))
            # cor.fin.gen
            # rmate.t0 
}
# test the function
#getObsExp(x)

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
	#(vg.obs <- delta.mat %*% k2.matrix %*% t(delta.mat)) #Check
	#matrix(c(vg.obs1,covg.obs,covg.obs,vg.obs2),nrow=2) #Check

	#Create lower a matrix
	(vg.lat1 <- vg1*(prop.h2.latent1))
	(vg.lat2 <- vg2*(prop.h2.latent2))
	(covg.lat <- rg*sqrt(vg.lat1)*sqrt(vg.lat2))
	a11 <- sqrt(vg.lat1)
	#a21 <- covg.lat/a11
	a21 <- 0 #new way
	a22 <- sqrt(vg.lat2 - a21^2)
	(a.mat <- matrix(c(a11,0,a21,a22),byrow=T,nrow=2))
	#(vg.lat <- a.mat %*% k2.matrix %*% t(a.mat)) #Check
	matrix(c(vg.lat1,covg.lat,covg.lat,vg.lat2),nrow=2) #Check

	#Find the total genetic covariance
	#(covg <- vg.obs[1,2] + vg.lat[1,2])
	#(covg.mat <- vg.obs + vg.lat)
	#a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) #check

	#Create E lower matrix
	ve1 <- 1 - vg1
	ve2 <- 1 - vg2
	(cove <- re*sqrt(ve1*ve2))
	(cove.mat <- matrix(c(ve1,cove,cove,ve2),nrow=2,byrow=T))
    (f.mat <- matrix(c(f11,f12,f21,f22),nrow=2,byrow=T))

    data_path <- paste0(save_dir_data, "/", conditionNames[condition])
    l_files <- list.files(data_path, pattern = ".rds")
    l_finalGen <- list()
    for (loop in 1: length(l_files)){
        cat("loop", loop, "\n")
        # get the number of .rds files in the folder
        # read all the .rds file one by one and save them as .txt files
        # load the data
        x <- readRDS(paste0(data_path, "/", l_files[loop]))
        exp <- getObsExp(x, delta = delta.mat, a = a.mat, f = f.mat, k = k2.matrix/2)
        l_finalGen[[loop]] <- exp
    }
    saveRDS(l_finalGen, file = paste0(save_path_expected, "/", conditionNames[condition], "_finalGen.rds"))

}
    #