conditionNames <- c("onlyAM", "onlyVT", "bothAMVT")

startingParamList1 <- list(vg1 = rep(.64,3),
						   vg2 = rep(.36,3),
						   am11 = rep(0.4,3),
						   am12 = c(.2,0,.2),
						   am21 = c(.2,0,.2), 
						   am22 = rep(0.3,3),
						   f11 = rep(0.15,3),
						   f12 = c(0,.1,.1),
						   f21 = c(0,.05,.05), 
						   f22 = rep(0.1,3),
						   Nfam = rep(8e4, 3),
						   rg = rep(0,3),
						   re = rep(0,3),
						   prop.h2.latent1 = rep(.60/.64,3),
						   prop.h2.latent2 = rep(.8,3))
condition = 3

   pathName = "testdata"
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

# save the covariance matrix to a file
save_path <- "PaperScripts/expCov.csv"
write.csv( CMatrix, file = save_path, row.names = TRUE)

# remove all environment variables

#rm(list = ls())
save_path <- "PaperScripts/expCov.csv"
CMatrix <- read.csv(save_path, row.names = 1)
Mean <- rep(0, nrow(CMatrix))
names(Mean) <- colnames(CMatrix)
library(MASS)
# generate multivariate normal samples
 # for reproducibility

#source("PaperScripts/04-OpenMxFunctions.R")
fitBiSEMPGS_m2_tol <- function(data_path,feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F, Mean = NULL){
        # load the packages
    library(OpenMx)
    library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))

    #mxOption(NULL,"Number of Threads","4")
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    #mxOption(NULL,"Analytic Gradients","No")

        Example_Data  <- fread(data_path, header = T)

        #cov(Example_Data, use="pairwise.complete.obs")

     # Create variables and define the algebra for each variables
        #VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=FALSE, values=NA, label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        #VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.20,0.06,0.06,.04), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.1) # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,.06,0.06,.4), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")
        #VE_Algebra <- mxAlgebra(VY - (2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra), name="VE_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        #VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
        #VE_Constraint    <- mxConstraint(VE == VE_Algebra,       name='VE_Constraint')

    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.8,.55), label=c("a11", "a22"),    name="a", lbound = c(.2,.2))     # Effect of latent PGS on phen
        k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.02,0.02,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.05)     # PGS variance (if no AM)
        j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.03,0.03,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.05)     # Latent PGS variance (if no AM)
        Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.05) # Within-person PGS-Phen covariance
        Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.10,0.15,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.05) # Within-person latent PGS-Phen covariance

        Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
        Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

        Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
        Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
        
        adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
        adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    # add a constraint to make k=j
        j_Algebra <- mxAlgebra(k, name = "j_Algebra")
        j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects:
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0,0.05,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.2) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.05,0.01,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.1)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.1)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.1)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.1)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.05,0.04,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.1) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.03,0.036,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.1)
        ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.05,0.05,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.1) 
        itlo_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Omega, name="itlo_Algebra") # E.g., cov(TPO, TML)
        itol_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Gamma, name="itol_Algebra") # E.g., cov(TPL, TMO)
        ic_Algebra <- mxAlgebra(.25 * (itlo + t(itlo) + itol + t(itol)), name="ic_Algebra") # ic should be symmetric
        
        

        gt_constraint <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
        ht_constraint <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
        gc_constraint <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
        hc_constraint <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
        gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
        itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
        itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
        ic_constraint <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
        f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.04,0.11,.09), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.05) # Vertical Transmission effect
        w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.05) # Genetic nurture
        v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.07,.08), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.05) # Latent nurture
        w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
        v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
        wv_constraint_algebra <- mxAlgebra((w * sqrt(2*delta%*%k%*%t(delta))/sqrt(2*a%*%j%*%t(a))), name='wv_constraint_algebra')

        v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')
        w_constraint <- mxConstraint(w == w_Algebra, name='w_constraint')
        wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')
    # Between-people covariances
        thetaNT <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + .5 * w, name="thetaNT")
        thetaT  <- mxAlgebra(1 * delta %*% k + thetaNT, name="thetaT")
        Yp_PGSm <- mxAlgebra(VY %*% mu %*% Omega, name="Yp_PGSm")
        Ym_PGSp <- mxAlgebra(VY %*% t(mu) %*% Omega, name="Ym_PGSp")
        Yp_Ym   <- mxAlgebra(VY %*% mu %*% VY,    name="Yp_Ym")
        Ym_Yp   <- mxAlgebra(VY %*% t(mu) %*% VY, name="Ym_Yp")
        Yo_Yp   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY, name = "Yo_Yp")
        Yo_Ym   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY, name = "Yo_Ym")
    # Expected covariances matrix
        CovMatrix <- mxAlgebra(rbind(
            #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
            cbind(VY,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
            cbind(Ym_Yp,     VY,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
            cbind(Yo_Yp,     Yo_Ym,     VY,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc)),   #NTm1 NTm2
            dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )
        #Example_Data_Mx <- mxData(observed=cov(Example_Data), type="cov", numObs=8e4, means=Mean) 

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, 
                    #VF, 
                    VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                    VY_Algebra, 
                    VF_Algebra,  
                    #VE_Algebra,
                    Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    #VF_Constraint, 
                    #VE_Constraint,
                    #Omega_Constraint, 
                    Gamma_Constraint, 
                    #adelta_Constraint,
                    j_constraint,
                    #gt_constraint, 
                    ht_constraint, 
                    #gc_constraint, 
                    hc_constraint, 
                    #gchc_constraint, 
                    itlo_constraint, 
                    itol_constraint, 
                    ic_constraint, 
                    v_constraint, 
                    w_constraint,
                    #wv_constraint,
                    thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                    CovMatrix, Means, ModelExpectations, FitFunctionML)
    # Create the model:
        options(warning.length = 8000)
        Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)

        fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, showInits = F, exhaustive = exhaustive, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
        #fitModel1 <- mxRun(Model1)
        #return(fitModel1)
        return(summary(fitModel1, verbose = TRUE))

}

num_samples <- 6.4e4  # number of samples to generate


# for (iter in 1:500){
#     set.seed(123+iter) 
#     samples <- mvrnorm(n = num_samples, mu = rep(0, nrow(CMatrix)), Sigma = CMatrix, empirical = FALSE)
#     # Save each sample to a file as a tsv
#     sample_file <- paste0("Data/Paper/MVN/samples_iter_empi_", iter, ".tsv")
#     write.table(samples, file = sample_file, sep = "\t", row.names = FALSE, col.names = colnames(CMatrix), quote = FALSE)

# }
# summary_list <- list()
# for (i in 1:500){
#     sample_file <- paste0("Data/Paper/MVN/samples_iter_empi_", i, ".tsv")
#     fit <- fitBiSEMPGS_m2_tol(sample_file, 
#                               feaTol = 1e-6, 
#                               optTol = 1e-9,
#                               jitterMean = 0.5,
#                               jitterVar = .1,
#                               exhaustive = FALSE,
#                               extraTries = 5,
#                               Mean = Mean)
#     summary_list[[paste0("iter_", i)]] <- fit
#     cat("\nMVN samples iteration", i, "has been fitted\n")
# }

# # save the list of summaries
# save_path <- "Analysis/Paper/MVN/paper_MVN_summary_list_500.rds"
# saveRDS(summary_list, save_path)

# load the summary list
summary_list <- readRDS("Analysis/Paper/MVN/paper_MVN_summary_list_500.rds")
# extract all the status code of openmx and put them into a vector
status_codes <- sapply(summary_list, function(x) x$statusCode)
summary(status_codes)
# extract all the estimates in the list and put each parameter as a column in a data frame
# Initialize an empty 78 column data frame
df <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
colnames(df) <- summary_list[[2]]$parameters$name
#colnames(df) 
# Loop over the elements in the summary_list
for(i in 1:length(summary_list)) {
    for(j in 1:nrow(summary_list[[i]]$parameters)){
        if (!is.null(summary_list[[i]]$parameters$Estimate[j])) {
            df[i,j] <- summary_list[[i]]$parameters$Estimate[j]
        } else {
            print(paste("NULL value at i =", i, "and j =", j))
        }
    }
}
df$status_codes <- status_codes
psych::describe(df, trim = 0) |> print(digits = 4)


# Code to create a true value named vector
break_into_named_vector <- function(mean_params, param) {
  named_vector <- c()
  matrix <- mean_params
  for (i in 1:nrow(matrix)) {
    for (p in 1:ncol(matrix)) {
      new_name <- paste0(param, i, p)
      named_vector[new_name] <- matrix[i, p]
    }
  }
  return(named_vector)
}



v_VY <- break_into_named_vector(VY, "VY")
v_f <- break_into_named_vector(f, "f")
v_VF <- break_into_named_vector(VF, "VF")
v_delta <- break_into_named_vector(delta, "delta")
v_a <- break_into_named_vector(a, "a")
v_k <- break_into_named_vector(k, "k")
v_j <- break_into_named_vector(k, "j")
v_Omega <- break_into_named_vector(Omega, "Omega")
v_Gamma <- break_into_named_vector(Gamma, "Gamma")
v_mu <- break_into_named_vector(mu, "mu")
v_gt <- break_into_named_vector(gt, "gt")
v_ht <- break_into_named_vector(ht, "ht")
v_gc <- break_into_named_vector(gc, "gc")
v_hc <- break_into_named_vector(hc, "hc")
v_itlo <- break_into_named_vector(itlo, "itlo")
v_itol <- break_into_named_vector(itol, "itol")
v_ic <- break_into_named_vector(ic, "ic")
v_w <- break_into_named_vector(w, "w")
v_v <- break_into_named_vector(v, "v")
v_VGO <- break_into_named_vector(VGO, "VGO")
v_VGL <- break_into_named_vector(VGL, "VGL")
v_COVLO <- break_into_named_vector(COVLO, "COVLO")

v_true_values <- c(v_VY, v_f, v_VF, v_delta, v_a, v_k, v_j, v_Omega, v_Gamma, v_mu, v_gt, v_ht, v_gc, v_hc, v_itlo, v_itol, v_ic, v_w, v_v, v_VGO, v_VGL, v_COVLO)

v_true_values


# chat's version of the permutation test
bootstrap_test_median <- function(data, m0, n_boot = 10000, seed = NULL) {
  # data:   Numeric vector of observations
  # m0:     Hypothesized median under H0
  # n_boot: Number of bootstrap iterations
  # seed:   Optional seed for reproducibility
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  data <- as.numeric(data)
  n <- length(data)
  
  # 1. Observed median and observed test statistic
  observed_median <- median(data, na.rm = TRUE)
  observed_stat <- observed_median - m0  # difference from hypothesized median
  observed_abs_stat <- abs(observed_stat)
  
  # 2. Center the data to reflect H0: shift so that its median = m0
  #    (remove the observed median, then add the hypothesized median)
  centered_data <- data - observed_median + m0
  
  # 3. Bootstrap replicates
  count_extreme <- 0
  for (i in seq_len(n_boot)) {
    # Sample with replacement from the centered data
    b_sample <- sample(centered_data, size = n, replace = TRUE)
    b_median <- median(b_sample, na.rm = TRUE)
    
    # Test statistic for the bootstrap sample
    b_stat <- b_median - m0
    if (abs(b_stat) >= observed_abs_stat) {
      count_extreme <- count_extreme + 1
    }
  }
  
  # 4. Two-sided p-value
  p_value <- count_extreme / n_boot
  
  return(p_value)
}

# do a permutation test to see if the means of the estimates are significantly different from the true values
permutationTest <- function(X, null_median, n_iter=1e5) {
  # Step 1: Compute the observed median
  m_obs <- median(X)
  # Step 2: Center the data so the median becomes the null median
  X_centered <- X - m_obs + null_median
  # Step 3: Generate the null distribution of the median using bootstrap resampling
  boot_medians <- replicate(n_iter, {
    sample_data <- sample(X_centered, size = length(X_centered), replace = TRUE)
    median(sample_data)
    #print(median(sample_data))
  })
  # Step 4: Calculate the two-tailed p-value
  p_lower <- mean(boot_medians <= m_obs)
  p_upper <- mean(boot_medians >= m_obs)
  p_value <- 2 * min(p_lower, p_upper)
  #p_value <- min(p_value, 1)  # Ensure p-value does not exceed 1
  # Return the results as a list
#   return(list(
#     observed_median = m_obs,
#     p_value = p_value
#     #,boot_medians = boot_medians
#   ))
    # Return the p value
    return(p_value)
}

# a function to get descriptive statistics for one parameter and return a vector
getDescriptive <- function(df, param, file_tv){
    med <- median(df[[param]], na.rm = TRUE)
    MAD <- mad(df[[param]], na.rm = TRUE)
    trueValue <- file_tv[names(file_tv) == param]
    #cat(param, "\t", med, "\t", MAD, "\t", trueValue, "\n")
    # significance test if median is significantly different from the true value
    p_value <- bootstrap_test_median( df[[param]], m0 = trueValue, n_boot = 10000)
    #p_value <- wilcox.test(df[[param]], mu = trueValue, alternative = "two.sided")$p.value

    # compute the total variance and the variance from randomness
    var_total <- var(df[[param]], na.rm = TRUE)
    var_random <- sum((df[[param]] - trueValue)^2)/(length(df[[param]]) - 1)
    var_systematic <- var_random - var_total
    p_systematic <- var_systematic/var_random

    final_v <- c(med, MAD, trueValue, p_value, p_systematic)
    
    names(final_v) <- c("median", "MAD", "trueValue", "p_value", "proportion_systematic")
    return(final_v)

}

# filter the data frame to only include the parameters that are in the df
v_true_values <- v_true_values[which(names(v_true_values) %in% colnames(df))]
v_true_values |> length()
 
 # for each parameter, get the descriptive statistics and put them into a data frame
descriptive_df <- data.frame(matrix(ncol = 5, nrow = length(v_true_values)))
colnames(descriptive_df) <- c("median", "MAD", "trueValue", "p_value", "proportion_systematic")

# fill the data frame with the descriptive statistics
for (i in 1:length(v_true_values)) {
    param <- names(v_true_values)[i]
    # get the true value from the vector
    trueValue <- v_true_values[i]
    # get the descriptive statistics for the parameter
    descriptive_df[i,] <- getDescriptive(df, param, v_true_values)
}
rownames(descriptive_df) <- names(v_true_values)

# make plots for each parameter
library(ggplot2)
library(patchwork)
# create a histogram for each parameter and combine them using the wrap plot function, also add the true value line, the median line, and the p value

# create a function to plot the histogram
plot_histogram <- function(df, param, trueValue, median, MAD, p_value) {
      # Calculate the maximum density value
    #     max_density <- max(density(df[[param]], na.rm = TRUE)$y)
  
    #   # Set the y-axis upper limit with a 5% expansion
    #     upper_limit <- max_density * 1.05
    p <- ggplot(df, aes_string(x = param)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "#1f77b4", alpha = 0.9) +
        geom_vline(xintercept = trueValue, color = "#d62728", linetype = "dashed", size = 1) +
        geom_vline(xintercept = median, color = "#b5cf6b", linetype = "dashed", size = 1) +
        #geom_vline(xintercept = median + MAD, color = "orange", linetype = "dashed", size = 1) +
       # geom_vline(xintercept = median - MAD, color = "orange", linetype = "dashed", size = 1) +
        scale_x_continuous(labels = function(x) sprintf("%.2f", x)) + 
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(title = paste(param),
             x = param,
             y = NULL) + 
             theme_minimal() +
        theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) +
            annotate("text", x = Inf, y = Inf, label = paste("p = ", round(p_value, 3)), color = "black", 
                    hjust = 1.1, vjust = 1.5)
    return(p)
}

# create a function to combine plots for specific parameters
combine_plots <- function(df, params, ncol=4) {
    plots <- list()
    for (param in params) {
        trueValue <- descriptive_df[param,"trueValue"]
        median <- descriptive_df[param,"median"]
        MAD <- descriptive_df[param,"MAD"]
        p_value <- descriptive_df[param,"p_value"]
        plots[[param]] <- plot_histogram(df, param, trueValue, median, MAD, p_value)
    }
    combined_plot <- wrap_plots(plots, ncol = ncol)
}

params1 <-  c("f11", "mu11", "a11", "delta11", "w11", "v11", "VY11","gc11")

plots1 <- combine_plots(df, params1, ncol = 4)
plots1

ggsave(paste0("Analysis/Paper/", "MVNfigure11.png") , plots1, width = 10, height = 6, type = "cairo-png", dpi = 600)

params2 <- c("f12", "mu12",  "w12", "v12", "VY12","gc12")

plots2 <- combine_plots(df, params2, ncol = 3)
plots2

ggsave(paste0("Analysis/Paper/", "MVNfigure12.png") , plots2, width = 10, height = 6, type = "cairo-png", dpi = 600)
