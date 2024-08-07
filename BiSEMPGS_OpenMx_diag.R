# this is the OpenMx script to fit the Bivariate SEM-PGS model
# Author: Xuanyu Lyu
# Date: 05/31/2024

# load the packages
library(OpenMx)
library(data.table)
library(stringr)

# Specify Options:
    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")
    #mxOption(NULL,"mxByRow","TRUE")

# Load the simulated data for this demonstration:
    Example_Data  <- fread("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/normal_full/loop1.txt", header = T)

    cov(Example_Data, use="pairwise.complete.obs")

# Create variables and define the algebra for each variables

    VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,F,F,T), values=c(3,0,0,1), label=c("VY11", "VY12", "VY12","VY22"), name="VY") # Phenotypic variance
	VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,F,F,T), values=c(1.3,0,0,.7), label=c("VF11", "VF12", "VF12","VF22"), name="VF") # Variance due to VT
    VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.7,.0,0,.5), label=c("VE11", "VE12", "VE12","VE22"), name="VE") # Residual variance

    VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")

    VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
    VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
# Genetic effects:
    delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta") # Effect of PGS on phen
    a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.2,.1), label=c("a11", "a22"),    name="a")     # Effect of latent PGS on phen
    k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0,0,.5), label=c("k11", "k12", "k12","k22"),    name="k")     # PGS variance (if no AM)
    j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0,0,.5), label=c("j11", "j12", "j12","j22"),    name="j")     # Latent PGS variance (if no AM)
    Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.6,0,0,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega") # Within-person PGS-Phen covariance
    Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.5,0,0,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma") # Within-person latent PGS-Phen covariance

    Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
    Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

	Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
# Assortative mating effects:
    mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.15,0,0,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu") # AM co-path coefficient
    gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.2,0,0,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt")  # Increase in cross-mate PGS (co)variances from AM
    ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.15,0,0,.5), label=c("ht11", "ht21", "ht12","ht22"),  name="ht")  # Increase in cross-mate latent PGS (co)variances from AM
    gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.4,0,0,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc")  # Increase in within-mate PGS (co)variances from AM
    hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.7,0,0,.5), label=c("hc11", "hc12", "hc12","hc22"),  name="hc")  # Increase in within-mate latent PGS (co)variances from AM
    gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
    ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
    gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
    hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
    gchc_constraint_Algebra <- mxAlgebra(solve(delta) %*% a %*% hc %*% t(a) %*% solve(t(delta)), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

    itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.25,0,0,.3), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo") 
    itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.55,0,0,.13), label=c("itol11", "itol21", "itol12","itol22"), name="itol")
    ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.47,0,0,.35), label=c("ic11", "ic12", "ic12","ic22"), name="ic") 
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
    f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.4,0,0,.3), label=c("f11", "f21","f12","f22"),  name="f") # Vertical Transmission effect
    w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(.3,0,0,.51), label=c("w11", "w12", "w21","w22"),  name="w") # Genetic nurture
    v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,F,F,T), values=c(2,0,0,.54), label=c("v11", "v12", "v21","v22"),  name="v") # Latent nurture
    w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
    v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
    
    v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')

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

# Create fit function:
    FitFunctionML <- mxFitFunctionML()
    # Create fit function:
    #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
# Specify what parameters we're going to be including in our model:
    Params <- list(
                VY, VF, VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v, 
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, 
                VY_Constraint, VF_Constraint, Omega_Constraint, Gamma_Constraint, gt_constraint, ht_constraint, gc_constraint, hc_constraint, gchc_constraint, itlo_constraint, itol_constraint, ic_constraint, v_constraint, 
                thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                CovMatrix, Means, ModelExpectations, FitFunctionML)
# Create the model:
    options(warning.length = 8000)
    Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)

    fitModel1 <- mxRun(Model1, intervals=F, silent=F)

    summary(fitModel1)


# covariance =  matrix(c(    # 14x14
#    3,   0, 1.35,   0, 2.233,   0, 0.6,   0, 0.6,   0, 0.27,   0, 0.27,   0
# ,   0,   5,   0, 7.5,   0, 9.175,   0, 0.5,   0, 0.5,   0, 0.75,   0, 0.75
# , 1.35,   0,   3,   0, 2.233,   0, 0.27,   0, 0.27,   0, 0.6,   0, 0.6,   0
# ,   0, 7.5,   0,   5,   0, 9.175,   0, 0.75,   0, 0.75,   0, 0.5,   0, 0.5
# , 2.233,   0, 2.233,   0,   3,   0, 0.858, 0.05, 0.658, 0.05, 0.858, 0.05, 0.658, 0.05
# ,   0, 9.175,   0, 9.175,   0,   5, 0.05, 0.955, 0.05, 0.655, 0.05, 0.955, 0.05, 0.655
# , 0.6,   0, 0.27,   0, 0.858, 0.05, 0.9,   0, 0.4,   0, 0.2,   0, 0.2,   0
# ,   0, 0.5,   0, 0.75, 0.05, 0.955,   0, 0.6,   0, 0.1,   0, 0.1,   0, 0.1
# , 0.6,   0, 0.27,   0, 0.658, 0.05, 0.4,   0, 0.9,   0, 0.2,   0, 0.2,   0
# ,   0, 0.5,   0, 0.75, 0.05, 0.655,   0, 0.1,   0, 0.6,   0, 0.1,   0, 0.1
# , 0.27,   0, 0.6,   0, 0.858, 0.05, 0.2,   0, 0.2,   0, 0.9,   0, 0.4,   0
# ,   0, 0.75,   0, 0.5, 0.05, 0.955,   0, 0.1,   0, 0.1,   0, 0.6,   0, 0.1
# , 0.27,   0, 0.6,   0, 0.658, 0.05, 0.2,   0, 0.2,   0, 0.4,   0, 0.9,   0
# ,   0, 0.75,   0, 0.5, 0.05, 0.655,   0, 0.1,   0, 0.1,   0, 0.1,   0, 0.6), byrow=TRUE, nrow=14, ncol=14)
# #[c(5:14),c(5:14)]
# eigen(covariance)$values
 eigen(mxGetExpected(Model1, "covariance")[c(2,4,6,8,10,12,14),c(2,4,6,8,10,12,14)])$values

 eigen(mxGetExpected(Model1, "covariance")[c(1,3,5,7,9,11,13),c(1,3,5,7,9,11,13)])$values

# eigen(mxGetExpected(Model1, "covariance")[c(2,4,6,8,10,12,14),c(2,4,6,8,10,12,14)])$values

eigen(mxGetExpected(Model1, "covariance"))$values
