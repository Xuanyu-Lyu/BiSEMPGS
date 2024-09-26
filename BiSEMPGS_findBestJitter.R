# This is a script to find the best possible distribution for the mxTryHard function in OpenMx for BoSEMPGS model. 

CMatrix <- read.csv("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/ExpectedCMatrix.csv", header = TRUE, row.names = 1) |> as.matrix()
colnames(CMatrix) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
rownames(CMatrix) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
#CMatrix <- read.csv("ExpectedCMatrix.csv", header = TRUE, row.names = 1) |> as.matrix()

Mean <- rep(0, 14)
names(Mean) <- colnames(CMatrix)


# the openmx model we want to use
fitBiSEMPGS_m2_cov <- function(cov, Mean, jitterMean = .5, jitterVar = .1){
        # load the packages
    library(OpenMx)
    #library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance","1e-7")
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())
    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    #mxOption(NULL,"Analytic Gradients","No")

        cMatrix <- cov

        #cov(Example_Data, use="pairwise.complete.obs")

    # Create variables and define the algebra for each variables

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.20,0.06,0.06,.04), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.1) # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,.06,0.06,.84), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.34), label=c("a11", "a22"),    name="a", lbound = c(.1,.1))     # Effect of latent PGS on phen
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
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0.15,0.1,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.1) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.05,0.01,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.05)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.05)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.05)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.05)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.05,0.04,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.05) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.03,0.036,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.05)
        ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.05,0.05,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.05) 
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
            dimnames=list(colnames(cMatrix),colnames(cMatrix)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=cMatrix, type="cov", numObs=4.8e4, means=Mean) 

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, VF, VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                    VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    VF_Constraint, 
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

        fitModel1 <- mxTryHardWideSearch(Model1, extraTries = 30, OKstatuscodes = c(0,1), intervals=T, silent=F, showInits = F, exhaustive = F, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
        return(summary(fitModel1))

}

v_jitterMean <- c(seq(0,3, by = .1))
v_jitterVar <- c(seq(0,1, by = .05))

# v_jitterMean <- c(seq(0,1, by = .5))
# v_jitterVar <- c(seq(0,1, by = .5))

# a loop to run through all the possible combinations of jitterMean and jitterVar and save the results in a list
summary_list <- list()
for (i in 1:length(v_jitterMean)){
    for (j in 1:length(v_jitterVar)){
        summary_list[[paste0("jitterMean_", v_jitterMean[i])]][[paste0("_jitterVar_", v_jitterVar[j])]] <- fitBiSEMPGS_m2_cov(CMatrix, Mean, jitterMean = v_jitterMean[i], jitterVar = v_jitterVar[j])
    }
}


conditionNames <- c("Full_Model", "MeasurePgs30", "MeasurePgs10", "MeasurePgsFully", 
"f11-decrease", "f12-decrease", "f11.12.21.22-decrease", "am11-decrease", "am12-decrease", "am11.12.21.22-decrease", "Full_Model_.5latent")
save_path <- paste0("/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Analysis/", conditionNames[1])
data_pattern <- c( "_48000.txt", "_32000.txt", "_64000.txt")
save_pattern <- c("_48000", "_32000", "_64000")
model_type <- "m2cov"
mxSetup <- "testJitter_VF-.1_a.1"
n_models <- "30by20"

if (!dir.exists(save_path)){
    dir.create(save_path)}

# save the summary list
saveRDS(summary_list, paste0(save_path, "/", model_type,mxSetup, save_pattern[1], "_nModel", n_models, "_summary_list.rds"))