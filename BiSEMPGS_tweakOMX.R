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
    
    mxOption(NULL,"Feasibility tolerance","1e-5")
    mxOption(NULL,"Number of Threads",parallel::detectCores())
    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
# Load the simulated data for this demonstration:
# Load the simulated data for this demonstration:
    Example_Data  <- fread("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/testdata/loop103.rds_16000.txt", header = T)
    Example_Data2 <- fread("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/testdata/loop104.rds_16000.txt", header = T)
    Example_Data3 <- fread("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/testdata/loop105.rds_16000.txt", header = T)
    Example_Data4 <- fread("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/testdata/loop106.rds_16000.txt", header = T)
    Example_Data  <- rbind(Example_Data, Example_Data2, Example_Data3, Example_Data4)

    cov(Example_Data, use="pairwise.complete.obs")

# change the type of the first column of the df

# Create variables and define the algebra for each variables

 
    VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.2,.2,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
	VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.3,0.3,.2), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.05) # Variance due to VT
    VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.7,.1,0.1,.5), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

    VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")

    VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
    VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
# Genetic effects:
    delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
    a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.2,.1), label=c("a11", "a22"),    name="a", lbound = -.05)     # Effect of latent PGS on phen
    k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.04,0.04,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.05)     # PGS variance (if no AM)
    j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.1,0.1,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.05)     # Latent PGS variance (if no AM)
    Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.05) # Within-person PGS-Phen covariance
    Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.15,0.1,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.05) # Within-person latent PGS-Phen covariance

    Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
    Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

	Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
# Assortative mating effects:
    mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0.15,0.1,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.05) # AM co-path coefficient
    gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.15,0.1,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.05)  # Increase in cross-mate PGS (co)variances from AM
    ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.15,0.1,.05), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.05)  # Increase in cross-mate latent PGS (co)variances from AM
    gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.05)  # Increase in within-mate PGS (co)variances from AM
    hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.17,0.05,0.05,.15), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.05)  # Increase in within-mate latent PGS (co)variances from AM
    gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
    ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
    gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
    hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
    
    gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta
    #gchc_constraint_Algebra <- mxAlgebra( solve(t(chol(solve(2* delta %*% k %*% t(delta))))) %*% t(chol(solve(2* a %*% j %*% t(a)))) %*% hc %*% chol(solve(2* a %*% j %*% t(a))) %*% solve(chol(solve(2* delta %*% k %*% t(delta)))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

    itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.15,0.1,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.05) 
    itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.15,0.1,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.05)
    ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.1,0.1,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.05) 
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
    f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.15,0.1,.08), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.05) # Vertical Transmission effect
    w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.05) # Genetic nurture
    v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.15,.1), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.05) # Latent nurture
    w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
    v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
    wv_constraint_algebra <- mxAlgebra((w * (2*delta%*%k%*%t(delta))/(2*a%*%j%*%t(a))), name='wv_constraint_algebra')
    #wv_constraint_algebra <- mxAlgebra(solve(t(chol(solve(2* delta %*% k %*% t(delta))))) %*% t(chol(solve(2* a %*% j %*% t(a)))) %*% v %*% chol(solve(2* a %*% j %*% t(a))) %*% solve(chol(solve(2* delta %*% k %*% t(delta)))), name = "wv_constraint_algebra") # w and v are equally proportional to a and delta
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

# Create fit function:
    FitFunctionML <- mxFitFunctionML()
    # Create fit function:
    #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
# Specify what parameters we're going to be including in our model:
    Params <- list(
                VY, VF, VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v, 
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                VY_Constraint, 
                VF_Constraint, 
                #Omega_Constraint, 
                Gamma_Constraint, 
                #gt_constraint, 
                ht_constraint, 
                #gc_constraint, 
                #hc_constraint, 
                gchc_constraint, 
                itlo_constraint, 
                itol_constraint, 
                #ic_constraint, 
                v_constraint, 
                w_constraint,
                #wv_constraint,
                thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                CovMatrix, Means, ModelExpectations, FitFunctionML)
# Create the model:
    options(warning.length = 8000)
    Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)
    #fitModel1 <- mxRun(Model1, intervals=T, silent=F)
    fitModel1 <- mxTryHard(Model1, extraTries = 5, intervals=T, silent=F)

    summary(fitModel1)



# some code check the constraints
a <- matrix(c(sqrt(.49*.5),0,0,sqrt(.16*.7)),nrow=2,byrow=T)
j <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.49*.5),0,0,sqrt(.16*.3)),nrow=2,byrow=T)
k <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
gt <- matrix(c(0.042433374,0.012050565,0.007746527,0.006565526),nrow=2,byrow=T)
gc <- (gt + t(gt))/2
print(gc)
ht <- matrix(c(0.043228259,0.01599563,0.009422301,0.01348196),nrow=2,byrow=T)
hc <- (ht + t(ht))/2
gc/(2* delta %*% k %*% t(delta))
hc/(2* a %*% j %*% t(a))
print(gc)

#a = delta
w <- matrix(c(0.20374558,0.07627376,0.08259022,0.05076416),nrow=2,byrow=T)
w/(2* delta %*% k %*% t(delta))
v <- matrix(c(0.20766181,0.10233118,0.08560629,0.07192283),nrow=2,byrow=T)
v/(2* a %*% j %*% t(a))




# some code check the constraints
a <- matrix(c(sqrt(.8*.2),0,0,sqrt(.2*.7)),nrow=2,byrow=T)
j <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.8*.8),0,0,sqrt(.2*.3)),nrow=2,byrow=T)
k <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
gt <- matrix(c(0.1040189,0.01429702,0.02546207,0.01548387),nrow=2,byrow=T)
gc <- (gt + t(gt))/2
print(gc)

gc/hc
ht <- matrix(c(0.02739692,0.01081928,0.01081928,0.03019162),nrow=2,byrow=T)
hc <- (ht + t(ht))/2
print(hc)
#a = delta
gc/(2* delta %*% k %*% t(delta))
hc/(2* a %*% j %*% t(a))
gc * solve(2* delta %*% k %*% t(delta))
hc * solve(2* a %*% j %*% t(a))
t(chol(solve(2* delta %*% k %*% t(delta)))) %*% gc %*% chol(solve(2* delta %*% k %*% t(delta)))
t(chol(solve(2* a %*% j %*% t(a)))) %*% hc %*% chol(solve(2* a %*% j %*% t(a)))
#print(gc)

#a = delta
w <- matrix(c(0.3707127,0.13154873,0.1488213,0.08749667),nrow=2,byrow=T)
w/(2* delta %*% k %*% t(delta))
w %*% solve(2* delta %*% k %*% t(delta))
v <- matrix(c(0.19777621,0.1516655,0.08412774,0.1142102),nrow=2,byrow=T)
v/(2* a %*% j %*% t(a))
v %*% solve(2* a %*% j %*% t(a))
t(chol(solve(sqrt(2* delta %*% k %*% t(delta))))) %*% (w) %*% chol(solve(sqrt(2* delta %*% k %*% t(delta))))
t(chol(solve(sqrt(2* a %*% j %*% t(a))))) %*% (v) %*% chol(solve(sqrt(2* a %*% j %*% t(a))))





# use pairs of delta and a to solve for alpha, beta and tau
library(matlib)
#1
a <- matrix(c(sqrt(.8*.2),0,0,sqrt(.2*.7)),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.8*.8),0,0,sqrt(.2*.3)),nrow=2,byrow=T)

g11 <- 0.1038053
g12 <- 0.01982926
g22 <- 0.009278857

h11 <- 0.02724393
h12 <- 0.009351237
h22 <- 0.01543556

x11 = delta[1,1]^2/g11 - a[1,1]^2/h11
x12 = delta[1,1]*delta[2,2]/g11 - a[1,1]*a[2,2]/h11
x13 = delta[2,2]^2/g11 - a[2,2]^2/h11

y11 = delta[1,1]*delta[2,2]/g11
y12 = -a[1,1]*a[2,2]/h11

b1 = a[1,1]^2/h11 - delta[1,1]^2/g11


#2
a <- matrix(c(sqrt(.6*.3),0,0,sqrt(.7*.4)),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.6*.7),0,0,sqrt(.7*.6)),nrow=2,byrow=T)

g11 <- 0.07267078
g12 <- 0.0172996
g22 <- 0.04569342


h11 <- 0.03141719
h12 <- 0.009036486
h22 <- 0.03007751

g11/(delta[1,1]^2+.5*.1*delta[1,1]*delta[2,2])
h11/(a[1,1]^2+.5*.1*a[1,1]*a[2,2])

x21 = delta[1,1]^2/g11 - a[1,1]^2/h11
x22 = delta[1,1]*delta[2,2]/g11 - a[1,1]*a[2,2]/h11
x23 = delta[2,2]^2/g11 - a[2,2]^2/h11

y21 = delta[1,1]*delta[2,2]/g11
y22 = -a[1,1]*a[2,2]/h11
b2 = a[1,1]^2/h11 - delta[1,1]^2/g11

A = matrix(c(y11,y12,y21,y22),nrow=2,byrow=T)
b = c(b1,b2)
showEqn(A,b)
Solve(A,b)

#3 
a <- matrix(c(sqrt(.3*.2),0,0,sqrt(.4*.3)),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.3*.8),0,0,sqrt(.4*.7)),nrow=2,byrow=T)

g11 <- 0.04209733
g12 <- 0.01624357
g22 <- 0.03215703

h11 <- 0.01078906
h12 <- 0.005306831
h22 <- 0.01343727

x31 = delta[1,1]^2/g11 - a[1,1]^2/h11
x32 = delta[1,1]*delta[2,2]/g11 - a[1,1]*a[2,2]/h11
x33 = delta[2,2]^2/g11 - a[2,2]^2/h11


# try to solve for alpha, beta and tau
A <- matrix(c(x11,x12,x13,x21,x22,x23,x31,x32,x33),nrow=3,byrow=T)
b <- c(0,0,0)
showEqn(A,b)
plotEqn3d(A,b)
Solve(A,b)

# some code check the constraints
j <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
k <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
gc <- matrix(c(g11,g12,g12,g22),nrow=2,byrow=T)
hc <- matrix(c(h11,h12,h12,h22),nrow=2,byrow=T)
#a = delta
gc/(2* delta %*% k %*% t(delta))
hc/(2* a %*% j %*% t(a))
gc * solve(2* delta %*% k %*% t(delta))
hc * solve(2* a %*% j %*% t(a))
t(chol(solve(2* delta %*% k %*% t(delta)))) %*% gc %*% chol(solve(2* delta %*% k %*% t(delta)))
t(chol(solve(2* a %*% j %*% t(a)))) %*% hc %*% chol(solve(2* a %*% j %*% t(a)))



#4
a <- matrix(c(sqrt(.5*.1),0,0,sqrt(.6*.5)),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.5*.9),0,0,sqrt(.6*.5)),nrow=2,byrow=T)

g11 <- 0.07815197
g12 <- 0.01936456
g22 <- 0.03428412

h11 <- 0.00955794
h12 <- 0.00697541
h22 <- 0.03200986

# some code check the constraints
j <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
k <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
gc <- matrix(c(g11,g12,g12,g22),nrow=2,byrow=T)
hc <- matrix(c(h11,h12,h12,h22),nrow=2,byrow=T)
#a = delta
gc/(2* delta %*% k %*% t(delta))
hc/(2* a %*% j %*% t(a))
gc * solve(2* delta %*% k %*% t(delta))
hc * solve(2* a %*% j %*% t(a))
t(chol(solve(2* delta %*% k %*% t(delta)))) %*% gc %*% chol(solve(2* delta %*% k %*% t(delta)))
t(chol(solve(2* a %*% j %*% t(a)))) %*% hc %*% chol(solve(2* a %*% j %*% t(a)))

#5
a <- matrix(c(sqrt(.4*.4),0,0,sqrt(.75*.6)),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.4*.6),0,0,sqrt(.75*.4)),nrow=2,byrow=T)

g11 <- 0.04196104
g12 <- 0.01209637
g22 <- 0.03269869

h11 <- 0.02870104
h12 <- 0.01234879
h22 <- 0.04802063

y11 = delta[1,1]*delta[2,2]/g11
y12 = -a[1,1]*a[2,2]/h11

b1 = a[1,1]^2/h11 - delta[1,1]^2/g11

# some code check the constraints
j <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
k <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
gc <- matrix(c(g11,g12,g12,g22),nrow=2,byrow=T)
hc <- matrix(c(h11,h12,h12,h22),nrow=2,byrow=T)
#a = delta
gc/(2* delta %*% k %*% t(delta))
hc/(2* a %*% j %*% t(a))
gc * solve(2* delta %*% k %*% t(delta))
hc * solve(2* a %*% j %*% t(a))
t(chol(solve(2* delta %*% k %*% t(delta)))) %*% gc %*% chol(solve(2* delta %*% k %*% t(delta)))
t(chol(solve(2* a %*% j %*% t(a)))) %*% hc %*% chol(solve(2* a %*% j %*% t(a)))

#6

a <- matrix(c(sqrt(.35*.65),0,0,sqrt(.2*.8)),nrow=2,byrow=T)
delta <- matrix(c(sqrt(.35*.35),0,0,sqrt(.2*.2)),nrow=2,byrow=T)

.05 * delta[1,1] * delta[2,2]

g11 <- 0.02121644
g12 <- 0.005760222
g22 <- 0.005132298

h11 <- 0.040252
h12 <- 0.01401233
h22 <- 0.018841

y21 = delta[1,1]*delta[2,2]/g11
y22 = -a[1,1]*a[2,2]/h11
b2 = a[1,1]^2/h11 - delta[1,1]^2/g11

A = matrix(c(y11,y12,y21,y22),nrow=2,byrow=T)
b = c(b1,b2)
showEqn(A,b)
Solve(A,b)

# some code check the constraints
j <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
k <- matrix(c(.5,.05,.05,.5),nrow=2,byrow=T)
gc <- matrix(c(g11,g12,g12,g22),nrow=2,byrow=T)
hc <- matrix(c(h11,h12,h12,h22),nrow=2,byrow=T)
#a = delta
gc/(2* delta %*% k %*% t(delta))
hc/(2* a %*% j %*% t(a))
gc %*% solve(2* delta %*% k %*% t(delta))
hc %*% solve(2* a %*% j %*% t(a))
t(chol(solve(2* delta %*% k %*% t(delta)))) %*% gc %*% chol(solve(2* delta %*% k %*% t(delta)))
t(chol(solve(2* a %*% j %*% t(a)))) %*% hc %*% chol(solve(2* a %*% j %*% t(a)))
