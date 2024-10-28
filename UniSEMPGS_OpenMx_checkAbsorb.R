# the script for fitting the OpenMx model with an univariate expected covariance matrix from the Bivariate Iterative Math
# adapted from Yongkong's function perform_SEM_model2_eq
# author: Xuanyu Lyu

###################################################################################
############ the openmx function
###################################################################################
fitCheckObs <- function(cov, max.cores = 2){
    library(OpenMx)
    library(data.table)
    library(stringr)

  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Optimality tolerance","1e-7")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- cov
  #dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel)]


  #Colnames
  #colnames(dat_SEM) <- c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")
  nv=1

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS_lat)=.5 at t0
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.18,label="F",name="f",lbound=-.05)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.6,label="EnvPath",name="e",lbound=-.05)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.01,label="hap_PRS_cov",name="g")
  h <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.01,label="latent_hap_PRS_cov",name="h")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.64),label="obs_coef",name="delta",lbound=-.05)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.4,label="latent_coef",name="a",lbound=-.05)


  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0.1,label="LatentF1",name="x1")
  w1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.15,label="covAandF",name="w1")
  v1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.2,label="covAand_lF",name="v1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.05,label="copath",name="mu1")
  i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.02,label="covOL",name="i1")
  sigma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=1.5,label="sigma",name="sigma1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.3,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.2, label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*sigma1*(1 + sigma1*mu1),name="x2")
  w2 <- mxAlgebra(2*f*Omega1*(1 + sigma1*mu1),name="w2")
  v2 <- mxAlgebra(2*f*Gamma1*(1 + sigma1*mu1),name="v2")
  i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(2*a*i1 + delta*k + 2*delta*g+.5*w1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra(2*a*h + a*j + 2*delta*i1 + .5*v1,name="Gamma2")
  #h <- mxAlgebra((g*a^2)/delta^2,name="h")
  h2 <- mxAlgebra(Gamma1^2*mu1,name="h2")


  # Equating nonlinear constraints and parameters
  xCon <- mxConstraint(x1==x2,name="xCon")
  wCon <- mxConstraint(w1==w2,name="wCon")
  vCon <- mxConstraint(v1==v2,name="vCon")
  iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  hCon <- mxConstraint(h==h2,name="hCon")
  sigmaCon <- mxConstraint(sigma1==sigma2,name="sigmaCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")


  # mxAlgebra for implied parameters and relative covariances
  sigma2 <- mxAlgebra(2*a*Gamma1 + 2*delta*Omega1 + a*v1 + delta*w1 + 2*f^2*sigma1*(1 + sigma1*mu1) + e^2,name="sigma2")
  ThetaNTM05 <- mxAlgebra((Omega1-delta*k),name="ThetaNTM05") #this is HALF theta_NTM in the math
  spsPRS <-mxAlgebra(sigma1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS
  PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigma1)*(1+mu1*sigma1),name="PO")
  spouse <-mxAlgebra(mu1*sigma1^2,name="spouse")

  v_Mean <- rep(0, 7)  
  names(v_Mean) <- c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp          Ym         Yp        Yo
    cbind(  k+g        ,g         ,g          ,g        ,Omega1     ,spsPRS    ,ThetaNTM05   ),    #NTm
    cbind(  g          ,k+g       ,g          ,g        ,Omega1     ,spsPRS    ,Omega1        ),    #Tm
    cbind(  g          ,g         ,k+g        ,g        ,spsPRS     ,Omega1    ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,g          ,k+g      ,spsPRS     ,Omega1    ,Omega1     ),    #Tp
    cbind(  Omega1     ,Omega1    ,spsPRS     ,spsPRS   ,sigma1     ,spouse    ,PO         ),    #Ym
    cbind(  spsPRS     ,spsPRS    ,Omega1     ,Omega1   ,spouse     ,sigma1    ,PO         ),    #Yp
    cbind(  ThetaNTM05 ,Omega1    ,ThetaNTM05 ,Omega1   ,PO         ,PO        ,sigma1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=7, free=TRUE,  values= rep(0.1,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYm","meanYp","meanYo"),dimnames=list(NULL,c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1"))

  params <- list(k ,j,  f , e , g , 
                 h, delta,a,
                 #x1 , 
                 w1 ,mu1 , i1, sigma1, Omega1 , Gamma1,
                 #x2,
                 w2, v1,v2, i2, g2,
                 h2,
                 Omega2 , Gamma2, sigma2,
                 #xCon, 
                 wCon, 
                 iCon, 
                 gCon, 
                 OmegaCon, 
                 GammaCon, 
                 vCon, 
                 hCon, 
                 sigmaCon,
                 ThetaNTM05, 
                 spsPRS, PO, spouse,
                 CVmat, MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="cov", numObs = 4.8e4, means=v_Mean)
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  AFE.Fit=mxRun(modelAFE,intervals=TRUE,silent=FALSE)
  
  # uncomment to print the openmx summary
  #return(summary(AFE.Fit))

  if(class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
    VAL=mxEval(a^2*2*(j+2*h), AFE.Fit$AFEmodel, T)
    VF=mxEval(x1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigma=mxEval(sigma1,AFE.Fit$AFEmodel,T)
    west=mxEval(w1,AFE.Fit$AFEmodel,T)
    vest=mxEval(v1,AFE.Fit$AFEmodel,T)
    gest=mxEval(g2,AFE.Fit$AFEmodel,T)
    hest=mxEval(h,AFE.Fit$AFEmodel,T)
    h2=(VAO+VAL)/sigma
    VFratio=VF/sigma
    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=mxEval(a,AFE.Fit$AFEmodel,T)
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results44=data.frame(method="model2_eq",VAO=VAO,VAL=VAL,VF=VF,VE=VE,VY=sigma,w=west,v=vest,g=gest,h=hest, h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)
  }else{
    results44=data.frame(method="model2_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,vest=NA,gest=NA,hest=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }
  return(results44)
}

###################################################################################
#### the iterative math to get the iterative covariance matrix
###################################################################################
startingParamList1 <- list(vg1 = .64,
						   vg2 = .49,
						   am11 = 0,
						   am12 = 0,
						   am21 = 0,
						   am22 = 0,
						   f11 = 0,
						   f12 = 0.5,
						   f21 = 0.1,
						   f22 = 0,
						   Nfam = 5e4,
						   rg = 0,
						   re = 0,
						   prop.h2.latent1 = 0,
						   prop.h2.latent2 = 0)


    condition =1
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

        (exp.ic[[it]] <- .5*(exp.itol[[it]] + t(exp.itlo[[it]])))
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

### create the expected covariance matrix
CMatrix <- rbind(
    #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
    cbind(VY,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
    cbind(Ym_Yp,     VY,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
    cbind(Yo_Yp,     Yo_Ym,     VY,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
    cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
    cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
    cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
    cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc))  #NTm1 NTm2

#write.csv(CMatrix, "ExpectedCMatrix.csv")
# compare with simulated covariance matrix
#CMatrix_obs <- read.csv("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/Full_Model/loop123.rds_64000.txt", sep="\t", header=TRUE) |> cov()
#CMatrix <- cov(read.csv("testDataFull.txt", sep="\t", header=TRUE))
# change the column names to what we defined in the OpenMx script
#colnames(data_df) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
#CMatrix <- read.csv("testDataFull.txt", sep="\t", header=TRUE)



colnames(CMatrix) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")
rownames(CMatrix) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")

# only get parts of the covariance matrix for trait 1
CMatrix_uni <- CMatrix[c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1"), c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")]

# switch the order and the rows of the covariance matrix to the order of the OpenMx script
CMatrix_uni <- CMatrix_uni[c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1"), c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")]



### check the estimates
print(fitCheckObs(CMatrix_uni))
