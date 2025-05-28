# A script to test the different delta simulation and model fitting
source("PaperScripts/Simulate.Multivariate.NOLD.AM.FUNCTIONS.differentDelta-xuanyu2.R")

save_dir_data <- "Data/Paper/differentDelta/"
if (!dir.exists(save_dir_data)) {
    dir.create(save_dir_data, recursive = TRUE)
}
save_dir_summary <- "Analysis/Paper/differentDelta/"
if (!dir.exists(save_dir_summary)) {
    dir.create(save_dir_summary, recursive = TRUE)
}
# simulation conditions setup
conditionNames <- c("onlyAM", "onlyVT", "bothAMVT", "bothAMVT_StrongCrossTrait", "OnlyUni")

startingParamList1 <- list(vg1 = rep(.64,5),
						   vg2 = rep(.36,5),
						   am11 = c(rep(0.4,3),.05,.4),
						   am12 = c(.2,0,.2,.4,0),
						   am21 = c(.2,0,.2,.4,0), 
						   am22 = c(rep(0.3,3),.1,.2),
						   f11 = c(rep(0.15,3),0,.3),
						   f12 = c(0,.1,.1,.3,0),
						   f21 = c(0,.05,.05,.25,0), 
						   f22 = c(rep(0.1,4),.1),
						   Nfam = rep(3e4, 5),
						   rg = rep(0,5),
						   re = rep(0,5),
						   prop.h2.latent1 = rep(.60/.64,5),
						   prop.h2.latent2 = rep(.8,5),
                           vg1_last = rep(.49,5),
                           vg2_last = rep(.36,5),
                           rg_last = rep(0,5),
                           prop.h2.latent1_last = rep(.3/.49,5),
                           prop.h2.latent2_last = rep(.8,5))

for (condition in 5:5){
	# WILDCARD parameters
	pop.size <-  startingParamList1["Nfam"][[1]][[condition]] #maybe something like 2e4, or 20000, when running for real
	num.cvs <- 50 #maybe 25
	seed <- 62*condition
	num.gen <-  15 #8 should be sufficient to get to AM equilibrium
	avoid.inb <- TRUE #avoid inbreeding?
	save.covariances <- TRUE #save the covariance matrices?
	save.history <- FALSE #save the data for each generation?
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
    # VG at last generation
    vg1_last <- startingParamList1["vg1_last"][[1]][[condition]] #trait 1 vg at last generation
    vg2_last <- startingParamList1["vg2_last"][[1]][[condition]] #trait 2 vg at last generation
    rg_last  <- startingParamList1["rg_last"][[1]][[condition]] 
    k2.matrix_last <- matrix(c(1,rg_last,rg_last,1),nrow=2,byrow=T) #k2 matrix is 2 * k matrix - i.e., genotypic (instead of haplotypic) var/covar at t0
    prop.h2.latent1_last <- startingParamList1["prop.h2.latent1_last"][[1]][[condition]] #  trait 1, e.g., height at last generation
    prop.h2.latent2_last <- startingParamList1["prop.h2.latent2_last"][[1]][[condition]] # trait 2, e.g., IQ at last generation


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
    #Create lower delta matrix for the last generation
    (vg.obs1_last <- vg1_last*(1-prop.h2.latent1_last))
    (vg.obs2_last <- vg2_last*(1-prop.h2.latent2_last))
    (covg.obs_last <- rg_last*sqrt(vg.obs1_last)*sqrt(vg.obs2_last))
    d11_last <- sqrt(vg.obs1_last)
    #d21_last <- covg.obs_last/d11_last
    d21_last <- 0 #new way
    d22_last <- sqrt(vg.obs2_last - d21_last^2)
    (delta.mat_last <- matrix(c(d11_last,0,d21_last,d22_last),byrow=T,nrow=2))
    (vg.obs_last <- delta.mat_last %*% k2.matrix_last %*% t(delta.mat_last)) #Check

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
    #Create lower a matrix for the last generation
    (vg.lat1_last <- vg1_last*(prop.h2.latent1_last))
    (vg.lat2_last <- vg2_last*(prop.h2.latent2_last))
    (covg.lat_last <- rg_last*sqrt(vg.lat1_last)*sqrt(vg.lat2_last))
    a11_last <- sqrt(vg.lat1_last)
    #a21_last <- covg.lat_last/a11_last
    a21_last <- 0 #new way
    a22_last <- sqrt(vg.lat2_last - a21_last^2)
    (a.mat_last <- matrix(c(a11_last,0,a21_last,a22_last),byrow=T,nrow=2))
    (vg.lat_last <- a.mat_last %*% k2.matrix_last %*% t(a.mat_last)) #Check

	#Find the total genetic covariance
	#(covg <- vg.obs[1,2] + vg.lat[1,2])
	(covg.mat <- vg.obs + vg.lat)
	a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) #check
    
    # covg for the last generation
    (covg.mat_last <- vg.obs_last + vg.lat_last)

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

    # last-gen alpha and causal variants
    seed <- 62*condition+1
    set.seed(seed)
    alphas.pre_last <- mvrnorm(num.cvs,c(0,0),matrix(c(1,rg,rg,1),nrow=2),empirical=TRUE)
    alphas_last <- alphas.pre_last * cbind(sqrt(1/((num.cvs-1)*gentp.var)),sqrt(1/((num.cvs-1)*gentp.var)))
    cv.info_last <- data.frame(maf=maf.vector,alpha1=alphas_last[,1],alpha2=alphas_last[,2]) #we'll use this for both the observed and latent

	ObjectsKeep <- as.character(ls())
	# write a loop to run the simulation 100 times and save all the summary data in a list
	#l.summaryLast <- list()
	#l.all <- list()
	for (i in 1:1){
	loop_index <- i
	# if the file is already there, skip
	if (file.exists(paste0(save_dir_data, conditionNames[condition],"/loop40k",loop_index,".rds"))){
		cat(conditionNames[condition],"/Simulation",loop_index,"already exists\n")
		next
	}
	#print the system time
	print(Sys.time())
	time <- Sys.time()
	AM.DATA <- AM.SIMULATE(CV.INFO=cv.info, NUM.GENERATIONS=num.gen, POP.SIZE=pop.size, AVOID.INB=avoid.inb, SAVE.EACH.GEN=save.history, SAVE.COVS=save.covariances, SEED=seed*loop_index, 
						cove.mat=cove.mat, fmat=f.mat, amat=a.mat, dmat=delta.mat, cor.list=am.list, covy=COVY, k2.matrix=k2.matrix, 
                        last_gen = TRUE, CV.INFO_last = cv.info_last, dmat_last = delta.mat_last, amat_last = a.mat_last,
                        k2.matrix_last = k2.matrix_last)
	SUMMARY.last <- AM.DATA$SUMMARY.RES[[num.gen]]
	# print the time it takes to run the simulation
	cat("Time to run simulation: ", Sys.time()-time, "\n")

	#l.summaryLast[[i]] <- SUMMARY.last
	#l.all[[i]] <- AM.DATA
	# test if a folder exist, if not, create one
	if (!dir.exists(paste0(save_dir_summary, conditionNames[condition]))){
		dir.create(paste0(save_dir_summary, conditionNames[condition]))
	}		
	if (!dir.exists(paste0(save_dir_data,conditionNames[condition]))){
		dir.create(paste0(save_dir_data, conditionNames[condition]))
	}
	# save the data
	saveRDS(SUMMARY.last, file=paste0(save_dir_summary, conditionNames[condition],"/loop40k",loop_index,".rds"))
	saveRDS(AM.DATA, file=paste0(save_dir_data, conditionNames[condition],"/loop",loop_index,".rds"))
	cat(conditionNames[condition],"/Simulation",loop_index,"done\n")
	rm(list = setdiff(ls(), c(ObjectsKeep, "i", "ObjectsKeep")))
}
}

# load a set and check the variance
# load the data
test_path <- "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/R-projects/BiSEMPGS/Data/Paper/differentDelta/OnlyUni"
test <- readRDS(test_path)
test$SUMMARY.RES[[16]]

data_prep <- function(data_path, save_path, target_n){
    # get the number of .rds files in the folder
    l_files <- list.files(data_path, pattern = ".rds")
    # read all the .rds file one by one and save them as .txt files
    for (i in 1:length(l_files)){
        
        # check if the data has been converted to txt
        if (file.exists(paste0(save_path, "/", l_files[i],"_",target_n, ".txt"))){
            cat("Data file", l_files[i], "has been converted to txt\n")
            next
        }
        # load the data
        data <- readRDS(paste0(data_path, "/", l_files[i]))

        data_df <- data$PHEN

        data_df <- data_df[,c("ID", "Father.ID", "Mother.ID",
                              "Y1P","Y2P","Y1M","Y2M","Y1","Y2",
                              "TPO1","TPO2","NTPO1","NTPO2",
                              "TMO1","TMO2","NTMO1","NTMO2")] |> as.data.frame()
        
        # remove the rows with the same Father.ID
        data_df <- data_df[!duplicated(data_df$Father.ID),]
        print(nrow(data_df))
        # remove the columns with ID
        data_df <- data_df[,-(1:3)]
        
        # change the column names to what we defined in the OpenMx script
        colnames(data_df) <- c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")

        # sample the data to the target_n
        data_df <- data_df[sample(1:nrow(data_df), target_n),]
        # write the data as a txt file
        if (!dir.exists(save_path)){
		    dir.create(save_path)
	    }
        write.table(data_df, file = paste0(save_path, "/", l_files[i],"_",target_n, ".txt"), sep = "\t", row.names = FALSE)
        cat("Data file", l_files[i], "has been saved as a txt file\n")
        #return(data_df)
    }
}

data_prep(test_path, "Data/Paper/differentDelta/OnlyUni", 12000)

#source("PaperScripts/10-lastDiffOpenMxFunctions.R")

#test_txt_path <- "PaperScripts/loop1.rds_4000.txt"
#test_summary <- fitBiSEMPGS_m2_tol(test_txt_path)
#test_summary
