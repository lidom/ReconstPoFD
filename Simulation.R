## R-packages 
library("devtools")
# install necessary packages:
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("fdapace")

library("ReconstPoFD")  # contains the function 'reconstuct()'


# library("doParallel")   # parallel-looping
# registerDoParallel(cores=6)
# getDoParWorkers()

## Location to store the simulation results:
setwd("/home/dom/ownCloud/Kneip_Liebl_Reconstruction/Simulation_Submission_2")

## #######################################
## Number of MC-Repetitions
B         <-  100
## #######################################

## #######################################
## Lower-value of the total domain 
a             <-   0             
## Upper-value of the total domain
b             <-   1
## Regular grid points
nRegGrid      <-  51 
## 
##
#determ_obs_interv <- c((a+(b-a)*0.35), (b-(b-a)*0.35))
determ_obs_interv <- c(0,.5)

## #######################################

(Start.Time <- Sys.time())


winsorize_x <- function(x, cut = 0.01){
  if(!any(is.na(x))){
  cut_point_top    <- quantile(x, 1 - cut)
  cut_point_bottom <- quantile(x,     cut)
  i <-  which(x >= cut_point_top);    x[i] <- cut_point_top
  j <-  which(x <= cut_point_bottom); x[j] <- cut_point_bottom
  }
  return(x)
}


# for(DGP in c('DGP1','DGP4','DGP5')[1]){
for(DGP in c('DGP1','DGP2','DGP3')[2]){
  ## Set seed
  ##
  for(n in c(50, 80)){
    # if(any(DGP==c('DGP1','DGP2'))){m_seq <- c(15, 50)}else{m_seq <- NA}
    if(DGP=='DGP1'){m_seq <- c(15, 50)}else{m_seq <- NA}
    for(m in m_seq){
      
      ## a <- 0; b <- 1; DGP <- 'DGP1'; n <- 50; m <- 15; nRegGrid <- 51; B <- 1
      ## a <- 0; b <- 1; DGP <- 'DGP2'; n <- 50; m <- NA; nRegGrid <- 51; B <- 10
      ## a <- 0; b <- 1; DGP <- 'DGP3'; n <- 50; m <- NA; nRegGrid <- 51; B <- 3
      
      ## #######################################################################
      cat(DGP,"n=",n,"m=",m,"\n")
      ## #######################################################################
      ## #######################################################################
      ## Preselecting partially observed *target* functions to be reconstructed
      ## #######################################################################
      ##
      set.seed(123)
      ##
      if(DGP=='DGP1'){SimDat <- ReconstPoFD::simuldata(n = 1, m = m, a = a, b = b, DGP='DGP1', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP2'){SimDat <- ReconstPoFD::simuldata(n = 1, m = m, a = a, b = b, DGP='DGP2', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP3'){SimDat <- ReconstPoFD::simuldata(n = 1, m = m, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      ##
      Y_target_true_mat  <- SimDat[['Y_true_mat']]
      U_target_true_mat  <- SimDat[['U_true_mat']]
      ##
      Y_target_mat       <- SimDat[['Y_mat']]
      U_target_mat       <- SimDat[['U_mat']]
      Y_target_list      <- SimDat[['Y_list']]
      U_target_list      <- SimDat[['U_list']]
      ##
      A_target_vec       <- SimDat[['A_vec']]
      B_target_vec       <- SimDat[['B_vec']]
      ##
      missings_of_target <- c(U_target_true_mat[,1] < A_target_vec[1] | U_target_true_mat[,1] > B_target_vec[1])
      ##
      Y_PS_FALSE_MC_mat  <- matrix(NA, nrow = nRegGrid, ncol = B)
      Y_PS_TRUE_MC_mat   <- matrix(NA, nrow = nRegGrid, ncol = B)
      Y_CEScores_MC_mat  <- matrix(NA, nrow = nRegGrid, ncol = B)
      Y_PACE_MC_mat      <- matrix(NA, nrow = nRegGrid, ncol = B)
      Y_Kraus_MC_mat     <- matrix(NA, nrow = nRegGrid, ncol = B)
      ##
      ## #######################################################################
      for(repet in 1:B){ # repet <- 1
        ##
        if(DGP=='DGP1'){SimDat <- ReconstPoFD::simuldata(n = n-1, m = m, a = a, b = b, DGP='DGP1', nRegGrid = nRegGrid)}
        if(DGP=='DGP2'){SimDat <- ReconstPoFD::simuldata(n = n-1, m = m, a = a, b = b, DGP='DGP2', nRegGrid = nRegGrid)}
        if(DGP=='DGP3'){SimDat <- ReconstPoFD::simuldata(n = n-1, m = m, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid)}
        ##
        Y_mat      <- cbind(Y_target_mat, SimDat[['Y_mat']])
        U_mat      <- cbind(U_target_mat, SimDat[['U_mat']])
        ##
        Y_list     <- c(Y_target_list,    SimDat[['Y_list']])
        U_list     <- c(U_target_list,    SimDat[['U_list']])
        ##
        ## Reconstruction Operator 'without Pre-Smoothing'
        result_PS_FALSE <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                    Lu           = U_list,
                                                    K            = NULL,
                                                    method       = "PS_FALSE",
                                                    reconst_fcts = 1, 
                                                    nRegGrid     = nRegGrid)
        Y_PS_FALSE_mat <- matrix(unlist(result_PS_FALSE[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
        Y_PS_FALSE_MC_mat[,repet] <- Y_PS_FALSE_mat
        ##
        ## Reconstruction Operator 'with Pre-Smoothing'
        result_PS_TRUE <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                   Lu           = U_list,
                                                   K            = NULL,
                                                   method       = "PS_TRUE",
                                                   reconst_fcts = 1,
                                                   nRegGrid     = nRegGrid)
        Y_PS_TRUE_mat <- matrix(unlist(result_PS_TRUE[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
        Y_PS_TRUE_MC_mat[,repet] <- Y_PS_TRUE_mat
        ## 
        ## Reconstruction Operator 'with PACE-scores'
        if(DGP=='DGP1'){
          result_CEScores <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                      Lu           = U_list,
                                                      K            = NULL,
                                                      method       = "CEScores",
                                                      BwMu         = h.mu,
                                                      BwCov        = h.cov,
                                                      reconst_fcts = 1, 
                                                      nRegGrid     = nRegGrid)
          Y_CEScores_mat <- matrix(unlist(result_CEScores[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
          Y_CEScores_MC_mat[,repet] <- Y_CEScores_mat
        }
        ##
        ## PACE of Yao, Müller & Wang (2005, JASA) 
        if(DGP=='DGP1'){
        if(DGP=='DGP1'){
          opt_list <- list(
            "dataType"       = "Sparse", "kernel"         = "gauss",
            "methodMuCovEst" = "smooth", "error"          = TRUE,
            "userBwMu"       = h.mu,     "userBwCov"      = h.cov,
            "nRegGrid"       = nRegGrid)
        }else{
          opt_list <- list(
            "dataType"       = "Dense",  "kernel"         = "gauss",
            "methodMuCovEst" = "smooth", "error"          = FALSE,    
            "userBwMu"       = h.mu,     "userBwCov"      = h.cov,
            "nRegGrid"       = nRegGrid) 
        }
        result_PACE <- fdapace::FPCA(Ly = Y_list, Lt = U_list, optns = opt_list)
        ##
        Y_PACE_mat <- t(fitted(result_PACE))[,1]
        Y_PACE_MC_mat[,repet] <- Y_PACE_mat
        }
        ##
        if(any(DGP==c('DGP2','DGP3'))){
          ## Reconstruction Operator of Kraus (2015, JRSSB)
          result_Kraus            <- ReconstPoFD::reconstructKraus(X_mat = Y_mat, reconst_fcts = 1)
          Y_Kraus_mat             <- result_Kraus[['X_reconst_mat']]
          Y_Kraus_MC_mat[,repet]  <- Y_Kraus_mat
        }
        ##
        ## ##################################################################
        if(repet %% 100 == 0) cat("repet/B=",repet,"/",B,"\n")
        ## ##################################################################
      } ## End of B-loop
      ##
      slct_M  <- missings_of_target
      ##
      # for(i in 1:30){
      # par(mfrow=c(2,3))
      # ##
      # plot(Y_PS_FALSE_MC_mat[,i], type="b", ylim=range(Y_PS_FALSE_MC_mat[,i],Y_target_true_mat[slct_M,1]),main="PS_FALSE")
      # lines(Y_target_true_mat[,1]); points(y=Y_PS_FALSE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
      # ##
      # plot(Y_PS_TRUE_MC_mat[,i], type="b", ylim=range(Y_PS_TRUE_MC_mat[,i],Y_target_true_mat[slct_M,1]),main="PS_TRUE")
      # lines(Y_target_true_mat[,1]); points(y=Y_PS_TRUE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
      # ##
      # if(DGP=='DGP1'){
      # plot(Y_CEScores_MC_mat[,i], type="b", ylim=range(Y_CEScores_MC_mat[,i],Y_target_true_mat[slct_M,1]),main="CEScores")
      # lines(Y_target_true_mat[,1]); points(y=Y_CEScores_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
      # #}
      # ##
      # plot(Y_PACE_MC_mat[,i], type="b", ylim=range(Y_PACE_MC_mat[,i],Y_target_true_mat[slct_M,1]),main="PACE")
      # lines(Y_target_true_mat[,1]); points(y=Y_PACE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
      # }
      # ##
      # if(any(DGP==c('DGP2','DGP3'))){
      # plot(Y_Kraus_MC_mat[,i], type="b", ylim=range(Y_Kraus_MC_mat[,i],Y_target_true_mat[slct_M,1]),main="Kraus")
      # lines(Y_target_true_mat[,1]); points(y=Y_Kraus_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
      # }
      # par(mfrow=c(1,1))
      # Sys.sleep(1)
      # }
      ##
      cut <- 0 # (if cut=0, no winsorization)
      ##
      PS_FALSE_Mean <- apply(Y_PS_FALSE_MC_mat[slct_M, ], 1, function(x){mean(winsorize_x(x,cut=cut))})
      PS_TRUE_Mean  <- apply(Y_PS_TRUE_MC_mat[ slct_M, ], 1, function(x){mean(winsorize_x(x,cut=cut))})
      CEScores_Mean <- apply(Y_CEScores_MC_mat[slct_M, ], 1, function(x){mean(winsorize_x(x,cut=cut))})
      PACE_Mean     <- apply(Y_PACE_MC_mat[    slct_M, ], 1, function(x){mean(winsorize_x(x,cut=cut))})
      Kraus_Mean    <- apply(Y_Kraus_MC_mat[   slct_M, ], 1, function(x){mean(winsorize_x(x,cut=cut))})
      ##
      PS_FALSE_Int_BiasSq <- sum( (PS_FALSE_Mean - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
      PS_TRUE_Int_BiasSq  <- sum( (PS_TRUE_Mean  - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
      CEScores_Int_BiasSq <- sum( (CEScores_Mean - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
      PACE_Int_BiasSq     <- sum( (PACE_Mean     - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
      Kraus_Int_BiasSq    <- sum( (Kraus_Mean    - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
      ## \int_{Missing} Var(t) dt ('Integrated variance')
      PS_FALSE_Int_Var    <- sum(apply(Y_PS_FALSE_MC_mat[slct_M, ], 1, function(x){var(winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
      PS_TRUE_Int_Var     <- sum(apply(Y_PS_TRUE_MC_mat[ slct_M, ], 1, function(x){var(winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
      CEScores_Int_Var    <- sum(apply(Y_CEScores_MC_mat[slct_M, ], 1, function(x){var(winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
      PACE_Int_Var        <- sum(apply(Y_PACE_MC_mat[    slct_M, ], 1, function(x){var(winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
      Kraus_Int_Var       <- sum(apply(Y_Kraus_MC_mat[   slct_M, ], 1, function(x){var(winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
      #}
      ##
      BiasSq_vec        <- c(PS_FALSE_Int_BiasSq,PS_TRUE_Int_BiasSq,CEScores_Int_BiasSq,PACE_Int_BiasSq,Kraus_Int_BiasSq)
      Var_vec           <- c(PS_FALSE_Int_Var,   PS_TRUE_Int_Var,   CEScores_Int_Var,   PACE_Int_Var,   Kraus_Int_Var)
      names(BiasSq_vec) <- c("PS_FALSE","PS_TRUE","CEScores","PACE","Kraus")
      names(Var_vec)    <- c("PS_FALSE","PS_TRUE","CEScores","PACE","Kraus")
      ##
      ##
      # par(mfrow=c(1,3))
      # barplot(c(BiasSq+Var)[!is.na(BiasSq)]/max(c(BiasSq+Var)[!is.na(BiasSq)]), main="", names.arg = c("NoPS","YesPS","CES","PACE","Kraus")[!is.na(BiasSq)])
      # mtext(text = paste0("MSPE (", round(max(c(BiasSq+Var)[!is.na(BiasSq)]),2),")"), side = 3, line = 1)
      # barplot(c(BiasSq    )[!is.na(BiasSq)]/max(c(BiasSq+Var)[!is.na(BiasSq)]), main="", names.arg = c("NoPS","YesPS","CES","PACE","Kraus")[!is.na(BiasSq)])
      # mtext(text = "Squared Bias", side = 3, line = 1)
      # barplot(c(Var       )[!is.na(BiasSq)]/max(c(BiasSq+Var)[!is.na(BiasSq)]), main="", names.arg = c("NoPS","YesPS","CES","PACE","Kraus")[!is.na(BiasSq)])
      # mtext(text = "Variance", side = 3, line = 1)
      # par(mfrow=c(1,1))      
      ## Save results:
      if(DGP=='DGP1'){
        save(BiasSq_vec, Var_vec, file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
      }
      if(any(DGP==c('DGP2','DGP3'))){
        save(BiasSq_vec, Var_vec, file = paste0(DGP,"_n",n,"_simResults.RData"))
      }
    }
  }
}
##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------

## rm(list=ls())

DGP <- c('DGP1','DGP2','DGP3')[2]
m   <- c(15,  50)[1] 
n   <- c(50,  80)[2] 

## Load results:
if(any(DGP==c('DGP1'))){
  load(file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
}
if(any(DGP==c('DGP2','DGP3'))){
  load(file = paste0(DGP,"_n",n,"_simResults.RData"))
  BiasSq_vec_orig <- BiasSq_vec
  #BiasSq_vec[4] <- NA
}
##
col_slct   <- !is.na(BiasSq_vec)
BiasSq_vec <- BiasSq_vec[col_slct]
Var_vec    <- Var_vec[col_slct]
MSE_vec    <- BiasSq_vec + Var_vec 
##
##
MSE_norm_vec     <- c(MSE_vec)    / max(MSE_vec)
BiasSq_norm_vec  <- c(BiasSq_vec) / max(MSE_vec)
Var_norm_vec     <- c(Var_vec)    / max(MSE_vec)

par(mfrow=c(1,3))
barplot(MSE_norm_vec, main="",    names.arg = c("NoPS", "YesPS", "CES", "PACE", "Kraus")[col_slct], ylim = c(0,1))
mtext(text = paste0("MSRE (", round(max(MSE_vec),2),")"), side = 3, line = 1)
barplot(BiasSq_norm_vec, main="", names.arg = c("NoPS", "YesPS", "CES", "PACE", "Kraus")[col_slct], ylim = c(0,1))
mtext(text = "Squared Bias", side = 3, line = 1)
barplot(Var_norm_vec, main="",    names.arg = c("NoPS", "YesPS", "CES", "PACE", "Kraus")[col_slct], ylim = c(0,1))
mtext(text = "Variance", side = 3, line = 1)
par(mfrow=c(1,1))      


# Interpretation: 
# 1. High variance in Kraus is due to (i) connection-points and (ii) non-smooth covariance estimate
# 2. CE-Scores (PACE-procedure) are extremely bad when estimated without a noise component as the noise component 
# has a ridge-type regularization-effect. 



