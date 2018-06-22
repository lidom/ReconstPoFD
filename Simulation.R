## R-packages 
library("devtools")
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("doParallel", "fdapace")
library("ReconstPoFD")  # contains the function 'reconstuct()'
#library("doParallel")   # parallel-looping

## Location to store the simulation results:
setwd("/home/dom/ownCloud/Kneip_Liebl_Reconstruction/Simulation_Submission_2")

## #######################################
## Set seed
set.seed(873)
## Number of MC-Repetitions
B         <-  500
## #######################################

## #######################################
## Lower-value of the total domain 
a         <-   0             
## Upper-value of the total domain
b         <-   1
## Regular grid points
nRegGrid  <-  61 
##
n_target_fcts <- 1
##
determ_obs_interv <- c((a+(b-a)*0.35), (b-(b-a)*0.35))
## #######################################

## #######################################
## Number of cores used for the outer loop
## Parallel loop might not work for windows!
## => use %do% instead of %dopar% (below) 
## for a non-parallel loop"
#registerDoParallel(cores=4)
#getDoParWorkers()
## #######################################


(Start.Time <- Sys.time())


for(DGP in c('DGP1','DGP4','DGP5')[1]){
  for(n in c(50, 100)){
    if(any(DGP==c('DGP1','DGP2'))){m_seq <- c(15, 50)}else{m_seq <- NA}
    for(m in m_seq){
      
      ## a <- 0; b <- 1; DGP <- 'DGP5'; n <- 50; m <- 50; nRegGrid <- 61; B <- 30
      
      ## #######################################################################
      cat(DGP,"n=",n,"m=",m,"\n")
      ## #######################################################################
      ## #######################################################################
      ## Preselecting partially observed target functions for the reconstructions
      ## #######################################################################
      if(DGP=='DGP1'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, DGP='DGP1', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP2'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, DGP='DGP2', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP3'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP4'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, DGP='DGP4', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP5'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, DGP='DGP5', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      ## 
      target_fcts       <- 1:n_target_fcts
      ##
      Y_target_true_mat <- SimDat[['Y_true_mat']]
      U_target_true_mat <- SimDat[['U_true_mat']]
      ##
      Y_target_mat      <- SimDat[['Y_mat']]
      U_target_mat      <- SimDat[['U_mat']]
      Y_target_list     <- SimDat[['Y_list']]
      U_target_list     <- SimDat[['U_list']]
      ##
      A_target_vec      <- SimDat[['A_vec']]
      B_target_vec      <- SimDat[['B_vec']]
      ##
      missings_target_mat <- matrix(NA, nrow = nRegGrid, ncol = n_target_fcts)
      for(i in 1:n_target_fcts){
        missings_target_mat[,i] <- c(U_target_true_mat[,i] < A_target_vec[i] | U_target_true_mat[,i] > B_target_vec[i])
      }       
      ##
      Y_PS_FALSE_MC_mat <- matrix(NA, nrow = nRegGrid, ncol = n_target_fcts*B)
      Y_PS_TRUE_MC_mat  <- matrix(NA, nrow = nRegGrid, ncol = n_target_fcts*B)
      Y_CEScores_MC_mat <- matrix(NA, nrow = nRegGrid, ncol = n_target_fcts*B)
      Y_PACE_MC_mat     <- matrix(NA, nrow = nRegGrid, ncol = n_target_fcts*B)
      Y_Kraus_MC_mat    <- matrix(NA, nrow = nRegGrid, ncol = n_target_fcts*B)
      K_PS_FALSE_MC_vec <- rep(NA, B)
      K_PS_TRUE_MC_vec  <- rep(NA, B)
      K_CEScores_MC_vec <- rep(NA, B)
      K_PACE_MC_vec     <- rep(NA, B)
      df_Kraus_MC_vec   <- rep(NA, B)
      ##
      ## #######################################################################
      ## Start 'foreach'-loop use "%do%" for a non-parallel loop
      ## #######################################################################
      ## sim.results <- foreach(repet=1:B, .combine=cbind)  %dopar% { 
      for(repet in 1:B){
        ##
        if(DGP=='DGP1'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, DGP='DGP1', nRegGrid = nRegGrid)}
        if(DGP=='DGP2'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, DGP='DGP2', nRegGrid = nRegGrid)}
        if(DGP=='DGP3'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid)}
        if(DGP=='DGP4'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, DGP='DGP4', nRegGrid = nRegGrid)}
        if(DGP=='DGP5'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, DGP='DGP5', nRegGrid = nRegGrid)}
        ##
        Y_mat      <- cbind(Y_target_mat, SimDat[['Y_mat']])
        U_mat      <- cbind(U_target_mat, SimDat[['U_mat']])
        ##
        Y_list     <- c(Y_target_list,    SimDat[['Y_list']])
        U_list     <- c(U_target_list,    SimDat[['U_list']])
        ##
        if(!is.na(m_seq)){
        if(n==50 & m<75){
          h.mu  <- (max(c(unlist(U_list))) - min(c(unlist(U_list)))) * 0.05
          h.cov <- (max(c(unlist(U_list))) - min(c(unlist(U_list)))) * 0.10
        }
        if(n==100 | m==75){
          h.mu  <- (max(c(unlist(U_list))) - min(c(unlist(U_list)))) * 0.05/2
          h.cov <- (max(c(unlist(U_list))) - min(c(unlist(U_list)))) * 0.10/2
        }}
        ##
        ## Reconstruction Operator 'without Pre-Smoothing'
        result_PS_FALSE <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                    Lu           = U_list,
                                                    K            = NULL,
                                                    K_max        = 3,#4
                                                    method       = "PS_FALSE",
                                                    reconst_fcts = target_fcts, 
                                                    nRegGrid     = nRegGrid)
        Y_PS_FALSE_mat <- matrix(unlist(result_PS_FALSE[['Y_reconst_list']]), nrow = nRegGrid, ncol = n_target_fcts) 
        Y_PS_FALSE_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_PS_FALSE_mat
        K_PS_FALSE_MC_vec[repet]  <- result_PS_FALSE[['K']]
        ##
        ## Reconstruction Operator 'with Pre-Smoothing'
        result_PS_TRUE <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                   Lu           = U_list,
                                                   K            = NULL,
                                                   K_max        = 3,#4
                                                   method       = "PS_TRUE",
                                                   reconst_fcts = target_fcts,
                                                   nRegGrid     = nRegGrid)
        Y_PS_TRUE_mat <- matrix(unlist(result_PS_TRUE[['Y_reconst_list']]), nrow = nRegGrid, ncol = n_target_fcts) 
        Y_PS_TRUE_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_PS_TRUE_mat
        K_PS_TRUE_MC_vec[repet]  <- result_PS_TRUE[['K']]
        ## 
        ## Reconstruction Operator 'without Pre-Smoothing'
        if(any(DGP==c('DGP1','DGP2'))){
          result_CEScores <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                      Lu           = U_list,
                                                      K            = NULL,
                                                      K_max        = 3,#4
                                                      method       = "CEScores",
                                                      BwMu         = h.mu,
                                                      BwCov        = h.cov,
                                                      reconst_fcts = target_fcts, 
                                                      nRegGrid     = nRegGrid)
          Y_CEScores_mat <- matrix(unlist(result_CEScores[['Y_reconst_list']]), nrow = nRegGrid, ncol = n_target_fcts) 
          Y_CEScores_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_CEScores_mat
          K_CEScores_MC_vec[repet]  <- result_CEScores[['K']]
        }
          ##
          ## PACE of Yao, Müller, Wang (2005, JASA) if error == true
          result_PACE <- fdapace::FPCA(Ly    = Y_list, 
                                       Lt    = U_list, 
                                       optns = list(
                                         "dataType"       = "Sparse", 
                                         "kernel"         = "gauss",
                                         "methodMuCovEst" = "smooth",
                                         "error"          = ifelse(any(DGP==c('DGP1','DGP2')),TRUE,FALSE),
                                         "userBwMu"       = h.mu,
                                         "userBwCov"      = h.cov,
                                         "nRegGrid"       = nRegGrid
                                       ))
          Y_PACE_mat <- t(fitted(result_PACE))[,target_fcts]
          Y_PACE_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_PACE_mat
          K_PACE_MC_vec[repet]  <- length(result_PACE$lambda)
        #}
        ##
        if(any(DGP==c('DGP3','DGP4','DGP5'))){
          ## Reconstruction Operator of Kraus (2015, JRSSB)
          result_Kraus            <- ReconstPoFD::reconstructKraus(X_mat = Y_mat, reconst_fcts = target_fcts)
          Y_Kraus_mat             <- result_Kraus[['X_reconst_mat']]
          Y_Kraus_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_Kraus_mat
          df_Kraus_MC_vec[repet]  <- result_Kraus[['df_median']]
        }
        ##
        ## ##################################################################
        if(repet %% 100 == 0) cat("repet/B=",repet,"/",B,"\n")
        ## ##################################################################
      } ## End of B-loop
      ##
      ## Integrated squared bias:
      PS_FALSE_Int_BiasSq_vec <- rep(NA, n_target_fcts)
      PS_TRUE_Int_BiasSq_vec  <- rep(NA, n_target_fcts)
      CEScores_Int_BiasSq_vec <- rep(NA, n_target_fcts)
      PACE_Int_BiasSq_vec     <- rep(NA, n_target_fcts)
      Kraus_Int_BiasSq_vec    <- rep(NA, n_target_fcts)
      ## Integrated variance
      PS_FALSE_Int_Var_vec     <- rep(NA, n_target_fcts)
      PS_TRUE_Int_Var_vec      <- rep(NA, n_target_fcts)
      CEScores_Int_Var_vec     <- rep(NA, n_target_fcts)
      PACE_Int_Var_vec         <- rep(NA, n_target_fcts)
      Kraus_Int_Var_vec        <- rep(NA, n_target_fcts)
      ##
      for(i in 1:n_target_fcts){
        ## i <- 1
        slct_MC_fcts                <- seq(from = i, to = n_target_fcts*B, by=n_target_fcts)
        slct_M                      <- missings_target_mat[,i]
        ##
        par(mfrow=c(2,2))
        ##
        plot(Y_PS_FALSE_MC_mat[,i], type="b", ylim=range(Y_PS_FALSE_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="PS_FALSE")
        lines(Y_target_true_mat[,i]); points(y=Y_PS_FALSE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        plot(Y_PS_TRUE_MC_mat[,i], type="b", ylim=range(Y_PS_TRUE_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="PS_TRUE")
        lines(Y_target_true_mat[,i]); points(y=Y_PS_TRUE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        if(any(DGP==c('DGP1','DGP2'))){
        plot(Y_CEScores_MC_mat[,i], type="b", ylim=range(Y_CEScores_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="CEScores")
        lines(Y_target_true_mat[,i]); points(y=Y_CEScores_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        }
        ##
        plot(Y_PACE_MC_mat[,i], type="b", ylim=range(Y_PACE_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="PACE")
        lines(Y_target_true_mat[,i]); points(y=Y_PACE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        if(any(DGP==c('DGP3','DGP4','DGP5'))){
        plot(Y_Kraus_MC_mat[,i], type="b", ylim=range(Y_Kraus_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="Kraus")
        lines(Y_target_true_mat[,i]); points(y=Y_Kraus_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        }
        par(mfrow=c(1,1))
        ##
        ## \int_{Missing} (Bias(t))^2 dt ('Integrated squared bias')
        PS_FALSE_Int_BiasSq_vec[i] <- sum(c(rowMeans(Y_PS_FALSE_MC_mat[slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        PS_TRUE_Int_BiasSq_vec[i]  <- sum(c(rowMeans(Y_PS_TRUE_MC_mat[ slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        CEScores_Int_BiasSq_vec[i] <- sum(c(rowMeans(Y_CEScores_MC_mat[slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        PACE_Int_BiasSq_vec[i]     <- sum(c(rowMeans(Y_PACE_MC_mat[    slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        Kraus_Int_BiasSq_vec[i]    <- sum(c(rowMeans(Y_Kraus_MC_mat[   slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        ## \int_{Missing} Var(t) dt ('Integrated variance')
        PS_FALSE_Int_Var_vec[i]     <- sum(apply(Y_PS_FALSE_MC_mat[slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        PS_TRUE_Int_Var_vec[i]      <- sum(apply(Y_PS_TRUE_MC_mat[ slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        CEScores_Int_Var_vec[i]     <- sum(apply(Y_CEScores_MC_mat[slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        PACE_Int_Var_vec[i]         <- sum(apply(Y_PACE_MC_mat[    slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        Kraus_Int_Var_vec[i]        <- sum(apply(Y_Kraus_MC_mat[   slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
      }
      ##
      BiasSq <- c(mean(PS_FALSE_Int_BiasSq_vec),mean(PS_TRUE_Int_BiasSq_vec),mean(CEScores_Int_BiasSq_vec),mean(PACE_Int_BiasSq_vec),mean(Kraus_Int_BiasSq_vec))
      Var    <- c(mean(PS_FALSE_Int_Var_vec),   mean(PS_TRUE_Int_Var_vec),   mean(CEScores_Int_Var_vec),   mean(PACE_Int_Var_vec),   mean(Kraus_Int_Var_vec))
      ##
      par(mfrow=c(1,3))
      barplot(c(BiasSq+Var)[!is.na(BiasSq)]/max(c(BiasSq+Var)[!is.na(BiasSq)]), main="", names.arg = c("NoPS","YesPS","CES","PACE","Kraus")[!is.na(BiasSq)])
      mtext(text = paste0("MSPE (", round(max(c(BiasSq+Var)[!is.na(BiasSq)]),2),")"), side = 3, line = 1)
      barplot(c(BiasSq    )[!is.na(BiasSq)]/max(c(BiasSq+Var)[!is.na(BiasSq)]), main="", names.arg = c("NoPS","YesPS","CES","PACE","Kraus")[!is.na(BiasSq)])
      mtext(text = "Squared Bias", side = 3, line = 1)
      barplot(c(Var       )[!is.na(BiasSq)]/max(c(BiasSq+Var)[!is.na(BiasSq)]), main="", names.arg = c("NoPS","YesPS","CES","PACE","Kraus")[!is.na(BiasSq)])
      mtext(text = "Variance", side = 3, line = 1)
      par(mfrow=c(1,1))      
      ## Save results:
      if(any(DGP==c('DGP1','DGP2'))){
        save(BiasSq, Var, file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
      }
      if(any(DGP==c('DGP3','DGP4','DGP5'))){
        save(BiasSq, Var, file = paste0(DGP,"_n",n,"_simResults.RData"))
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

DGP <- c('DGP1','DGP4','DGP5')[3]
m   <- c(15,  50, 75)[1] # c(15,25,50)[3]
n   <- c(50, 100)[2]     # c(100, 200)[1]

## Load results:
if(any(DGP==c('DGP1','DGP2'))){
  load(file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
  ## Relative Bias^2 and Var (basis: MSE)
  bias2_vec  <- c(BiasSq)[!is.na(BiasSq)] / max(c(BiasSq+Var)[!is.na(BiasSq)])
  var_vec    <- c(Var   )[!is.na(BiasSq)] / max(c(BiasSq+Var)[!is.na(BiasSq)])
}
if(any(DGP==c('DGP3','DGP4','DGP5'))){
  load(file = paste0(DGP,"_n",n,"_simResults.RData"))
  ## Relative Bias^2 and Var (basis: MSE)
  BiasSq[4]  <- NA
  bias2_vec  <- c(BiasSq)[!is.na(BiasSq)] / max(c(BiasSq+Var)[!is.na(BiasSq)])
  var_vec    <- c(Var   )[!is.na(BiasSq)] / max(c(BiasSq+Var)[!is.na(BiasSq)])
}
##
winsorize_x <- function(x, cut = 0.01){
  cut_point_top    <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x,     cut, na.rm = T)
  i <-  which(x >= cut_point_top);    x[i] <- cut_point_top
  j <-  which(x <= cut_point_bottom); x[j] <- cut_point_bottom
  return(x)
}
##
cut <- 0.0#5
##
bias2_cut_vec <- winsorize_x(bias2_vec, cut=cut)
var_cut_vec   <- winsorize_x(var_vec,   cut=cut)
mse_cut_vec   <- bias2_vec + var_vec
##
par(mfrow=c(1,3))
barplot(mse_cut_vec, main="",   names.arg = c("NoPS", "YesPS", "CES", "PACE", "Kraus")[!is.na(BiasSq)], ylim = c(0,1))
mtext(text = paste0("MSRE (", round(max(c(BiasSq+Var)[!is.na(BiasSq)]),2),")"), side = 3, line = 1)
barplot(bias2_cut_vec, main="", names.arg = c("NoPS", "YesPS", "CES", "PACE", "Kraus")[!is.na(BiasSq)], ylim = c(0,1))
mtext(text = "Squared Bias", side = 3, line = 1)
barplot(var_cut_vec, main="",   names.arg = c("NoPS", "YesPS", "CES", "PACE", "Kraus")[!is.na(BiasSq)], ylim = c(0,1))
mtext(text = "Variance", side = 3, line = 1)
par(mfrow=c(1,1))      


