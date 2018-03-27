## R-packages 
library("devtools")
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("doParallel", "fdapace")
library("ReconstPoFD")  # contains the function 'reconstuct()'
library("doParallel")   # parallel-looping


## #######################################
## Set seed
set.seed(873)
## Number of MC-Repetitions
B         <-  10
## #######################################

## #######################################
## Lower-value of the total domain 
a         <-   0             
## Upper-value of the total domain
b         <-   1
## Regular grid points
nRegGrid  <-  51 
##
n_target_fcts <- 5
##
determ_obs_interv <- c((a+(b-a)*0.33), (b-(b-a)*0.33))
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


for(DGP in c('DGP1','DGP2','DGP3','DGP4')){
  for(n in c(100, 200)){
    if(any(DGP==c('DGP1','DGP2'))){m_seq <- c(15,25,50)}else{m_seq <- NA}
    for(m in m_seq){
      
      ## a <- 0; b <- 1; DGP <- 'DGP3'; n <- 70; m <- 25; nRegGrid <- 51; B <- 10
      
      ## #######################################################################
      cat(DGP,"n=",n,"m=",m,"\n")
      ## #######################################################################
      ## #######################################################################
      ## Preselecting partially observed target functions for the reconstructions
      ## #######################################################################
      if(DGP=='DGP1'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, n_basis = 10, DGP='DGP1', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP2'){SimDat <- ReconstPoFD::simuldata(n = n_target_fcts, m = m, a = a, b = b, n_basis = 10, DGP='DGP2', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP3'){SimDat <- ReconstPoFD::simuldataKraus(n = n_target_fcts, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
      if(DGP=='DGP4'){SimDat <- ReconstPoFD::simuldataKraus(n = n_target_fcts, a = a, b = b, DGP='DGP4', nRegGrid = nRegGrid, determ_obs_interv = determ_obs_interv)}
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
      K_PACE_MC_vec     <- rep(NA, B)
      df_Kraus_MC_vec   <- rep(NA, B)
      ##
      ## #######################################################################
      ## Start 'foreach'-loop use "%do%" for a non-parallel loop
      ## #######################################################################
      ## sim.results <- foreach(repet=1:B, .combine=cbind)  %dopar% { 
      for(repet in 1:B){
        ##
        if(DGP=='DGP1'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, n_basis = 10, DGP='DGP1', nRegGrid = nRegGrid)}
        if(DGP=='DGP2'){SimDat <- ReconstPoFD::simuldata(n = n-n_target_fcts, m = m, a = a, b = b, n_basis = 10, DGP='DGP2', nRegGrid = nRegGrid)}
        if(DGP=='DGP3'){SimDat <- ReconstPoFD::simuldataKraus(n = n-n_target_fcts, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid)}
        if(DGP=='DGP4'){SimDat <- ReconstPoFD::simuldataKraus(n = n-n_target_fcts, a = a, b = b, DGP='DGP4', nRegGrid = nRegGrid)}
        ##
        Y_mat      <- cbind(Y_target_mat, SimDat[['Y_mat']])
        U_mat      <- cbind(U_target_mat, SimDat[['U_mat']])
        Y_list     <- c(Y_target_list,    SimDat[['Y_list']])
        U_list     <- c(U_target_list,    SimDat[['U_list']])
        ##
        ## Reconstruction Operator 'without Pre-Smoothing'
        result_PS_FALSE <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                    Lu           = U_list,
                                                    K            = NULL,
                                                    K_max        = 5,
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
                                                   K_max        = 5,
                                                   method       = "PS_TRUE",
                                                   reconst_fcts = target_fcts,
                                                   nRegGrid     = nRegGrid)
        Y_PS_TRUE_mat <- matrix(unlist(result_PS_TRUE[['Y_reconst_list']]), nrow = nRegGrid, ncol = n_target_fcts) 
        Y_PS_TRUE_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_PS_TRUE_mat
        K_PS_TRUE_MC_vec[repet]  <- result_PS_TRUE[['K']]
        ## 
        ## Reconstruction Operator 'without Pre-Smoothing'
        result_CEScores <- ReconstPoFD::reconstruct(Ly           = Y_list, 
                                                    Lu           = U_list,
                                                    K            = NULL,
                                                    K_max        = 5,
                                                    method       = "CEScores",
                                                    reconst_fcts = target_fcts, 
                                                    nRegGrid     = nRegGrid)
        Y_CEScores_mat <- matrix(unlist(result_CEScores[['Y_reconst_list']]), nrow = nRegGrid, ncol = n_target_fcts) 
        Y_CEScores_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_CEScores_mat
        ## K is as in PACE
        ##
        ## PACE of Yao, Müller, Wang (2005, JASA)
        result_PACE <- fdapace::FPCA(Ly    = Y_list, 
                                     Lt    = U_list, 
                                     optns = list(
                                       "dataType"       = "Sparse", 
                                       "kernel"         = "gauss",
                                       "methodMuCovEst" = "smooth",
                                       "error"          = TRUE,#ifelse(any(DGP==c('DGP1','DGP2')),TRUE,FALSE),
                                       "nRegGrid"       = nRegGrid
                                     ))
        Y_PACE_mat <- t(fitted(result_PACE))[,target_fcts]
        Y_PACE_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_PACE_mat
        K_PACE_MC_vec[repet]  <- length(result_PACE$lambda)
        ##
        if(any(DGP==c('DGP3','DGP4'))){
          ## Reconstruction Operator of Kraus (2015, JRSSB)
          result_Kraus            <- ReconstPoFD::reconstructKraus(X_mat = Y_mat, reconst_fcts = target_fcts)
          Y_Kraus_mat             <- result_Kraus[['X_reconst_mat']]
          Y_Kraus_MC_mat[,((repet-1)*n_target_fcts+1):(repet*n_target_fcts)] <- Y_Kraus_mat
          df_Kraus_MC_vec[repet]  <- result_Kraus[['df_median']]
        }
        ##
        ## ##################################################################
        if(repet %% 10) cat("repet/B=",repet,"/",B,"\n")
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
        ## i <- 2
        slct_MC_fcts                <- seq(from = i, to = n_target_fcts*B, by=n_target_fcts)
        slct_M                      <- missings_target_mat[,i]
        ##
        par(mfrow=c(2,3))
        plot(Y_PS_FALSE_MC_mat[,i], type="b", ylim=range(Y_PS_FALSE_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="PS_FALSE")
        lines(Y_target_true_mat[,i]); points(y=Y_PS_FALSE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        plot(Y_PS_TRUE_MC_mat[,i], type="b", ylim=range(Y_PS_TRUE_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="PS_TRUE")
        lines(Y_target_true_mat[,i]); points(y=Y_PS_TRUE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        plot(Y_CEScores_MC_mat[,i], type="b", ylim=range(Y_CEScores_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="CEScores")
        lines(Y_target_true_mat[,i]); points(y=Y_CEScores_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        plot(Y_PACE_MC_mat[,i], type="b", ylim=range(Y_PACE_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="PACE")
        lines(Y_target_true_mat[,i]); points(y=Y_PACE_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        ##
        plot(Y_Kraus_MC_mat[,i], type="b", ylim=range(Y_Kraus_MC_mat[,i],Y_target_true_mat[slct_M,i]),main="Kraus")
        lines(Y_target_true_mat[,i]); points(y=Y_Kraus_MC_mat[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
        par(mfrow=c(1,1))
        ##
        PS_FALSE_Int_BiasSq_vec[i] <- sum(c(rowMeans(Y_PS_FALSE_MC_mat[slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        PS_TRUE_Int_BiasSq_vec[i]  <- sum(c(rowMeans(Y_PS_TRUE_MC_mat[ slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        CEScores_Int_BiasSq_vec[i] <- sum(c(rowMeans(Y_CEScores_MC_mat[slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        PACE_Int_BiasSq_vec[i]     <- sum(c(rowMeans(Y_PACE_MC_mat[    slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        Kraus_Int_BiasSq_vec[i]    <- sum(c(rowMeans(Y_Kraus_MC_mat[   slct_M, slct_MC_fcts]) - Y_target_true_mat[slct_M, i])^2) * (b-a)/nRegGrid
        ##
        PS_FALSE_Int_Var_vec[i]     <- sum(apply(Y_PS_FALSE_MC_mat[slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        PS_TRUE_Int_Var_vec[i]      <- sum(apply(Y_PS_TRUE_MC_mat[ slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        CEScores_Int_Var_vec[i]     <- sum(apply(Y_CEScores_MC_mat[slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        PACE_Int_Var_vec[i]         <- sum(apply(Y_PACE_MC_mat[    slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
        Kraus_Int_Var_vec[i]        <- sum(apply(Y_Kraus_MC_mat[   slct_M, slct_MC_fcts], 1, var)) * (b-a)/nRegGrid
      }
      ##
      # ## Save results:
      # if(any(DGP==c('DGP1','DGP2'))){
      #   save(sim.results, file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
      # }
      # if(any(DGP==c('DGP3','DGP4'))){
      #   save(sim.results, file = paste0(DGP,"_n",n,"_simResults.RData"))
      # }
    }
  }
}
##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------