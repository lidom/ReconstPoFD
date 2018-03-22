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
nRegGrid  <-  75  
## #######################################

## #######################################
## Number of cores used for the outer loop
## Parallel loop might not work for windows!
## => use %do% instead of %dopar% (below) 
## for a non-parallel loop"
registerDoParallel(cores=4)
getDoParWorkers()
## #######################################


(Start.Time <- Sys.time())


for(DGP in c('DGP1','DGP2','DGP3','DGP4')){
  for(n in c(100, 200)){
    if(any(DGP==c('DGP1','DGP2'))){m_seq <- c(15,25,50)}else{m_seq <- NA}
    for(m in m_seq){
      ## #######################################################################
      cat(DGP,"n=",n,"m=",m,"\n")
      ## #######################################################################
      
      ## #######################################################################
      ## Start 'foreach'-loop use "%do%" for a non-parallel loop
      ## #######################################################################
      sim.results <- foreach(repet=1:B, .combine=cbind)  %dopar% { 
        ##
        if(DGP=='DGP1'){
          SimDat <- ReconstPoFD::simuldata(n = n, m = m, a = a, b = b, n_basis = 10, DGP='DGP1', nRegGrid = nRegGrid)
        }
        if(DGP=='DGP2'){
          SimDat <- ReconstPoFD::simuldata(n = n, m = m, a = a, b = b, n_basis = 10, DGP='DGP2', nRegGrid = nRegGrid)
        }
        if(DGP=='DGP3'){
          SimDat <- ReconstPoFD::simuldataKraus(n = n, a = a, b = b, DGP='DGP3', nRegGrid = nRegGrid)
        }
        if(DGP=='DGP4'){
          SimDat <- ReconstPoFD::simuldataKraus(n = n, a = a, b = b, DGP='DGP4', nRegGrid = nRegGrid)
        }
        ##
        Y_true_mat <- SimDat[['Y_true_mat']]
        U_true_mat <- SimDat[['U_true_mat']]
        ##
        Y_mat      <- SimDat[['Y_mat']]
        U_mat      <- SimDat[['U_mat']]
        Y_list     <- SimDat[['Y_list']]
        U_list     <- SimDat[['U_list']]
        ##
        A_vec      <- SimDat[['A_vec']]
        B_vec      <- SimDat[['B_vec']]
        ##
        ## Reconstruction Operator 'without Pre-Smoothing'
        result_PS_FALSE <- ReconstPoFD::reconstruct(Ly         = Y_list, 
                                                    Lu         = U_list,
                                                    K          = NULL,
                                                    K_max      = 10,
                                                    pre_smooth = FALSE,
                                                    nRegGrid   = nRegGrid)
        Y_PS_FALSE_mat <- matrix(unlist(result_PS_FALSE[['Y_reconst_list']]), nrow = nRegGrid, ncol = n) 
        U_PS_FALSE_mat <- matrix(unlist(result_PS_FALSE[['U_reconst_list']]), nrow = nRegGrid, ncol = n) 
        K_PS_FALSE     <- result_PS_FALSE[['K_AIC']]
        ##
        ## Reconstruction Operator 'with Pre-Smoothing'
        result_PS_TRUE <- ReconstPoFD::reconstruct(Ly         = Y_list, 
                                                   Lu         = U_list,
                                                   K          = NULL,
                                                   K_max      = 10,
                                                   pre_smooth = TRUE,
                                                   nRegGrid   = nRegGrid)
        Y_PS_TRUE_mat <- matrix(unlist(result_PS_TRUE[['Y_reconst_list']]), nrow = nRegGrid, ncol = n) 
        U_PS_TRUE_mat <- matrix(unlist(result_PS_TRUE[['U_reconst_list']]), nrow = nRegGrid, ncol = n) 
        K_PS_TRUE     <- result_PS_TRUE[['K_AIC']]
        ## 
        ## PACE of Yao, Müller, Wang (2005, JASA)
        result_PACE <- fdapace::FPCA(Ly    = Y_list, 
                                     Lt    = U_list, 
                                     optns = list(
                                       "dataType"       = "Sparse", 
                                       "kernel"         = "gauss",
                                       "methodMuCovEst" = "smooth",
                                       "nRegGrid"       = nRegGrid
                                     ))
        Y_PACE_mat <- t(fitted(result_PACE))
        U_PACE_mat <- matrix(result_PACE$workGrid, nrow = nRegGrid, ncol = n)
        K_PACE     <- length(result_PACE$lambda)
        ##
        if(any(DGP==c('DGP3','DGP4'))){
          ## Reconstruction Operator of Kraus (2015, JRSSB)
          result_Kraus         <- ReconstPoFD::reconstructKraus(X_mat = Y_mat)
          Y_Kraus_mat          <- result_Kraus[['X_reconst_mat']]
          df_gcv_median_Kraus  <- result_Kraus[['df_gcv_median']]
        }else{df_gcv_median_Kraus <- NA}
        
        
        ## Reconstruction errors
        L2_PS_FALSE_vec  <- rep(NA, n)
        L2_PS_TRUE_vec   <- rep(NA, n)
        L2_PACE_vec      <- rep(NA, n)
        L2_Kraus_vec     <- rep(NA, n)
        for(i in 1:n){
          ##
          missing_vec <- rep(FALSE, nRegGrid)
          missing_vec[U_true_mat[,i] < A_vec[i]] <- TRUE
          missing_vec[U_true_mat[,i] > B_vec[i]] <- TRUE
          ## Integrated squared errors
          L2_PS_FALSE_vec[i] <- sum(c(Y_PS_FALSE_mat[missing_vec,i] - Y_true_mat[missing_vec,i])^2) * (b-a)/nRegGrid
          L2_PS_TRUE_vec[i]  <- sum(c(Y_PS_TRUE_mat[missing_vec,i]  - Y_true_mat[missing_vec,i])^2) * (b-a)/nRegGrid
          L2_PACE_vec[i]     <- sum(c(Y_PACE_mat[missing_vec,i]     - Y_true_mat[missing_vec,i])^2) * (b-a)/nRegGrid
          if(any(DGP==c('DGP3','DGP4'))){
            L2_Kraus_vec[i]  <- sum(c(Y_Kraus_mat[missing_vec,i]    - Y_true_mat[missing_vec,i])^2) * (b-a)/nRegGrid
          }
        }
        ## ##################################################################
        if(B %% 10) cat("repet/B=",repet,"/",B,"\n")
        ## ##################################################################
        ## Return the results (Mean integrated squared errors and K)
        c(mean(L2_PS_FALSE_vec), mean(L2_PS_TRUE_vec), mean(L2_PACE_vec), mean(L2_Kraus_vec), K_PS_FALSE, K_PS_TRUE, K_PACE, df_gcv_median_Kraus)
        ## ##################################################################
      } ## End of foreach-loop
      ## Row-names:
      rownames(sim.results) <- c("ML2_PS_FALSE", "ML2_PS_TRUE", "ML2_PACE", "ML2_Kraus", "K_PS_FALSE", "K_PS_TRUE", "K_PACE", "df_gcv_median_Kraus")
      
      ## Save results:
      if(any(DGP==c('DGP1','DGP2'))){
        save(sim.results, file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
      }
      if(any(DGP==c('DGP3','DGP4'))){
        save(sim.results, file = paste0(DGP,"_n",n,"_simResults.RData"))
      }
    }
  }
}
##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------