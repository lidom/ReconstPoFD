## R-packages 
library("devtools")
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("doParallel", "fdapace")
library("ReconstPoFD")  # contains the function 'reconstuct()'
library("doParallel")   # parallel-looping


## #######################################
## Set seed
set.seed(873)
## DGPs
DGP       <-  c('DGP1','DGP2','DGP3','DGP4')
## Number of discretization points for DGP1 and DGP2 
m         <-  c(15,25,50,75)[1]  
## Number of functions 
n         <-  c(100,200)[1]
## Number of MC-Repetitions
B         <-  2   
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
registerDoParallel(cores=2)
getDoParWorkers()
## #######################################

## #######################################################################
## Start 'foreach'-loop use "%do%" for a non-parallel loop
## #######################################################################
sim.results <- foreach(repet=1:B, .combine=cbind)  %dopar% { 
  ##
  if(DGP=='DGP1'){
    SimDat <- ReconstPoFD::simuldata(n = n, m = m, a = a, b = b, n_basis = 10, DGP='DGP1')
  }
  if(DGP=='DGP2'){
    SimDat <- ReconstPoFD::simuldata(n = n, m = m, a = a, b = b, n_basis = 10, DGP='DGP2')
  }
  if(DGP=='DGP3'){
    SimDat <- ReconstPoFD::simuldataKraus(n = n, a = a, b = b, DGP='DGP3')
  }
  if(DGP=='DGP4'){
    SimDat <- ReconstPoFD::simuldataKraus(n = n, a = a, b = b, DGP='DGP4')
  }
  Y_list   <- SimDat[['Y_list']]
  Y_mat    <- SimDat[['Y_mat']]
  U_list   <- SimDat[['U_list']]
  U_mat    <- SimDat[['U_mat']]
  ##
  result_PS_FALSE <- ReconstPoFD::reconstruct(Ly         = Y_list, 
                                              Lu         = U_list,
                                              K          = NULL,
                                              K_max      = 15,
                                              pre_smooth = FALSE,
                                              nRegGrid   = nRegGrid)
  Y_PS_FALSE_mat <- matrix(unlist(result_PS_FALSE[['Y_reconst_list']]), nrow = nRegGrid, ncol = n) 
  U_PS_FALSE_mat <- matrix(unlist(result_PS_FALSE[['U_reconst_list']]), nrow = nRegGrid, ncol = n) 
  ##
  if(any(DGP==c('DGP1','DGP2'))){
    result_PS_TRUE <- ReconstPoFD::reconstruct(Ly         = Y_list, 
                                               Lu         = U_list,
                                               K          = NULL,
                                               K_max      = 15,
                                               pre_smooth = TRUE,
                                               nRegGrid   = nRegGrid)
  }else{result_PS_TRUE <- NULL}
  ##
  if(any(DGP==c('DGP3','DGP4'))){
    reconst_Kraus <- ReconstPoFD::reconstructKraus(X_mat = Y_mat)
  }else{result_PS_TRUE <- NULL}
  
  ## K chosen by AIC-type criterion:
  K.norm <- reconst_norm_obj[['K_AIC']]
  K.expd <- reconst_expd_obj[['K_AIC']]
  
  ## Reconstruction errors
  MaxAE_norm_vec     <- rep(NA, n)
  MaxAE_expd_vec     <- rep(NA, n)
  for(i in 1:n){
    y_norm_tmp <- reconst_norm_obj[['y_reconst_list']][[i]]
    x_norm_tmp <- reconst_norm_obj[['x_reconst_list']][[i]]
    ##
    y_expd_tmp <- reconst_expd_obj[['y_reconst_list']][[i]]
    x_expd_tmp <- reconst_expd_obj[['x_reconst_list']][[i]]
    ## Max absolute error:
    MaxAE_norm_vec[i] <- max(abs(y_norm_tmp - True_norm_curve_fun(u=x_norm_tmp, i=i)))
    MaxAE_expd_vec[i] <- max(abs(y_expd_tmp - True_expd_curve_fun(u=x_expd_tmp, i=i)))
  }
  
  ## ##################################################################
  cat("repet/B=",repet,"/",B,"\n")
  ## ##################################################################
  ## Return the results:
  c(mean(MaxAE_norm_vec), mean(MaxAE_expd_vec), K.norm, K.expd)
  ## ##################################################################
} ## End of foreach-loop.

## Row-names:
rownames(sim.results) <- c("MAE_norm", "MAE_expd", "K.norm", "K.expd")

## sim.results

## Save results:
## save(sim.results, file = paste0("n",n,"_m",m,"_simResults.RData"))
