## R-packages 
library("devtools")
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("doParallel", "fdapace")
library("ReconstPoFD")  # contains the function 'reconstuct()'
library("doParallel")   # parallel-looping


## #######################################
## Set seed
set.seed(873)
## Number of discretization points
m         <-  c(15,25,50,75)[1]  
## Number of functions 
n         <-  c(100,200)[1]
## Number of MC-Repetitions
B         <-  2   
## #######################################

## #######################################
## Lower-value of the total domain 
a         <-    0             
## Upper-value of the total domain
b         <-   10
## Regular grid points
nRegGrid  <-  101  
## #######################################

## #######################################
## Number of cores used for the outer loop
## Parallel loop might not work for windows!
## => use %do% instead of %dopar% (below) 
## for a non-parallel loop"
registerDoParallel(cores=2)
getDoParWorkers()
## #######################################

## ##################################
## Variance-parameters: #############
lambda1  <- 4# This one 
lambda2  <- 0#
lambda3  <- 3# This one 
lambda4  <- 0#
lambda5  <- 0#
lambda6  <- 0#
## 
lam_vec <- c(lambda1,lambda2,lambda3,lambda4,lambda5,lambda6)
(K_true  <- length(lam_vec[lam_vec>0]))
##
eps_var  <- .20
## ##################################


## #########################################################
## mean-/eigen-functions and eigenvalues ###################
mean_fun <- function(u){return(u+sin(u))}
e1_fun   <- function(u){return(-cos(1*pi*u/(b-a))/sqrt(5))}# This one 
e2_fun   <- function(u){return( sin(1*pi*u/(b-a))/sqrt(5))}
e3_fun   <- function(u){return(-cos(2*pi*u/(b-a))/sqrt(5))}# This one 
e4_fun   <- function(u){return( sin(2*pi*u/(b-a))/sqrt(5))}
e5_fun   <- function(u){return(-cos(3*pi*u/(b-a))/sqrt(5))}
e6_fun   <- function(u){return( sin(3*pi*u/(b-a))/sqrt(5))}
## #########################################################

## check-plot:
par(mfrow=c(1,3))
plot(y=mean_fun(u=seq(a,b,len=25)), x=seq(a,b,len=25), type="l")
plot(y=e1_fun(  u=seq(a,b,len=25)), x=seq(a,b,len=25), type="l")
plot(y=e2_fun(  u=seq(a,b,len=25)), x=seq(a,b,len=25), type="l")
par(mfrow=c(1,1))
dev.off()


## #########################################################
## covariance function
gamma_fun <- function(u,v){
  result <- lambda1 * e1_fun(u) * e1_fun(v) + 
    lambda2 * e2_fun(u) * e2_fun(v) + 
    lambda3 * e3_fun(u) * e3_fun(v) +
    lambda4 * e4_fun(u) * e4_fun(v) + 
    lambda5 * e5_fun(u) * e5_fun(v) + 
    lambda6 * e6_fun(u) * e6_fun(v)
  return(result)
}
## ########################################################


## #######################################################################
## Start 'foreach'-loop use "%do%" for a non-parallel loop
## #######################################################################
sim.results <- foreach(repet=1:B, .combine=cbind)  %dopar% { 
  
  ## Additional measurement errors
  eps_norm_mat <- matrix(data = rnorm(n=n*m,  mean = 0, sd=sqrt(eps_var)), m, n)
  
  ## Random subintervals
  A_i_vec       <- runif(n=n, min = a, max = (a+ (b-a) * 0.25))
  B_i_vec       <- runif(n=n, min = (b- (b-a) * 0.25), max = b)
  
  ## Generation of prediction points U over (proper) random subdomains
  U_mat         <- matrix(NA, m, n)
  for(i in 1:n){
    ## sampling from the total grid
    U_mat[,i] <- runif(n=m, min = A_i_vec[i], max = B_i_vec[i])
    ## ordering
    U_mat[,i] <- unique(U_mat[,i][order(U_mat[,i])])
  }
  
  ## ####################################################
  ## Generate PC-scores: ################################
  ## Normal
  xi1_norm_vec <- rnorm(n=n, mean=0, sd=sqrt(lambda1))
  xi2_norm_vec <- rnorm(n=n, mean=0, sd=sqrt(lambda2))
  xi3_norm_vec <- rnorm(n=n, mean=0, sd=sqrt(lambda3))
  xi4_norm_vec <- rnorm(n=n, mean=0, sd=sqrt(lambda4))
  xi5_norm_vec <- rnorm(n=n, mean=0, sd=sqrt(lambda5))
  xi6_norm_vec <- rnorm(n=n, mean=0, sd=sqrt(lambda6))
  
  ## Exp
  xi1_expd_vec <- c(rexp(n=n, rate = 1/sqrt(lambda1)) - sqrt(lambda1))
  xi2_expd_vec <- c(rexp(n=n, rate = 1/sqrt(lambda2)) - sqrt(lambda2))
  xi3_expd_vec <- c(rexp(n=n, rate = 1/sqrt(lambda3)) - sqrt(lambda3))
  xi4_expd_vec <- c(rexp(n=n, rate = 1/sqrt(lambda4)) - sqrt(lambda4))
  xi5_expd_vec <- c(rexp(n=n, rate = 1/sqrt(lambda5)) - sqrt(lambda5))
  xi6_expd_vec <- c(rexp(n=n, rate = 1/sqrt(lambda6)) - sqrt(lambda6))
  
  ## #########################################################################
  ## Generate Y-data: ########################################################
  Y_norm_mat     <- matrix(NA, m, n)
  Y_expd_mat     <- matrix(NA, m, n)
  for(i in 1:n){
    Y_norm_mat[,i] <- mean_fun(u=U_mat[,i]) + 
      xi1_norm_vec[i] * e1_fun(u=U_mat[,i]) + 
      xi2_norm_vec[i] * e2_fun(u=U_mat[,i]) + 
      xi3_norm_vec[i] * e3_fun(u=U_mat[,i]) + 
      xi4_norm_vec[i] * e4_fun(u=U_mat[,i]) + 
      xi5_norm_vec[i] * e5_fun(u=U_mat[,i]) + 
      xi6_norm_vec[i] * e6_fun(u=U_mat[,i]) + 
      eps_norm_mat[,i]
    ##
    Y_expd_mat[,i] <- mean_fun(u=U_mat[,i]) + 
      xi1_expd_vec[i] * e1_fun(u=U_mat[,i]) + 
      xi2_expd_vec[i] * e2_fun(u=U_mat[,i]) + 
      xi3_expd_vec[i] * e3_fun(u=U_mat[,i]) + 
      xi4_expd_vec[i] * e4_fun(u=U_mat[,i]) + 
      xi5_expd_vec[i] * e5_fun(u=U_mat[,i]) + 
      xi6_expd_vec[i] * e6_fun(u=U_mat[,i]) + 
      eps_norm_mat[,i]
  }
  ## ##########################################################
  
  ## ##########################################################
  ## True curve function: #####################################
  True_norm_curve_fun <- function(u, i){
    mean_fun(u=u) + 
      xi1_norm_vec[i] * e1_fun(u=u) + 
      xi2_norm_vec[i] * e2_fun(u=u) +
      xi3_norm_vec[i] * e3_fun(u=u) + 
      xi4_norm_vec[i] * e4_fun(u=u) +
      xi5_norm_vec[i] * e5_fun(u=u) + 
      xi6_norm_vec[i] * e6_fun(u=u) 
  }
  ##
  True_expd_curve_fun <- function(u, i){
    mean_fun(u=u) + 
      xi1_expd_vec[i] * e1_fun(u=u) + 
      xi2_expd_vec[i] * e2_fun(u=u) +
      xi3_expd_vec[i] * e3_fun(u=u) + 
      xi4_expd_vec[i] * e4_fun(u=u) +
      xi5_expd_vec[i] * e5_fun(u=u) + 
      xi6_expd_vec[i] * e6_fun(u=u) 
  }
  ## #############################################################
  ## Discretizing the true mean and cov
  true_mu_vec  <- mean_fun(seq(a,b,len=nRegGrid))
  true_cov_mat <- matrix(NA, nRegGrid, nRegGrid)
  for(i in 1:nRegGrid){
    for(j in 1:nRegGrid){
      true_cov_mat[i,j] <- gamma_fun(seq(a,b,len=nRegGrid)[i], seq(a,b,len=nRegGrid)[j])
    }
  }
  ## ##############################################################
  ## From matrix to list:
  U_list      <- split(U_mat,      rep(1:ncol(U_mat),      each = nrow(U_mat)))
  Y_norm_list <- split(Y_norm_mat, rep(1:ncol(Y_norm_mat), each = nrow(Y_norm_mat)))
  Y_expd_list <- split(Y_expd_mat, rep(1:ncol(Y_expd_mat), each = nrow(Y_expd_mat)))
  
  ##
  reconst_norm_obj <- ReconstPoFD::reconstruct(Ly         = Y_norm_list, 
                                               Lu         = U_list,
                                               K          = NULL,
                                               K_max      = 4,
                                               pre_smooth = FALSE,
                                               nRegGrid   = nRegGrid)
  ##
  reconst_expd_obj <- ReconstPoFD::reconstruct(Ly         = Y_expd_list, 
                                               Lu         = U_list,
                                               K          = NULL,
                                               K_max      = 4,
                                               pre_smooth = FALSE,
                                               nRegGrid   = nRegGrid)
  
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
