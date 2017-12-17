## R-packages 
library("devtools")
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("doParallel", "fdapace")
library("ReconstPoFD")  # contains the function 'reconst_fun()'
library("doParallel")   # parallel-looping
library("fdapace")      # for estimating mean- and covariance-function


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
(K_irue  <- length(lam_vec[lam_vec>0]))
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
  U_list      <- split(U_mat, 
                       rep(1:ncol(U_mat), 
                           each = nrow(U_mat)))
  Y_norm_list <- split(Y_norm_mat, 
                       rep(1:ncol(Y_norm_mat), 
                           each = nrow(Y_norm_mat)))
  Y_expd_list <- split(Y_expd_mat, 
                       rep(1:ncol(Y_expd_mat), 
                           each = nrow(Y_expd_mat)))
  
  ## Normal #############################################################################
  FPCA_norm_obj <- FPCA(Ly    = Y_norm_list, 
                        Lt    = U_list, 
                        optns = list(
                          "dataType"      = "Sparse", 
                          "kernel"        = "gauss",
                          "error"         = TRUE,
                          "nRegGrid"      = nRegGrid))
  ## Mixture ############################################################################
  FPCA_expd_obj <- FPCA(Ly     = Y_expd_list, 
                        Lt     = U_list, 
                        optns = list(
                          "dataType"      = "Sparse", 
                          "kernel"        = "gauss",
                          "error"         = TRUE,
                          "nRegGrid"      = nRegGrid))
  ## 
  cov_norm_la_mat     <- FPCA_norm_obj$smoothedCov
  cov_expd_la_mat     <- FPCA_expd_obj$smoothedCov
  grid_la_vec         <- FPCA_norm_obj$workGrid
  
  
  ## Centering the data ###############################################
  mu_norm_est_fun <- splinefun(y= FPCA_norm_obj$mu, x=grid_la_vec)
  mu_expd_est_fun <- splinefun(y= FPCA_expd_obj$mu, x=grid_la_vec)
  ## 
  Y_cent_norm_mat <- matrix(NA, m, n)
  Y_cent_expd_mat <- matrix(NA, m, n)
  for(i in 1:n){
    Y_cent_norm_mat[,i]  <- Y_norm_mat[,i] - mu_norm_est_fun(U_mat[,i])
    Y_cent_expd_mat[,i]  <- Y_expd_mat[,i] - mu_expd_est_fun(U_mat[,i])
  }
  ## ###################################################################
  
  ## Aligning Y_cent according to 'grid_la_vec' ########################
  Y.c.align.norm.mat    <- matrix(NA, nRegGrid, n)
  Y.c.align.expd.mat    <- matrix(NA, nRegGrid, n)
  U.align.mat           <- matrix(NA, nRegGrid, n)
  for(i in 1:n){
    Y.c.norm.tmp    <- c(na.omit(Y_cent_norm_mat[,i]))
    Y.c.expd.tmp    <- c(na.omit(Y_cent_expd_mat[,i]))
    U.tmp           <- c(na.omit(U_mat[,i]))
    for(j in 1:length(c(na.omit(U_mat[,i])))){
      loc                        <- order(abs(grid_la_vec - U.tmp[j]))[1]
      Y.c.align.norm.mat[loc,i]  <- Y.c.norm.tmp[j]
      Y.c.align.expd.mat[loc,i]  <- Y.c.expd.tmp[j]
    }
    na.loc                  <- is.na(Y.c.align.norm.mat[,i])
    U.align.mat[!na.loc,i]  <- grid_la_vec[!na.loc]       
  }
  ## ###################################################################
  
  
  ## ############################## ############################## ############################
  ## Selecting K
  ## ############################## ############################## ############################
  ## #########################################################
  ## Nonparametric variance estimation 
  ## Gasser, Stroka, Jennen-Steinmetz (1986, Biometrika)
  ## #########################################################
  sig2.GSJ.eps.norm.vec <- NULL
  sig2.GSJ.eps.expd.vec <- NULL
  for(i in 1:n){
    x.vec      <- c(na.omit(U.align.mat[,i]))
    y.norm.vec <- c(na.omit(Y.c.align.norm.mat[,i]))
    y.expd.vec <- c(na.omit(Y.c.align.expd.mat[,i]))
    ##
    a.seq             <- (diff(x.vec, lag=1)[-1]/diff(x.vec, lag=2))
    b.seq             <- (diff(x.vec, lag=1)[-length(diff(x.vec, lag=1))]/diff(x.vec, lag=2))
    c.sq              <- 1/(a.seq^2+b.seq^2+1)
    ##
    pseudo.norm.eps   <- a.seq * y.norm.vec[-c( length(y.norm.vec)-1,  length(y.norm.vec))] +
      b.seq * y.norm.vec[-c(1,2)] - y.norm.vec[-c(1,length(y.norm.vec))]
    pseudo.expd.eps   <- a.seq * y.expd.vec[-c( length(y.expd.vec)-1,  length(y.expd.vec))] +
      b.seq * y.expd.vec[-c(1,2)] - y.expd.vec[-c(1,length(y.expd.vec))]
    sig2.GSJ.norm     <- mean(c(pseudo.norm.eps^2*c.sq))
    sig2.GSJ.expd     <- mean(c(pseudo.expd.eps^2*c.sq))
    ##
    sig2.GSJ.eps.norm.vec <- c(sig2.GSJ.eps.norm.vec, sig2.GSJ.norm)
    sig2.GSJ.eps.expd.vec <- c(sig2.GSJ.eps.expd.vec, sig2.GSJ.expd)
  }
  sig2.GSJ.norm.eps <- mean(sig2.GSJ.eps.norm.vec)
  sig2.GSJ.expd.eps <- mean(sig2.GSJ.eps.expd.vec)
  
  
  ## ############################
  ## Selecting K via AIC 
  ## ############################
  K.max    <- 4
  AIC.norm <- rep(NA,K.max)
  AIC.expd <- rep(NA,K.max)
  ##
  for(K in 1:K.max){
    RSS.norm    <- rep(NA,n)
    RSS.expd    <- rep(NA,n)
    L.norm.vec  <- rep(NA,n)
    L.expd.vec  <- rep(NA,n)
    for(i in 1:n){ 
      Y_cent_sm.norm_i <- c(na.omit(Y.c.align.norm.mat[,i]))
      Y_cent_sm.expd_i <- c(na.omit(Y.c.align.expd.mat[,i]))
      U_sm_i           <- c(na.omit(U.align.mat[,i]))
      lo.half          <- 1:floor(length(U_sm_i)/2)
      up.half          <- (floor(length(U_sm_i)/2)+1):length(U_sm_i)
      ##
      List.norm    <- reconst_fun(cov_la_mat     = cov_norm_la_mat, 
                                      domain_grid    = grid_la_vec, 
                                      Y_cent_sm_i    = Y_cent_sm.norm_i[lo.half], 
                                      U_sm_i         = U_sm_i[lo.half], 
                                      K              = K)
      List.expd    <- reconst_fun(cov_la_mat     = cov_expd_la_mat, 
                                      domain_grid    = grid_la_vec, 
                                      Y_cent_sm_i    = Y_cent_sm.expd_i[lo.half], 
                                      U_sm_i         = U_sm_i[lo.half], 
                                      K              = K)
      ##
      y.norm.fit  <- List.norm[['y_reconst']] + mu_norm_est_fun(List.norm[['x_reconst']])
      y.expd.fit  <- List.expd[['y_reconst']] + mu_norm_est_fun(List.expd[['x_reconst']])
      ##
      RSS.norm[i] <- sum((y.norm.fit[!is.na(U.align.mat[,i])][up.half] - (Y_cent_sm.norm_i[up.half] + mu_norm_est_fun(U_sm_i)[up.half]))^2)
      RSS.expd[i] <- sum((y.expd.fit[!is.na(U.align.mat[,i])][up.half] - (Y_cent_sm.expd_i[up.half] + mu_norm_est_fun(U_sm_i)[up.half]))^2)
      ##
      L.norm.vec[i] <- -(n*log(2*pi)/2)-(n*log(sig2.GSJ.norm.eps)/2)-(RSS.norm[i]/(2*sig2.GSJ.norm.eps))
      L.expd.vec[i] <- -(n*log(2*pi)/2)-(n*log(sig2.GSJ.expd.eps)/2)-(RSS.expd[i]/(2*sig2.GSJ.expd.eps))
    }
    AIC.norm[K] <- -sum(L.norm.vec) + K 
    AIC.expd[K] <- -sum(L.expd.vec) + K
  }
  K.norm <- which.min(AIC.norm)
  K.expd <- which.min(AIC.expd)
  
  ## Fitted Covariance:
  if(K.norm>1){
    cov_norm_la_mat <- FPCA_norm_obj$phi[,1:K.norm,drop=FALSE] %*% diag(FPCA_norm_obj$lambda[1:K.norm]) %*%  t(FPCA_norm_obj$phi[,1:K.norm,drop=FALSE])
  }
  if(K.norm==1){
    cov_norm_la_mat <- FPCA_norm_obj$phi[,1:K.norm,drop=FALSE]  %*%  t(FPCA_norm_obj$phi[,1:K.norm,drop=FALSE]) * FPCA_norm_obj$lambda[1:K.norm]
  }
  if(K.expd>1){
    cov_expd_la_mat <- FPCA_expd_obj$phi[,1:K.expd,drop=FALSE] %*% diag(FPCA_expd_obj$lambda[1:K.expd]) %*% t(FPCA_expd_obj$phi[,1:K.expd,drop=FALSE])
  }
  if(K.expd==1){
    cov_expd_la_mat <- FPCA_expd_obj$phi[,1:K.expd,drop=FALSE] %*% t(FPCA_expd_obj$phi[,1:K.expd,drop=FALSE]) * FPCA_expd_obj$lambda[1:K.expd]
  }
  
  ## ##########################  ## ################################  ## ################################
  ## Reconstructing all functions
  ## ##########################  ## ################################  ## ################################
  Y.reconst.norm.list  <- vector("list", n)
  Y.reconst.expd.list  <- vector("list", n)
  X.reconst.norm.list  <- vector("list", n)
  X.reconst.expd.list  <- vector("list", n)
  ##
  MaxAE_norm_vec     <- rep(NA, n)
  MaxAE_expd_vec     <- rep(NA, n)
  for(i in 1:n){
    tmp.norm  <- reconst_fun(cov_la_mat  = cov_norm_la_mat, 
                             domain_grid = grid_la_vec, 
                             Y_cent_sm_i = c(na.omit(Y.c.align.norm.mat[,i])), 
                             U_sm_i      = c(na.omit(U.align.mat[,i])), 
                             K           = K.norm)
    
    tmp.expd  <- reconst_fun(cov_la_mat  = cov_expd_la_mat, 
                             domain_grid = grid_la_vec, 
                             Y_cent_sm_i = c(na.omit(Y.c.align.expd.mat[,i])), 
                             U_sm_i      = c(na.omit(U.align.mat[,i])), 
                             K           = K.expd)
    
    x_algo <- tmp.norm[['x_reconst']]
    y_algo <- tmp.norm[['y_reconst']]
    y_algo <- y_algo + mu_norm_est_fun(x_algo)
    Y.reconst.norm.list[[i]]  <- y_algo
    X.reconst.norm.list[[i]]  <- x_algo
    ##
    x_algo <- tmp.expd[['x_reconst']]
    y_algo <- tmp.expd[['y_reconst']]
    y_algo <- y_algo + mu_expd_est_fun(x_algo)
    Y.reconst.expd.list[[i]]  <- y_algo
    X.reconst.expd.list[[i]]  <- x_algo
    
    ## Max absolute error:
    MaxAE_norm_vec[i] <- max(abs(Y.reconst.norm.list[[i]] - True_norm_curve_fun(u=X.reconst.norm.list[[i]], i=i)))
    MaxAE_expd_vec[i] <- max(abs(Y.reconst.expd.list[[i]] - True_expd_curve_fun(u=X.reconst.expd.list[[i]], i=i)))
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

## Save results:
## save(sim.results, file = paste0("n",n,"_m",m,"_simResults.RData"))
