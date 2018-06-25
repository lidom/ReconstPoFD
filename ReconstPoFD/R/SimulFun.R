#' Simulate Data 
#'
#' This function allows to simulate partially observed functional data. 
#' @param n                 Number of functions
#' @param m                 Number of discretization points 
#' @param a                 Lower interval boundary
#' @param b                 Upper interval boundary
#' @param DGP               Data Generating Process. DGP1: Gaussian scores. DGP2: Exponential scores. 
#' @param nRegGrid          Number of grid-points used for the equidistant 'workGrid'.
#' @param determ_obs_interv Set a deterministic interval for the observed part. Default (determ_obs_interv = NULL) means random intervals.
#' @export simuldata
#' @examples
#' DGP=c('DGP1','DGP2','DGP3')[1]
#' SimDat        <- simuldata(n = 50, m = 15, a = 0, b = 1, DGP=DGP)
#' Y_mat         <- SimDat[['Y_mat']]
#' U_mat         <- SimDat[['U_mat']]
#' Y_true_mat    <- SimDat[['Y_true_mat']]
#' U_true_mat    <- SimDat[['U_true_mat']]
#' mean_true_vec <- SimDat[['mean_true_vec']]
#' ##
#' par(mfrow=c(2,1))
#' matplot(x=U_mat[,1:5], y=Y_mat[,1:5], col=gray(.5), type="l", main="Missings")
#' lines(  x=U_true_mat[,1], y=mean_true_vec, col="red")
#' matplot(x=U_true_mat[,1:5], y=Y_true_mat[,1:5], col=gray(.5), type="l",
#' main=ifelse(DGP=="DGP1", "No Missings & No Noise", "No Missings"))
#' lines(  x=U_true_mat[,1],   y=mean_true_vec, col="red")
#' par(mfrow=c(1,1))
simuldata <- function(n = 100, m = 15, a = 0, b = 1, DGP=c('DGP1','DGP2','DGP3')[1], nRegGrid = 51, determ_obs_interv = NULL){
  if(DGP=="DGP1"){SimDat <- simuldata_1(n=n,m=m,a=a,b=b,nRegGrid=nRegGrid,determ_obs_interv=determ_obs_interv)}
  if(DGP=="DGP2"){SimDat <- simuldata_2(n=n,    a=a,b=b,nRegGrid=nRegGrid,determ_obs_interv=determ_obs_interv)}
  if(DGP=="DGP3"){SimDat <- simuldata_3(n=n,    a=a,b=b,nRegGrid=nRegGrid,determ_obs_interv=determ_obs_interv)}
  return(SimDat)
}

##-------------------------------------------------------------------------------------
#simuldataDGP_1_2 <- function(n = 100, m = 15, a = 0, b = 1, DGP=c('DGP1','DGP2')[1], nRegGrid = 51, determ_obs_interv = NULL){
simuldata_1 <- function(n = 100, m = 15, a = 0, b = 1, nRegGrid = 51, determ_obs_interv = NULL){
  ##
  ## meanfunction
  mean_fun <- function(u){return( ((u-a)/(b-a)) + 1*sin(2*pi*((u-a)/(b-a))) )}
  eps_var  <- .025
  n_basis  <-  50
  ##
  ## Generation of prediction points U
  U_true_mat    <- matrix(seq(a,b,len=nRegGrid), nRegGrid, n)
  U_mat         <- matrix(NA, m, n)
  A_vec         <- rep(NA, n)
  B_vec         <- rep(NA, n)
  ##
  for(i in 1:n){
    if(is.null(determ_obs_interv)){
      ## Random observed interval
      if(1 == stats::rbinom(n = 1, size = 1, prob = .75)){
        A_vec[i]  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * .45))
        B_vec[i]  <- stats::runif(n = 1, min = (b- (b-a) * .45), max = b)
      }else{
        A_vec[i] = a
        B_vec[i] = b
      }
    }else{
      A_vec[i] = determ_obs_interv[1]
      B_vec[i] = determ_obs_interv[2]
    }
    ##
    ## sampling 
    U_mat[,i] <- stats::runif(n=m, min = A_vec[i], max = B_vec[i])
    ## ordering
    U_mat[,i] <- unique(U_mat[,i][order(U_mat[,i])])
    U_mat[,i][which.min(U_mat[,i])] <- A_vec[i]
    U_mat[,i][which.max(U_mat[,i])] <- B_vec[i]
  }
  ##
  Y_true_mat <- matrix(NA, nRegGrid, n)
  Y_mat      <- matrix(NA, m, n)
  k_vec      <- 1:n_basis
  for(i in 1:n){
#    if(DGP=="DGP1"){
      xi1 <- 25*sqrt(exp(-((k_vec-1)^2)/5)) * stats::rnorm(n=n_basis)
      xi2 <- 25*sqrt(exp(-((k_vec  )^2)/5)) * stats::rnorm(n=n_basis)
#    }
    # plot(y=c(sqrt(10-(10/(n_basis + 1))*(k_vec - 1))), x=k_vec)
#    if(DGP=="DGP2"){
#      xi1 <- c(stats::rexp(n=n_basis, rate=1/sqrt(10-(10/(n_basis + 1))*(k_vec-1))) - sqrt(10-(10/(n_basis + 1))*(k_vec-1)))
#      xi2 <- c(stats::rexp(n=n_basis, rate=1/sqrt(10-(10/(n_basis + 1))*(k_vec)))   - sqrt(10-(10/(n_basis + 1))*(k_vec))) 
#    }
    ##
    Y_true_mat[,i] <- c(c(rowMeans(
      sapply(k_vec, function(k){
        xi1[k] * -1 * cos(k*pi*(U_true_mat[,i]-a)/(b-a))/sqrt(5) +
          xi2[k] *      sin(k*pi*(U_true_mat[,i]-a)/(b-a))/sqrt(5) 
      }))) + mean_fun(U_true_mat[,i])) 
    ##
    Y_mat[,i] <- c(c(rowMeans(
      sapply(k_vec, function(k){
        xi1[k] * -1 * cos(k*pi*(U_mat[,i]-a)/(b-a))/sqrt(5) +
          xi2[k] *      sin(k*pi*(U_mat[,i]-a)/(b-a))/sqrt(5) 
      }))) + mean_fun(U_mat[,i]) + stats::rnorm(n=m,  mean = 0, sd=sqrt(eps_var)))
    ##
  }
  return(list("Y_mat"         = Y_mat, 
              "U_mat"         = U_mat,
              "Y_true_mat"    = Y_true_mat, 
              "U_true_mat"    = U_true_mat,
              ##
              "Y_list"        = split(Y_mat, rep(1:ncol(Y_mat), each = nrow(Y_mat))), 
              "U_list"        = split(U_mat, rep(1:ncol(U_mat), each = nrow(U_mat))),
              "Y_true_list"   = split(Y_true_mat, rep(1:ncol(Y_true_mat), each = nrow(Y_true_mat))), 
              "U_true_list"   = split(U_true_mat, rep(1:ncol(U_true_mat), each = nrow(U_true_mat))),
              ##
              "A_vec"         = A_vec,
              "B_vec"         = B_vec,
              ##
              "mean_true_vec" = mean_fun(U_true_mat[,1])
  ))
}

##-------------------------------------------------------------------------------------
# simuldataKraus <- function(n=100, a=0, b=1, DGP=c('DGP3','DGP4')[1], nRegGrid = 51, determ_obs_interv = NULL)
simuldata_2 <- function(n=100, a=0, b=1, nRegGrid = 51, determ_obs_interv = NULL)
{
  ##
  ## Number of grid points in [a,b]
  p <- nRegGrid
  ##
  Y_mat       <- matrix(NA, p, n)
  U_mat       <- matrix(NA, p, n)
  Y_true_mat  <- matrix(NA, p, n)
  U_true_mat  <- matrix(NA, p, n)
  ##
  Y_list      <- vector("list", n)
  U_list      <- vector("list", n)
  Y_true_list <- vector("list", n)
  U_true_list <- vector("list", n)
  ##
  k_vec <- 1:100
  ##
  A_vec         <- rep(NA, n)
  B_vec         <- rep(NA, n)
  ##
  for(i in 1:n){
#    if(DGP=='DGP3'){
#      mean_fun <- function(u){return(rep(1,length(u)))}
#      ##
#      xi1   <- 9*sqrt(3^(-2*(k_vec+1))) * stats::rnorm(n = length(k_vec))
#      xi2   <- 9*sqrt(3^(-2*(k_vec  ))) * stats::rnorm(n = length(k_vec))
#    }
#    if(DGP=='DGP4'){
      mean_fun <- function(u){return( ((u-a)/(b-a))^2 + cos(3*pi*((u-a)/(b-a))) )}
      ##
      xi1   <- 25*sqrt(exp(-((k_vec-1)^2)/5)) * stats::rnorm(n = length(k_vec))
      xi2   <- 25*sqrt(exp(-((k_vec  )^2)/5)) * stats::rnorm(n = length(k_vec))
#    } 
    ##
    U_vec  <- seq(a, b, len=p)
    ##
    Y_vec <- c(rowMeans(
      sapply(k_vec, function(k){
        xi1[k] * sqrt(2) * cos(2*pi*k*(U_vec-a)/(b-a)) +
          xi2[k] * sqrt(2) * sin(2*pi*k*(U_vec-a)/(b-a)) 
      }))) + mean_fun(U_vec)
    ##
    if(is.null(determ_obs_interv)){
      ## Random observed interval
      if(1 == stats::rbinom(n = 1, size = 1, prob = .75)){
        A_vec[i]  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * .45))
        B_vec[i]  <- stats::runif(n = 1, min = (b- (b-a) * .45), max = b)
      }else{
        A_vec[i] = a
        B_vec[i] = b
      }
    }else{
      A_vec[i] = determ_obs_interv[1]
      B_vec[i] = determ_obs_interv[2]
    }
    ##
    Y_true_vec              <- Y_vec
    Y_vec[U_vec < A_vec[i]] <- NA
    Y_vec[U_vec > B_vec[i]] <- NA
    ##
    U_true_vec                   <- U_vec
    U_vec[U_true_vec < A_vec[i]] <- NA
    U_vec[U_true_vec > B_vec[i]] <- NA
    ##------------------------- 
    # ## Kraus-Version (contains 'holes' that cannot be used with our method so far)
    # d=1.4; f=0.2
    # U_12 <- stats::runif(2)
    # C    <- d*sqrt(U_12[1])
    # E    <- f*U_12[2]
    # M_lo <- max(0,C-E)
    # M_up <- min(1,C+E)
    # X_tmp[M_lo <= t_vec & t_vec <= M_up] <- NA
    ##-------------------------
    Y_mat[,i]        <- Y_vec
    U_mat[,i]        <- U_vec
    Y_true_mat[,i]   <- Y_true_vec
    U_true_mat[,i]   <- U_true_vec
    ##
    Y_list[[i]]      <- c(stats::na.omit(Y_vec))
    U_list[[i]]      <- c(stats::na.omit(U_vec))
    Y_true_list[[i]] <- Y_true_vec
    U_true_list[[i]] <- U_true_vec
  }
  return(list("Y_mat"  = Y_mat, 
              "U_mat"  = U_mat,
              "Y_list" = Y_list, 
              "U_list" = U_list,
              ##
              "Y_true_mat"  = Y_true_mat, 
              "U_true_mat"  = U_true_mat,
              "Y_true_list" = Y_true_list, 
              "U_true_list" = U_true_list,
              ##
              "A_vec"       = A_vec,
              "B_vec"       = B_vec,
              ##
              "mean_true_vec" = mean_fun(U_true_mat[,1])
  ))
}

##-------------------------------------------------------------------------------------
# simuldataWBF <- function(n=100, a=0, b=1, DGP=c('DGP5'), nRegGrid = 51, determ_obs_interv = NULL)
simuldata_3 <- function(n=100, a=0, b=1, nRegGrid = 51, determ_obs_interv = NULL)  
{
  ##
  ## Number of grid points in [a,b]
  p <- nRegGrid
  ##
  Y_mat       <- matrix(NA, p, n)
  U_mat       <- matrix(NA, p, n)
  Y_true_mat  <- matrix(NA, p, n)
  U_true_mat  <- matrix(NA, p, n)
  ##
  Y_list      <- vector("list", n)
  U_list      <- vector("list", n)
  Y_true_list <- vector("list", n)
  U_true_list <- vector("list", n)
  ##
  A_vec         <- rep(NA, n)
  B_vec         <- rep(NA, n)
  ##
  for(i in 1:n){
#    if(DGP=='DGP5'){
      mean_fun <- function(u){return( 10 - 5*((u-a)/(b-a) - 0.5)^3  )}
      ##
      rand_vec <- stats::rnorm(n = 3, sd=1)
#    }
    ##
    loc_vec <- seq(a, b, len=length(rand_vec))
    U_vec   <- seq(a, b, len=p)
    ##
    Y_vec <- stats::spline(x = loc_vec, y = rand_vec, method = "natural", xout = U_vec)$y + mean_fun(U_vec)
    ##
    if(is.null(determ_obs_interv)){
      ## Random observed interval
      if(1 == stats::rbinom(n = 1, size = 1, prob = .75)){
        A_vec[i]  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * .45))
        B_vec[i]  <- stats::runif(n = 1, min = (b- (b-a) * .45), max = b)
      }else{
        A_vec[i] = a
        B_vec[i] = b
      }
    }else{
      A_vec[i] = determ_obs_interv[1]
      B_vec[i] = determ_obs_interv[2]
    }
    ##
    Y_true_vec              <- Y_vec
    Y_vec[U_vec < A_vec[i]] <- NA
    Y_vec[U_vec > B_vec[i]] <- NA
    ##
    U_true_vec                   <- U_vec
    U_vec[U_true_vec < A_vec[i]] <- NA
    U_vec[U_true_vec > B_vec[i]] <- NA
    ##
    Y_mat[,i]        <- Y_vec
    U_mat[,i]        <- U_vec
    Y_true_mat[,i]   <- Y_true_vec
    U_true_mat[,i]   <- U_true_vec
    ##
    Y_list[[i]]      <- c(stats::na.omit(Y_vec))
    U_list[[i]]      <- c(stats::na.omit(U_vec))
    Y_true_list[[i]] <- Y_true_vec
    U_true_list[[i]] <- U_true_vec
  }
  return(list("Y_mat"  = Y_mat, 
              "U_mat"  = U_mat,
              "Y_list" = Y_list, 
              "U_list" = U_list,
              ##
              "Y_true_mat"  = Y_true_mat, 
              "U_true_mat"  = U_true_mat,
              "Y_true_list" = Y_true_list, 
              "U_true_list" = U_true_list,
              ##
              "A_vec"       = A_vec,
              "B_vec"       = B_vec,
              ##
              "mean_true_vec" = mean_fun(U_true_mat[,1])
  ))
}



