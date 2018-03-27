#' Simulate Data 
#'
#' This function allows to simulate partially observed functional data. 
#' @param n                 Number of functions
#' @param m                 Number of discretization points 
#' @param a                 Lower interval boundary
#' @param b                 Upper interval boundary
#' @param n_basis           Number of basis functions
#' @param DGP               Data Generating Process. DGP1: Gaussian scores. DGP2: Exponential scores. 
#' @param nRegGrid          Number of grid-points used for the equidistant 'workGrid'.
#' @param determ_obs_interv Set a deterministic interval for the observed part. Default (determ_obs_interv = NULL) means random intervals.
#' @export simuldata
#' @examples  
#' a <- 0; b <- 1; m <- 15; n <- 100
#' mean_fun <- function(u){return( ((u-a)/(b-a)) + 2*sin(2*pi*((u-a)/(b-a)) ))}
#' SimDat   <- simuldata(n = n, m = m, a = a, b = b)
#' ## 
#' Y_mat       <- SimDat[['Y_mat']]
#' U_mat       <- SimDat[['U_mat']]
#' Y_true_mat  <- SimDat[['Y_true_mat']]
#' U_true_mat  <- SimDat[['U_true_mat']]
#' ##
#' par(mfrow=c(2,1))
#' matplot(x=U_mat[,1:5], y=Y_mat[,1:5], col=gray(.5), type="l", main="Missings & Noise")
#' lines(x=U_true_mat[,1], y=mean_fun(U_true_mat[,1]), col="red")
#' matplot(x=U_true_mat[,1:5], y=Y_true_mat[,1:5], col=gray(.5), type="l", main="NoMissings & NoNoise")
#' lines(x=U_true_mat[,1], y=mean_fun(U_true_mat[,1]), col="red")
#' par(mfrow=c(1,1))
simuldata <- function(n = 100, m = 15, a = 0, b = 1, n_basis = 10, DGP=c('DGP1','DGP2')[1], nRegGrid = 75, determ_obs_interv = NULL){
  ##
  ## meanfunction
  mean_fun <- function(u){return( ((u-a)/(b-a)) + 2*sin(2*pi*((u-a)/(b-a)) ))}
  eps_var  <- .20
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
      if(1 == stats::rbinom(n = 1, size = 1, prob = .6)){
        A_vec[i]  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * 0.33))
        B_vec[i]  <- stats::runif(n = 1, min = (b- (b-a) * 0.33), max = b)
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
  }
  ##
  Y_true_mat <- matrix(NA, nRegGrid, n)
  Y_mat      <- matrix(NA, m, n)
  k_vec      <- 1:n_basis
  for(i in 1:n){
    if(DGP=="DGP1"){
      xi1 <- sqrt(10-(10/(n_basis + 1))*(k_vec - 1)) * stats::rnorm(n=n_basis)
      xi2 <- sqrt(10-(10/(n_basis + 1))*(k_vec))     * stats::rnorm(n=n_basis)
    }
    if(DGP=="DGP2"){
      xi1 <- c(stats::rexp(n=n_basis, rate=1/sqrt(10-(10/(n_basis + 1))*(k_vec-1))) - sqrt(10-(10/(n_basis + 1))*(k_vec-1)))
      xi2 <- c(stats::rexp(n=n_basis, rate=1/sqrt(10-(10/(n_basis + 1))*(k_vec)))   - sqrt(10-(10/(n_basis + 1))*(k_vec))) 
    }
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
  return(list("Y_mat"       = Y_mat, 
              "U_mat"       = U_mat,
              "Y_true_mat"  = Y_true_mat, 
              "U_true_mat"  = U_true_mat,
              ##
              "Y_list"      = split(Y_mat, rep(1:ncol(Y_mat), each = nrow(Y_mat))), 
              "U_list"      = split(U_mat, rep(1:ncol(U_mat), each = nrow(U_mat))),
              "Y_true_list" = split(Y_true_mat, rep(1:ncol(Y_true_mat), each = nrow(Y_true_mat))), 
              "U_true_list" = split(U_true_mat, rep(1:ncol(U_true_mat), each = nrow(U_true_mat))),
              ##
              "A_vec"       = A_vec,
              "B_vec"       = B_vec
  ))
}

#' Simulate Data (Adapted from Kraus (2015))
#'
#' This function allows to simulate functional data as used in the simulation study of: 
#' 
#' However, the missingness process is adjusted to avoid missing 'holes'.
#' @param n          Number of functions
#' @param a          Lower interval boundary
#' @param b          Upper interval boundary
#' @param DGP        Data Generating Process. DGP3: The DGP of Kraus. DGP4: Similar to the DGP of Kraus, but with a mean function and slower decaying eigenvalues.
#' @param nRegGrid    Number of grid-points used for the equidistant 'workGrid'.
#' @param determ_obs_interv Set a deterministic interval for the observed part. Default (determ_obs_interv = NULL) means random intervals.
#' @export simuldataKraus
#' @references 
#' Kraus, D. (2015). Components and completion of partially observed functional data. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(4), 777-801. 
#' @examples  
#' a <- 0; b <- 1; n <- 100; DGP <- c('DGP3', 'DGP4')[2]
#' if(DGP=='DGP3'){
#'   mean_fun <- function(u){return(rep(0,length(u)))}
#' }
#' if(DGP=='DGP4'){
#'   mean_fun <- function(u){return( ((u-a)/(b-a))^2 + cos(3*pi*((u-a)/(b-a)) ))}
#' }
#' SimDat   <- simuldataKraus(n = n, a = a, b = b, DGP = DGP)
#' ## 
#' Y_mat       <- SimDat[['Y_mat']]
#' U_mat       <- SimDat[['U_mat']]
#' Y_true_mat  <- SimDat[['Y_true_mat']]
#' U_true_mat  <- SimDat[['U_true_mat']]
#' ##
#' par(mfrow=c(2,1))
#' matplot(x=U_mat[,1:5], y=Y_mat[,1:5], col=gray(.5), type="l", 
#' main="Missings & Noise", xlim=c(a,b))
#' lines(x=U_true_mat[,1], y=mean_fun(U_true_mat[,1]), col="red")
#' matplot(x=U_true_mat[,1:5], y=Y_true_mat[,1:5], col=gray(.5), type="l", 
#' main="NoMissings & NoNoise", xlim=c(a,b))
#' lines(x=U_true_mat[,1], y=mean_fun(U_true_mat[,1]), col="red")
#' par(mfrow=c(1,1))
simuldataKraus <- function(n=100, a=0, b=1, DGP=c('DGP3','DGP4')[1], nRegGrid = 75, determ_obs_interv = NULL)
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
    if(DGP=='DGP3'){
      mean_fun <- function(u){return(rep(0,length(u)))}
      ##
      xi1   <- sqrt(3^(-2*k_vec+1)) * stats::rnorm(n = length(k_vec))
      xi2   <- sqrt(3^(-2*k_vec))   * stats::rnorm(n = length(k_vec))
    }
    if(DGP=='DGP4'){
      mean_fun <- function(u){return( ((u-a)/(b-a))^2 + cos(3*pi*((u-a)/(b-a)) ))}
      ##
      xi1   <- 25*sqrt(25^(-(k_vec+1)/5)) * stats::rnorm(n = length(k_vec))
      xi2   <- 25*sqrt(25^(-(k_vec  )/5)) * stats::rnorm(n = length(k_vec))
    }
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
      if(1 == stats::rbinom(n = 1, size = 1, prob = .6)){
        A_vec[i]  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * 0.33))
        B_vec[i]  <- stats::runif(n = 1, min = (b- (b-a) * 0.33), max = b)
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
              "B_vec"       = B_vec
  ))
}



