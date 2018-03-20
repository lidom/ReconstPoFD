#' Simulate Data 
#'
#' This function allows to simulate partially observed functional data. 
#' @param n          Number of functions
#' @param m          Number of discretization points 
#' @param a          Lower interval boundary
#' @param b          Upper interval boundary
#' @param n_basis    Number of basis functions
#' @param DGP        Data Generating Process. DGP1: Gaussian scores. DGP2: Exponential scores. 
#' @export simuldata
#' @examples  
#' a <- 0; b <- 10; m <- 15; n <- 100
#' mean_fun <- function(u){return( ((u-a)/(b-a)) +sin((u-a)/(b-a)))}
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
simuldata <- function(n = 100, m = 15, a = 0, b = 1, n_basis = 10, DGP=c('DGP1','DGP2')[1]){
  ##
  ## meanfunction
  mean_fun <- function(u){return(((u-a)/(b-a)) + sin((u-a)/(b-a)))}
  eps_var  <- .20
  ##
  ## Generation of prediction points U
  U_true_mat    <- matrix(seq(a,b,len=75), 75, n)
  U_mat         <- matrix(NA, m, n)
  for(i in 1:n){
    ## Random observed interval
    if(1 == stats::rbinom(n = 1, size = 1, prob = .6)){
      A_i  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * 0.33))
      B_i  <- stats::runif(n = 1, min = (b- (b-a) * 0.33), max = b)
    }else{A_i = a; B_i = b}
    ##
    ## sampling from the total grid
    U_mat[,i] <- stats::runif(n=m, min = A_i, max = B_i)
    ## ordering
    U_mat[,i] <- unique(U_mat[,i][order(U_mat[,i])])
  }
  ##
  Y_true_mat <- matrix(NA, 75, n)
  Y_mat      <- matrix(NA, m, n)
  k_vec      <- 1:n_basis
  for(i in 1:n){
    if(DGP=="DGP1"){
      xi1 <- stats::rnorm(n=n_basis, mean=0, sd=sqrt(4-(4/(n_basis + 1))*(k_vec - 1)))
      xi2 <- stats::rnorm(n=n_basis, mean=0, sd=sqrt(4-(4/(n_basis + 1))*(k_vec)))
    }
    if(DGP=="DGP2"){
      xi1 <- c(stats::rexp(n=n_basis, rate=1/sqrt(4-(4/(n_basis + 1))*(k_vec-1)))-sqrt(4-(4/(n_basis + 1))*(k_vec-1)))
      xi2 <- c(stats::rexp(n=n_basis, rate=1/sqrt(4-(4/(n_basis + 1))*(k_vec)))  -sqrt(4-(4/(n_basis + 1))*(k_vec))) 
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
              "Y_list"      = split(Y_mat, rep(1:ncol(Y_mat), each = nrow(Y_mat))), 
              "U_list"      = split(U_mat, rep(1:ncol(U_mat), each = nrow(U_mat))),
              "Y_true_list" = split(Y_true_mat, rep(1:ncol(Y_true_mat), each = nrow(Y_true_mat))), 
              "U_true_list" = split(U_true_mat, rep(1:ncol(U_true_mat), each = nrow(U_true_mat)))
  ))
}

#' Simulate Data as in Kraus (JRSSB, 2015)
#'
#' This function allows to simulate functional data as used in the simulation study of: 
#' Kraus, D. (2015). Components and completion of partially observed functional data. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(4), 777-801. 
#' However, the missingness process is adjusted to avoid missing 'holes'.
#' @param p          Number of grid points in [a,b]
#' @param n          Number of functions
#' @param a          Lower interval boundary
#' @param b          Upper interval boundary
#' @export simuldataKraus
#' @examples  
#' a <- 0; b <- 10; p <- 51; n <- 100
#' SimDat   <- simuldataKraus(p = p, n = n, a = a, b = b)
#' ## 
#' Y_mat    <- SimDat[['Y_mat']]
#' ##
#' matplot(Y_mat[,1:10], col=gray(.5), type="l")
simuldataKraus <- function(p=100, n=100, a=0, b=1){
  ##
  u_vec  <- seq(a, b, len=p)
  ##
  Y_mat  <- matrix(NA, p, n)
  U_mat  <- matrix(NA, p, n)
  ##
  Y_list <- vector("list", n)
  U_list <- vector("list", n)
  ##
  for(i in 1:n){
    Y_vec <- c(rowMeans(
      sapply(1:100,function(k){
        sqrt(3^(-2*k+1)) * stats::rnorm(1) * sqrt(2) * cos(2*pi*k*(u_vec-a)/(b-a)) +
          sqrt(3^(-2*k)) * stats::rnorm(1) * sqrt(2) * sin(2*pi*k*(u_vec-a)/(b-a)) 
      })))
    ## Random subintervals
    if(1 == stats::rbinom(n = 1, size = 1, prob = .6)){
      A_i  <- stats::runif(n = 1, min = a, max = (a+ (b-a) * 0.25))
      B_i  <- stats::runif(n = 1, min = (b- (b-a) * 0.25), max = b)
    }else{A_i = a; B_i = b}
    Y_vec[u_vec < A_i] <- NA
    Y_vec[u_vec > B_i] <- NA
    ##
    U_vec <- seq(a, b, len=p)
    U_vec[u_vec < A_i] <- NA
    U_vec[u_vec > B_i] <- NA
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
    Y_mat[,i]   <- Y_vec
    U_mat[,i]   <- U_vec
    ##
    Y_list[[i]] <- c(stats::na.omit(Y_vec))
    U_list[[i]] <- c(stats::na.omit(U_vec))
  }
  return(list("Y_mat"  = Y_mat, 
              "U_mat"  = U_mat,
              "Y_list" = Y_list, 
              "U_list" = U_list))
}