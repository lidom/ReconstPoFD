#' Reconstruct functional data as proposed by Kraus (JRSSB, 2015)
#'
#' This function allows to reconstruct partially observed functional data as proposed in: 
#' Kraus, D. (2015). Components and completion of partially observed functional data. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(4), 777-801. 
#' @param X_mat      pxn matrix (p: number of discretization points, n=number of functions
#' @param alpha      Ridge parameter. If alpha = NULL (default), an optimal alpha is determined by GCV
#' @export reconstructKraus
#' @examples  
#' a <- 0; b <- 10; p <- 51; n <- 100
#' SimDat   <- simuldataKraus(p = p, n = n, a = a, b = b)
#' ## 
#' Y_mat <- SimDat[['Y_mat']]
#' ##
#' reconst_mat <- reconstructKraus(Y_mat)
#' ##
#' par(mfrow=c(2,1))
#' matplot(Y_mat[,1:5],       col=gray(.5), type="l", main="Original Data")
#' matplot(reconst_mat[,1:5], col=gray(.5), type="l")
#' par(mfrow=c(1,1))
reconstructKraus <- function(X_mat, alpha = NULL){
  ##
  mean_vec      <- meanKraus(X_mat)
  cov_mat       <- covKraus(X_mat)
  n             <- ncol(X_mat)
  X_reconst_mat <- X_mat
  ##
  NonNA_fcts    <- apply(X_mat,2,function(x)!any(is.na(x)))
  X_Compl_mat   <- X_mat[,NonNA_fcts]
  ##
  for(i in 1:n){
    X_tmp      <- X_mat[,i]
    ##
    M_bool_vec <- is.na(X_tmp)
    O_bool_vec <- !M_bool_vec
    ##
    if(is.null(alpha)){
      alpha <- optimize(f = function(alpha){gcvKraus(cov_mat     = cov_mat, 
                                                     mean_vec    = mean_vec, 
                                                     X_Compl_mat = X_Compl_mat, 
                                                     M_bool_vec  = M_bool_vec, 
                                                     alpha       = alpha)},
                        interval = c(.Machine$double.eps, sum(diag(cov_mat))*n), maximum = FALSE)$minimum
    }
    ##
    result_tmp <- reconstKraus_fun(cov_mat    = cov_mat, 
                                   X_cent_vec = c(X_tmp - mean_vec), 
                                   alpha      = alpha)
    ##
    X_reconst_mat[,i]  <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
  }
  return(X_reconst_mat)
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
#' @param missings   Generate functions with missing parts (NAs); default: missings = TRUE
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
        sqrt(3^(-2*k+1)) * rnorm(1) * sqrt(2) * cos(2*pi*k*(u_vec-a)/(b-a)) +
          sqrt(3^(-2*k)) * rnorm(1) * sqrt(2) * sin(2*pi*k*(u_vec-a)/(b-a)) 
      })))
    if(missings){
      ## Random subintervals
      if(1 == rbinom(n = 1, size = 1, prob = .6)){
        A_i  <- runif(n = 1, min = a, max = (a+ (b-a) * 0.25))
        B_i  <- runif(n = 1, min = (b- (b-a) * 0.25), max = b)
      }else{A_i = a; B_i = b}
      Y_vec[u_vec < A_i] <- NA
      Y_vec[u_vec > B_i] <- NA
      ##
      U_vec <- seq(a, b, len=p)
      U_vec[u_vec < A_i] <- NA
      U_vec[u_vec > B_i] <- NA
      ## 
      # ## Kraus-Version (contains 'holes' that cannot be used with our method so far)
      # d=1.4; f=0.2
      # U_12 <- runif(2)
      # C    <- d*sqrt(U_12[1])
      # E    <- f*U_12[2]
      # M_lo <- max(0,C-E)
      # M_up <- min(1,C+E)
      # X_tmp[M_lo <= t_vec & t_vec <= M_up] <- NA
      }
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


meanKraus <- function(X_mat){
  rowMeans(X_mat, na.rm = TRUE)
}

covKraus <- function(X_mat){
  p <- nrow(X_mat)
  n <- ncol(X_mat)
  ##
  X_cent_mat  <- X_mat - rowMeans(X_mat, na.rm = TRUE)
  ##
  covKraus_mat <- matrix(NA, ncol = p, nrow = p)	
  for(s in seq(1, p)){
    for(t in seq(s, p)){
      X_cent_s  <- X_cent_mat[s, ]
      X_cent_t  <- X_cent_mat[t, ]
      n_na      <- sum(is.na(c(X_cent_s * X_cent_t)))
      if(n-n_na == 0){
        covKraus_mat[s,t] <- NA
      }else{
        covKraus_mat[s,t] <- mean(X_cent_s * X_cent_t, na.rm = TRUE)
      }
      covKraus_mat[t,s] <- covKraus_mat[s,t]
    }
  }
  return(covKraus_mat)
}


reconstKraus_fun <- function(cov_mat, X_cent_vec, alpha=1e-4){
  ##
  M_bool_vec   <- is.na(X_cent_vec)
  O_bool_vec   <- !M_bool_vec
  p            <- nrow(cov_mat)
  ##
  covMO_mat    <- cov_mat[M_bool_vec, O_bool_vec] 
  covOO_mat    <- cov_mat[O_bool_vec, O_bool_vec]
  ##
  covOO_a_mat      <- covOO_mat + alpha * diag(p)[O_bool_vec, O_bool_vec]
  covOO_a_mat_inv  <- solve(covOO_a_mat)
  Aa_mat           <- covMO_mat %*% covOO_a_mat_inv
  ##
  X_M_fit_cent_vec <- Aa_mat %*% X_cent_vec[O_bool_vec]
  ##
  X_cent_reconst_vec             <- X_cent_vec
  X_cent_reconst_vec[M_bool_vec] <- X_M_fit_cent_vec
  ##
  df <- sum(diag(covOO_a_mat_inv %*% covOO_mat))
  ##
  return(list("X_cent_reconst_vec" = X_cent_reconst_vec,
              "df"                 = df))
}

gcvKraus <- function(cov_mat, mean_vec, X_Compl_mat, M_bool_vec, alpha){
  n_Compl  <- ncol(X_Compl_mat)
  rss_vec  <- rep(NA,n_Compl)
  ##
  for(j in 1:n_Compl){
    X_gcv             <-  X_Compl_mat[,j]
    X_gcv[M_bool_vec] <- NA
    ##
    result_tmp <- reconstKraus_fun(cov_mat    = cov_mat, 
                                   X_cent_vec = c(X_gcv - mean_vec), 
                                   alpha      = alpha)
    ##
    X_fit      <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
    ##
    rss_vec[j] <- sum((X_fit[M_bool_vec] - X_Compl_mat[M_bool_vec,j])^2)
  }
  gcv <- sum(rss_vec)/((1-result_tmp[['df']]/n_Compl)^2)
  ##
  return(gcv)
}
gcvKraus <- Vectorize(FUN = gcvKraus, vectorize.args = "alpha")

