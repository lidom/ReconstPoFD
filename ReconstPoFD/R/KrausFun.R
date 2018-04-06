#' Reconstruct functional data as proposed by Kraus (JRSSB, 2015)
#'
#' This function allows to reconstruct partially observed functional data as proposed in: 
#' Kraus, D. (2015). Components and completion of partially observed functional data. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(4), 777-801. 
#' @param X_mat        pxn matrix (p: number of discretization points, n=number of functions
#' @param alpha        Ridge parameter. If alpha = NULL (default), an optimal alpha is determined by GCV
#' @param reconst_fcts A vector specifying the list elements in Ly which need to be reconstructed. Default (reconst_fcts=NULL) will reconstruct all functions.
#' @export reconstructKraus
#' @examples  
#' SimDat       <- simuldata(n = 50, a = 0, b = 1, DGP="DGP3")
#' Y_mat        <- SimDat[['Y_mat']]
#' U_mat        <- SimDat[['U_mat']]
#' U_true_mat   <- SimDat[['U_true_mat']]
#' ##
#' result        <- reconstructKraus(X_mat = Y_mat)
#' Y_reconst_mat <- result[['X_reconst_mat']]
#' ##
#' par(mfrow=c(2,1))
#' matplot(x=U_mat[,1:5], y=Y_mat[,1:5], col=gray(.5), type="l", 
#' main="Original Data", ylab="", xlab="")
#' matplot(x=U_true_mat[,1:5], y=Y_reconst_mat[,1:5], col=gray(.5), 
#' type="l", main="Kraus (2015)", ylab="", xlab="")
#' par(mfrow=c(1,1))
reconstructKraus <- function(X_mat, 
                             alpha        = NULL,
                             reconst_fcts = NULL){
  ##
  mean_vec      <- meanKraus(X_mat)
  cov_mat       <- covKraus(X_mat)
  n             <- ncol(X_mat)
  if(is.null(reconst_fcts)){
    reconst_fcts <- 1:n
  }
  X_reconst_mat <- X_mat[, reconst_fcts, drop=FALSE]
  ##
  NonNA_fcts    <- apply(X_mat,2,function(x)!any(is.na(x)))
  X_Compl_mat   <- X_mat[,NonNA_fcts]
  alpha_vec     <- rep(NA, length(reconst_fcts)) 
  df_vec        <- rep(NA, length(reconst_fcts)) 
  ##
  for(i in 1:length(reconst_fcts)){
    X_tmp      <- X_mat[,reconst_fcts[i]]
    ##
    M_bool_vec <- is.na(X_tmp)
    O_bool_vec <- !M_bool_vec
    ##
    if(is.null(alpha)){
      alpha_vec[i] <- stats::optimize(f = function(alpha){gcvKraus(cov_mat     = cov_mat, 
                                                                   mean_vec    = mean_vec, 
                                                                   X_Compl_mat = X_Compl_mat, 
                                                                   M_bool_vec  = M_bool_vec, 
                                                                   alpha       = alpha)},
                                      interval = c(.Machine$double.eps, sum(diag(cov_mat))*n), maximum = FALSE)$minimum
    }else{
      alpha_vec[i] <- alpha
    }
    ##
    result_tmp <- reconstKraus_fun(cov_mat    = cov_mat, 
                                   X_cent_vec = c(X_tmp - mean_vec), 
                                   alpha      = alpha_vec[i])
    ##
    X_reconst_mat[,i]  <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
    df_vec[i]          <- result_tmp[['df']]
  }
  return(list("X_reconst_mat"    = X_reconst_mat,
              "alpha_median"     = stats::median(alpha_vec), 
              "df_median"        = stats::median(df_vec)
              ))
}


## -------------------------------------------------------------------------
meanKraus <- function(X_mat){
  rowMeans(X_mat, na.rm = TRUE)
}


## -------------------------------------------------------------------------
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


## -------------------------------------------------------------------------
reconstKraus_fun <- function(cov_mat, X_cent_vec, alpha=1e-4){
  ##
  M_bool_vec       <- is.na(X_cent_vec)
  O_bool_vec       <- !M_bool_vec
  p                <- nrow(cov_mat)
  ##
  covMO_mat        <- cov_mat[M_bool_vec, O_bool_vec] 
  covOO_mat        <- cov_mat[O_bool_vec, O_bool_vec]
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


## -------------------------------------------------------------------------
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


