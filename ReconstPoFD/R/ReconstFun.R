#' Reconstruct partially observed functions
#'
#' This function allows you to reconstruct the missing parts of a function given the observed parts.
#' @param Ly           List of Y-values. The ith (i=1,...,n) list-element contains \eqn{Y_{i1},\dots,Y_{im}}{Y_{i1},...,Y_{im}}
#' @param Lu           List of U-values. The ith (i=1,...,n) list-element contains \eqn{U_{i1},\dots,U_{im}}{U_{i1},...,U_{im}}
#' @param K            Truncation parameter. If K=NULL (default), K is determined using an AIC-type criterion.
#' @param K_max        Maximum K (used in the AIC-type criterion)
#' @param method       If method=PS_TRUE:  Pre-smoothing of the 'observed' part (Reconstruction operator: \eqn{L^*}{L*}). 
#'                     If method=PS_FALSE: FPCA-estimation of the 'observed' part (Reconstruction operator: \eqn{L}{L}).
#'                     If method=CEScores: FPCA-estimation of the 'observed' part with CEScores from the fdapace package.
#' @param reconst_fcts A vector specifying the list elements in Ly which need to be reconstructed. Default (reconst_fcts=NULL) will reconstruct all functions.
#' @param nRegGrid     Number of grid-points used for the equidistant 'workGrid'; needed for the fdapace::FPCA() function among others.
#' @param messages     Printing messages? (default: messages=FALSE)
#' @export reconstruct
#' @examples  
#' a <- 0; b <- 1; n <- 50
#' SimDat   <- simuldataKraus(n = n, a = a, b = b)
#' ## 
#' Y_list   <- SimDat[['Y_list']]; Y_mat <- SimDat[['Y_mat']]
#' U_list   <- SimDat[['U_list']]; U_mat <- SimDat[['U_mat']]
#' ##
#' reconst_result_1 <- reconstruct(Ly = Y_list, Lu = U_list, method = "PS_TRUE",
#' reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' reconst_result_2 <- reconstruct(Ly = Y_list, Lu = U_list, method = "PS_FALSE", 
#' reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' reconst_result_3 <- reconstruct(Ly = Y_list, Lu = U_list, method = "CEScores",
#' reconst_fcts = 1:3)
#' Y_reconst_mat_3  <- matrix(unlist(reconst_result_3[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_3  <- matrix(unlist(reconst_result_3[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(2,2))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="PS_TRUE", ylab="", xlab="", xlim=c(a,b))
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="PS_FALSE", ylab="", xlab="", xlim=c(a,b))
#' matplot(x=U_reconst_mat_3, y=Y_reconst_mat_3, col=gray(.5), 
#' type="l", main="CEScores", ylab="", xlab="", xlim=c(a,b))
#' par(mfrow=c(1,1))
reconstruct <- function(Ly,
                        Lu,
                        K            = NULL,
                        K_max        = 4,
                        method       = c("PS_TRUE", "PS_FALSE", "CEScores")[3],
                        reconst_fcts = NULL,
                        nRegGrid     = 51,
                        messages     = FALSE)
{
  ##
  n  <- length(Ly)
  if(is.null(reconst_fcts)){
    reconst_fcts <- 1:n
  }
  ##
  ## Estimate Mean and Covariance 
  fdapace_obj <- fdapace::FPCA(Ly    = Ly, 
                               Lt    = Lu, 
                               optns = list(
                                 "dataType"       = "Sparse", 
                                 "kernel"         = "gauss",
                                 "methodMuCovEst" = "smooth",
                                 "nRegGrid"       = nRegGrid
                                 ))
  ## Regular grid
  workGrid        <- fdapace_obj$workGrid
  ## Covariance
  cov_est_mat     <- fdapace_obj$smoothedCov
  ## Mean
  mu_est_fun      <- stats::splinefun(y= fdapace_obj$mu, x=workGrid)
  
  ## Centering the data ###############################################
  Ly_cent <- vector("list", n)
  for(i in 1:n){
    Ly_cent[[i]] <- Ly[[i]] - mu_est_fun(Lu[[i]])
  }
  ##
  ## Aligning Y_cent and U according to 'workGrid' ########################
  Y_cent_mat  <- matrix(data = NA, nrow = nRegGrid, ncol = n)
  U_mat       <- matrix(data = NA, nrow = nRegGrid, ncol = n)
  for(i in 1:n){
    Y_cent_tmp  <- c(stats::na.omit(Ly_cent[[i]]))
    U_tmp       <- c(stats::na.omit(Lu[[i]]))
    for(j in 1:length(U_tmp)){
      loc               <- order(abs(workGrid - U_tmp[j]))[1]
      Y_cent_mat[loc,i] <- Y_cent_tmp[j]
    }
    na_loc            <- is.na(Y_cent_mat[,i])
    U_mat[!na_loc,i]  <- workGrid[!na_loc]       
  }
  ## From matrix to list:
  Ly_align_cent <- split(Y_cent_mat, rep(1:n, each = nRegGrid))
  Lu_align      <- split(U_mat,      rep(1:n, each = nRegGrid))
  ##
  ## K AIC
  if(is.null(K) ){#& method!="CEScores"
    K_AIC <- K_aic_fun(Ly_cent     = Ly_align_cent,
                       Lu          = Lu_align,
                       cov_la_mat  = cov_est_mat,
                       workGrid    = workGrid,
                       K_max       = K_max,
                       messages    = messages)
    K <- K_AIC
  }
  # if(is.null(K) & method=="CEScores"){
  #   K <- length(fdapace_obj$lambda)
  # }
  ## ##################################################################
  ## Re-Fitted Covariance:
  if(method!="CEScores"){
    e_list        <- eigen(cov_est_mat, symmetric = TRUE)
    positiveInd   <- e_list[['values']] >= 0
    eval_vec      <- e_list[['values']][positiveInd]
    evec_mat      <- e_list[['vectors']][,positiveInd, drop=FALSE]
    if(K>1){
      cov_est_refit_mat <- evec_mat[,1:K,drop=FALSE] %*% diag(eval_vec[1:K]) %*% t(evec_mat[,1:K,drop=FALSE])
    }
    if(K==1){
      cov_est_refit_mat <- evec_mat[, 1,drop=FALSE]  %*% t(evec_mat[,1,drop=FALSE]) * eval_vec[1]
    }
  }
  ## ##################################################################
  ## Reconstructing all functions
  ## As list, since this facilitates a future generalization to 'random m'
  Y_reconst_list  <- vector("list", length(reconst_fcts))
  U_reconst_list  <- vector("list", length(reconst_fcts))
  ##
  for(i in 1:length(reconst_fcts)){ # i <- 1
    if(method!="CEScores"){
      tmp  <- reconst_fun(cov_la_mat  = cov_est_refit_mat, 
                          workGrid    = workGrid, 
                          Y_cent_sm_i = c(stats::na.omit(Y_cent_mat[,reconst_fcts[i]])), 
                          U_sm_i      = c(stats::na.omit(U_mat[,reconst_fcts[i]])), 
                          K           = K, 
                          pre_smooth  = ifelse(method=="PS_TRUE", TRUE, FALSE),
                          messages    = messages)
    }
    ##
    if(method=="CEScores"){
      tmp  <- reconst_use_CEScores_fun(Y_cent_mat  = Y_cent_mat,
                                       U_mat       = U_mat,
                                       reconst_fct = reconst_fcts[i],
                                       fdapace_obj = fdapace_obj,
                                       K           = K, 
                                       pre_smooth  = FALSE,
                                       messages    = messages)
    }
    ##
    x_tmp      <- tmp[['x_reconst']]
    y_cent_tmp <- tmp[['y_reconst']]
    y_tmp      <- y_cent_tmp + mu_est_fun(x_tmp)
    ##
    Y_reconst_list[[i]]  <- y_tmp
    U_reconst_list[[i]]  <- x_tmp
  }
  ##
  return(list(
    "Y_reconst_list"  = Y_reconst_list,
    "U_reconst_list"  = U_reconst_list,
    "K"               = K,
    "fdapace_obj"     = fdapace_obj
    ))
}

## ########################################################################
## ########################################################################
reconst_fun <- function(
  cov_la_mat,     
  workGrid,    
  Y_cent_sm_i,    
  U_sm_i,         
  K,              
  pre_smooth      = FALSE,
  messages        = FALSE  
){
  ##
  ## Extracting the [A_i,B_i]^2 part from the large cov-matrix:
  sm_compl_gridloc        <- workGrid>=min(U_sm_i, na.rm = TRUE) & workGrid<=max(U_sm_i, na.rm = TRUE)
  cov_sm_compl_mat        <- cov_la_mat[sm_compl_gridloc, sm_compl_gridloc]
  grid_sm_compl_vec       <- workGrid[sm_compl_gridloc]
  ##
  if(pre_smooth==TRUE){
    smooth.fit              <- stats::smooth.spline(y=Y_cent_sm_i, x=U_sm_i)
    Y_cent_sm_compl_fit_i   <- stats::predict(smooth.fit,grid_sm_compl_vec)$y
  }
  ##
  ## Compute 'small' eigenvalues and eigenfunctions
  e_sm_compl         <- eigen(cov_sm_compl_mat, symmetric = TRUE)
  positiveInd        <- e_sm_compl[['values']] >= 0
  eval_sm_compl      <- e_sm_compl[['values']][positiveInd]
  evec_sm_compl      <- e_sm_compl[['vectors']][,positiveInd, drop=FALSE]
  ##
  ## ##########################################
  ## Check and set K
  ## ##########################################
  if(length(eval_sm_compl)<K){
    if(messages){
      warning("Less non-negative eigenvalues than K. Therefore, K is set to the number of non-negative eigenvalues.")
    }
    K  <- length(eval_sm_compl)
  }
  ## Standardize direction of eigenfunctions
  for(k in 1:K){
    evec_sm_compl[,k] <- evec_sm_compl[,k]*sign(stats::cov(evec_sm_compl[,k], 1:length(evec_sm_compl[,k])))
  }
  ## #################################################
  ## 'Extrapolated/Reconstructive' eigenfunctions
  ela_reconst    <- matrix(NA, nrow=length(workGrid), ncol=K)
  for(k in 1:K){
    ela_reconst[,k] <- c(evec_sm_compl[,k] %*% cov_la_mat[sm_compl_gridloc,,drop=FALSE])
  }
  ## #################################################
  ## Estimating small PC-scores ######################
  ## Evaluate 'evec_sm_compl' at 'U_sm_i'
  evec_sm_compl_at_sm_i        <- matrix(NA, length(U_sm_i), K)
  for(k in 1:K){
    evec_sm_compl_at_sm_i[,k]  <- stats::spline(y=evec_sm_compl[,k], x=grid_sm_compl_vec, xout = U_sm_i)$y
  }
  ## 'Small' PC-scores (OLS-approach)
  xi_sm_i         <- unname(c(stats::lm(c(stats::na.omit(Y_cent_sm_i)) ~ -1 + evec_sm_compl_at_sm_i[,1:K,drop=FALSE])$coefficients))
  ## Refitting (only of effect in the more noisy first run)
  Y_cent_i_fit    <- c(evec_sm_compl[,1:K,drop=FALSE] %*% xi_sm_i)
  xi_sm_i_refit   <- try(unname(c(stats::lm(Y_cent_i_fit ~ -1 + evec_sm_compl[,1:K,drop=FALSE])$coefficients)))
  ## Use the refitted version only if no error was produced:
  if(!is.error(xi_sm_i_refit)){xi_sm_i <- xi_sm_i_refit} 
  ## #################################################
  ##
  ## Recovering: ###########################
  reconst_ls_vec      <- rep(0, length(workGrid))
  for(k in 1:K){
    reconst_ls_vec    <- c(reconst_ls_vec + c((xi_sm_i[k]/eval_sm_compl[k]) * ela_reconst[,k]))
  }
  y_reconst_vec <- reconst_ls_vec[!is.na(reconst_ls_vec)]
  x_reconst_vec <- workGrid[!is.na(reconst_ls_vec)]
  if(pre_smooth==TRUE){
    ## Aliment with observed part 'Y_cent_sm_compl_fit_i':
    sm.gr.loc                       <- c(1:length(workGrid))[sm_compl_gridloc]
    ## lower marginal point
    y_reconst_vec[1:min(sm.gr.loc)] <- y_reconst_vec[1:min(sm.gr.loc)] + Y_cent_sm_compl_fit_i[1] - y_reconst_vec[min(sm.gr.loc)]
    ## upper marginal point
    y_reconst_vec[max(sm.gr.loc):length(workGrid)] <- y_reconst_vec[max(sm.gr.loc):length(workGrid)] + Y_cent_sm_compl_fit_i[length(Y_cent_sm_compl_fit_i)] - y_reconst_vec[max(sm.gr.loc)]
    ## estimated ('observed') part
    y_reconst_vec[sm.gr.loc]                       <- Y_cent_sm_compl_fit_i
  }
  ## ######################
  return(list("y_reconst"  = c(y_reconst_vec),
              "x_reconst"  = c(x_reconst_vec)))
  ## ######################
}

## ###########################################################
## ###########################################################
reconst_use_CEScores_fun <- function(
  Y_cent_mat,
  U_mat, 
  reconst_fct,
  fdapace_obj,
  K,
  pre_smooth = FALSE,
  messages = FALSE  
){
  ##
  ## Regular grid
  workGrid        <- fdapace_obj$workGrid
  ## Covariance
  cov_la_mat      <- fdapace_obj$smoothedCov
  ## Mean
  mu_est_vec      <- fdapace_obj$mu
  ##
  U_sm_min_i              <- min(U_mat[,reconst_fct], na.rm = TRUE)
  U_sm_max_i              <- max(U_mat[,reconst_fct], na.rm = TRUE)
  sm_compl_gridloc        <- workGrid>=U_sm_min_i & workGrid<=U_sm_max_i
  #cov_sm_compl_mat        <- cov_la_mat[sm_compl_gridloc, sm_compl_gridloc]
  grid_sm_compl_vec       <- workGrid[sm_compl_gridloc]
  ##
  if(pre_smooth==TRUE){
    smooth.fit              <- stats::smooth.spline(y=Y_cent_mat[,reconst_fct], x=U_mat[,reconst_fct])
    Y_cent_sm_compl_fit_i   <- stats::predict(smooth.fit,grid_sm_compl_vec)$y
  }
  ##
  ## Compute 'small' eigenvalues and eigenfunctions
  # e_sm_compl         <- eigen(cov_sm_compl_mat, symmetric = TRUE)
  # positiveInd        <- e_sm_compl[['values']] >= 0
  # eval_sm_compl      <- e_sm_compl[['values']][positiveInd]
  # evec_sm_compl      <- e_sm_compl[['vectors']][,positiveInd, drop=FALSE]
  ##
  ## PACE on small:
  Ly_sm_i <- vector("list", ncol(Y_cent_mat))
  Lu_sm_i <- vector("list", ncol(Y_cent_mat))
  for(i in seq_len(ncol(Y_cent_mat))){
    Ly_sm_i[[i]] <- c(stats::na.omit(Y_cent_mat[,i][U_mat[,i]>=U_sm_min_i & U_mat[,i]<=U_sm_max_i]))
    Lu_sm_i[[i]] <- c(stats::na.omit(     U_mat[,i][U_mat[,i]>=U_sm_min_i & U_mat[,i]<=U_sm_max_i]))
  }
  ## 
  #grid_sm_compl_vec <- seq(from=range(c(unlist(Lu_sm_i)))[1], to=range(c(unlist(Lu_sm_i)))[2], length.out = length(grid_sm_compl_vec))
  
  fdapace_sm_obj <- fdapace::FPCA(Ly    = Ly_sm_i, 
                                  Lt    = Lu_sm_i, 
                                  optns = list(
                                    "dataType"       = "Sparse", 
                                    "kernel"         = "gauss",
                                    "methodMuCovEst" = "smooth",
                                    #"userCov"=list("t"=grid_sm_compl_vec, "cov"=cov_sm_compl_mat),
                                    #"userMu" =list("t"=grid_sm_compl_vec, "mu" =rep(0,length(grid_sm_compl_vec))),
                                    #"outPercent"=c(0,1)
                                    "methodSelectK" = K#length(fdapace_obj$lambda)
                                  ))
  ##
  K              <- length(fdapace_sm_obj$lambda)#min(length(fdapace_sm_obj$lambda), ncol(evec_sm_compl))
  CEScores_vec_i <- c(fdapace_sm_obj$xiEst[reconst_fct, seq_len(K), drop=FALSE])
  efcts_sm_pace  <- matrix(0, nrow=length(grid_sm_compl_vec), ncol=K)
  ##
  for(k in seq_len(K)){
    efcts_sm_pace[,k]  <- stats::spline(x = fdapace_sm_obj$workGrid, y=fdapace_sm_obj$phi[, k, drop=FALSE], xout = grid_sm_compl_vec)$y
  }
  ## matplot(cbind(evec_sm_compl[,1], efcts_sm_pace[,1]))
  ##
  ## #################################################
  ## 'Extrapolated/Reconstructive' eigenfunctions
  ela_reconst    <- matrix(NA, nrow=length(workGrid), ncol=K)
  for(k in seq_len(K)){
        ela_reconst[,k] <- apply(X      = cov_la_mat[sm_compl_gridloc,,drop=FALSE],
                                 MARGIN = 2,
                                 FUN    = function(x){pracma::trapz(x=grid_sm_compl_vec, y = efcts_sm_pace[,k] * x)})
    ##
    ## ela_reconst[,k] <- c(evec_sm_compl[,k] %*% cov_la_mat[sm_compl_gridloc,,drop=FALSE])
    # re-scale:
    #sign_tmp        <- sign(ela_reconst[sm_compl_gridloc,k] %*% efcts_sm_pace[,k])
    scale_tmp       <- sqrt(sum(efcts_sm_pace[,k]^2)) / sqrt(sum(ela_reconst[sm_compl_gridloc,k]^2)) #* sign_tmp
    ela_reconst[,k] <- ela_reconst[,k] * scale_tmp
  }
  ## k <- 2; matplot(cbind(ela_reconst[sm_compl_gridloc,k], efcts_sm_pace[,k]))
  ##
  Y_cent_fit_i <- numeric(length = length(workGrid))
  for(k in seq_len(K)){
    Y_cent_fit_i <- Y_cent_fit_i + CEScores_vec_i[k] * ela_reconst[,k]
  }
  ## ######################
  return(list("y_reconst"  = c(Y_cent_fit_i),
              "x_reconst"  = c(workGrid)))
  ## ######################
}

#   slct_vec        <- match(U_sm_i, workGrid)
#   ##
#   eigen_la_obj    <- eigen(cov_la_mat, symmetric = TRUE)
#   evals_la_vec    <- eigen_la_obj[['values']]
#   evec_la_mat     <- eigen_la_obj[['vectors']]
#   non_neg_la_ev   <- evals_la_vec>0
#   evals_la_vec    <- evals_la_vec[non_neg_sm_ev]
#   evec_la_mat     <- evec_la_mat[,non_neg_sm_ev]
#   K_la_MAX        <- length(evals_la_vec)
#   ##
#   cov_sm_mat      <- cov_la_mat[slct_vec, slct_vec]
#   eigen_sm_obj    <- eigen(cov_sm_mat, symmetric = TRUE)
#   evals_sm_vec    <- eigen_sm_obj[['values']]
#   evec_sm_mat     <- eigen_sm_obj[['vectors']]
#   non_neg_sm_ev   <- evals_sm_vec>0
#   evals_sm_vec    <- evals_sm_vec[non_neg_sm_ev]
#   evec_sm_mat     <- evec_sm_mat[,non_neg_sm_ev]
#   K_sm_MAX        <- length(evals_sm_vec)
#   ##
#   K               <- min(K_la_MAX, K_sm_MAX, K_AIC)
#   SIGMA_mat       <- cov_sm_mat + diag(sigma2, nrow = nrow(cov_sm_mat))
#   SIGMA_INV_mat   <- solve(SIGMA_mat)
#   ##
# 
#   ##
#   Y_fit_i <- mu_est_vec + Y_cent_fit_i
#   cbind(Y_fit_i, fitted(fdapace_obj)[1,])
#   matplot(y=cbind(Y_fit_i, fitted(fdapace_obj)[1,]), x=workGrid)
#   lines(y=SimDat$Y_true_list[[1]], x=SimDat$U_true_list[[1]])
#   
#   
#   
#   
#   
#   K          <- length(fdapace_obj$lambda)
#   workGrid   <- fdapace_obj$workGrid
#   cov_la_mat <- fdapace_obj$fittedCov
#   ##
#   ## Extracting the [A_i,B_i]^2 part from the large cov-matrix:
#   sm_compl_gridloc        <- workGrid>=min(U_sm_i, na.rm = TRUE) & workGrid<=max(U_sm_i, na.rm = TRUE)
#   cov_sm_compl_mat        <- cov_la_mat[sm_compl_gridloc, sm_compl_gridloc]
#   grid_sm_compl_vec       <- workGrid[sm_compl_gridloc]
#   ##
#   if(pre_smooth==TRUE){
#     smooth.fit              <- stats::smooth.spline(y=Y_cent_sm_i, x=U_sm_i)
#     Y_cent_sm_compl_fit_i   <- stats::predict(smooth.fit,grid_sm_compl_vec)$y
#   }
#   ##
#   ## Compute 'small' eigenvalues and eigenfunctions
#   # e_sm_compl         <- eigen(cov_sm_compl_mat, symmetric = TRUE)
#   # positiveInd        <- e_sm_compl[['values']] >= 0
#   # eval_sm_compl      <- e_sm_compl[['values']][positiveInd]
#   # evec_sm_compl      <- e_sm_compl[['vectors']][,positiveInd, drop=FALSE]
#   evec_sm_compl <- fdapace_obj$phi
#   ##
#   ## Standardize direction and L2norm of 'small' eigenfunctions
#   for(k in 1:K){
#     evec_sm_compl[,k] <- evec_sm_compl[,k]*sign(stats::cov(evec_sm_compl[,k], 1:length(evec_sm_compl[,k])))
#     scale             <- sqrt(pracma::trapz(x=grid_sm_compl_vec, y=evec_sm_compl[,k]^2))
#     evec_sm_compl[,k] <- evec_sm_compl[,k]/scale
#   }
#   ## #################################################
#   ## 'Extrapolated/Reconstructive' eigenfunctions
#   ela_reconst       <- matrix(NA, nrow=length(workGrid), ncol=K)
#   for(k in 1:K){
#     # ela_reconst[,k] <- c(evec_sm_compl[,k] %*% cov_la_mat[sm_compl_gridloc,,drop=FALSE])
#     ela_reconst[,k] <- apply(X      = cov_la_mat[sm_compl_gridloc,,drop=FALSE], 
#                              MARGIN = 2, 
#                              FUN    = function(x){pracma::trapz(x=grid_sm_compl_vec, evec_sm_compl[,k] * x)})
#   }
#   ## #################################################
#   xi_sm_i <- fdapace_obj$xiEst
#   ## #################################################
#   ##
#   ## Recovering: ###########################
#   reconst_ls_vec  <- rep(0, length(workGrid))
#   for(k in 1:K){
#     reconst_ls_vec    <- c(reconst_ls_vec + c((xi_sm_i[k]/fdapace_obj$lambda[k]) * ela_reconst[,k]))
#   }
#   y_reconst_vec <- reconst_ls_vec[!is.na(reconst_ls_vec)]
#   x_reconst_vec <- workGrid[!is.na(reconst_ls_vec)]
#   if(pre_smooth==TRUE){
#     ## Aliment with observed part 'Y_cent_sm_compl_fit_i':
#     sm.gr.loc                       <- c(1:length(workGrid))[sm_compl_gridloc]
#     ## lower marginal point
#     y_reconst_vec[1:min(sm.gr.loc)] <- y_reconst_vec[1:min(sm.gr.loc)] + Y_cent_sm_compl_fit_i[1] - y_reconst_vec[min(sm.gr.loc)]
#     ## upper marginal point
#     y_reconst_vec[max(sm.gr.loc):length(workGrid)] <- y_reconst_vec[max(sm.gr.loc):length(workGrid)] + Y_cent_sm_compl_fit_i[length(Y_cent_sm_compl_fit_i)] - y_reconst_vec[max(sm.gr.loc)]
#     ## estimated ('observed') part
#     y_reconst_vec[sm.gr.loc]                       <- Y_cent_sm_compl_fit_i
#   }
#   ## ######################
#   return(list("y_reconst"  = c(y_reconst_vec),
#               "x_reconst"  = c(x_reconst_vec)))
#   ## ######################
# }


## ###########################################################
## ###########################################################
K_aic_fun <- function(Ly_cent, 
                      Lu,
                      cov_la_mat,
                      workGrid,
                      K_max      = 4,
                      pre_smooth = FALSE,
                      messages   = FALSE)
{
  n <- length(Ly_cent)
  ## #########################################################
  ## Nonparametric variance estimation 
  ## Gasser, Stroka, Jennen-Steinmetz (1986, Biometrika)
  ## #########################################################
  sig2_GSJ_eps_vec <- NULL
  for(i in 1:n){
    y.vec <- c(stats::na.omit(Ly_cent[[i]]))
    x.vec <- c(stats::na.omit(Lu[[i]]))
    ##
    a.seq      <- (diff(x.vec, lag=1)[-1]/diff(x.vec, lag=2))
    b.seq      <- (diff(x.vec, lag=1)[-length(diff(x.vec, lag=1))]/diff(x.vec, lag=2))
    c.sq       <- 1/(a.seq^2+b.seq^2+1)
    ##
    pseudo.eps <- a.seq * y.vec[-c( length(y.vec)-1,  length(y.vec))] + b.seq * y.vec[-c(1,2)] - y.vec[-c(1,length(y.vec))]
    sig2.GSJ   <- mean(c(pseudo.eps^2*c.sq))
    sig2_GSJ_eps_vec <- c(sig2_GSJ_eps_vec, sig2.GSJ)
  }
  sig2_GSJ_eps  <- mean(sig2_GSJ_eps_vec)
  
  
  ## ############################
  ## Selecting K via AIC 
  ## ############################
  AIC.vec <- rep(NA,K_max)
  ##
  for(K_iter in 1:K_max){
    RSS.vec <- rep(NA,n)
    L.vec   <- rep(NA,n)
    for(i in 1:n){ 
      Y_cent_sm_i <- c(stats::na.omit(Ly_cent[[i]]))
      U_sm_i      <- c(stats::na.omit(Lu[[i]]))
      lo.half     <- 1:floor(length(U_sm_i)/2)
      up.half     <- (floor(length(U_sm_i)/2)+1):length(U_sm_i)
      ##
      List_obj    <- reconst_fun(cov_la_mat     = cov_la_mat, 
                                 workGrid       = workGrid, 
                                 Y_cent_sm_i    = Y_cent_sm_i[lo.half], 
                                 U_sm_i         = U_sm_i[lo.half], 
                                 K              = K_iter,
                                 pre_smooth     = pre_smooth, 
                                 messages       = messages)
      ##
      y_fit      <- List_obj[['y_reconst']] #+ mu_est_fun(List_obj[['x_reconst']])
      ## RSS.vec[i] <- sum((y_fit[!is.na(Lu[[i]])][up.half] - (Y_cent_sm_i[up.half] + mu_norm_est_fun(U_sm_i)[up.half]))^2)
      RSS.vec[i] <- sum((y_fit[!is.na(Lu[[i]])][up.half] - Y_cent_sm_i[up.half] )^2)
      L.vec[i]   <- -(n*log(2*pi)/2)-(n*log(sig2_GSJ_eps)/2)-(RSS.vec[i]/(2*sig2_GSJ_eps))
    }
    AIC.vec[K_iter] <- -sum(L.vec) + K_iter 
  }
  K_AIC <- which.min(AIC.vec)
  return(K_AIC)
}

#' Iterative Reconstruction Algoritm
#'
#' This function iteratively applies the function reconst_fun() in order to reconstruct the missing parts of a function given the observed parts. 
#' The iterative procedure allows to reconstruct functions when their covariance function cannot be estimated over the total domain. However, the covariance function must be estimated over a sufficiently large part of the domain.  
#' 
#' @param cov_la_mat  Discretized covariance function over \eqn{[a,b]\times[a,b]}{[a,b]x[a,b]}
#' @param workGrid    Equidistant discretization grid in \eqn{[a,b]}{[a,b]}
#' @param Y_cent_sm_i Centered function values of the ith function: \eqn{Y_{ij}-\hat\mu(U_{ij}), j=1,\dots,m}{Y_{ij}-\hat\mu(U_{ij}), j=1,...,m},
#' @param U_sm_i      Discretization points of the ith function: \eqn{U_{i1},\dots,U_{im}}{U_{i1},...,U_{im}}
#' @param K           Truncation parameter
#' @param fraction    A value between 0 and 1. (See details)
#' @param max_rep     Maximal number of iterations (default max_rep=5)
#' @param pre_smooth  If pre_smooth==TRUE:  Pre-smoothing of the 'observed' part.  (Reconstruction operator: \eqn{L^*}{L*}). If pre_smooth==FALSE (default): FPCA-estimation of the 'observed' part (Reconstruction operator: \eqn{L}{L})
#' @param messages    Printing messages if the algorithm stops (default messages=FALSE)
#' @details 
#' Idea of the procedure: In each iteration the observed function is devided in a upper and lower fragment which are used to reconstruct further upper and lower missing parts. 
#' 
#' The lengths of the upper and lower fragments are determined by the argument 'fraction'. 
#' Large values of 'fraction' lead to large (i.e., more informative) fragments which can improve the reconstructions, but will not allow to reconstruct large missing parts. 
#' Small values of 'fraction' lead to small (i.e., less informative) fragments which can worsen the reconstructions, but will allow to reconstruct large missing parts. 
#' @export iter_reconst_fun
iter_reconst_fun <- function(cov_la_mat, 
                             workGrid, 
                             Y_cent_sm_i, 
                             U_sm_i, 
                             K, 
                             fraction   = 0.15,
                             max_rep    = 5,
                             pre_smooth = FALSE,
                             messages   = FALSE
){
  
  ## extract data:
  grid.len <- length(workGrid)
  
  ## Checks #############################################################
  if(grid.len * fraction < K+1){stop("'fraction' too small.")}
  
  
  ## ####################################################################
  ## First Run
  pred.tmp     <- reconst_fun(cov_la_mat      = cov_la_mat, 
                                           workGrid        = workGrid, 
                                           Y_cent_sm_i     = Y_cent_sm_i, 
                                           U_sm_i          = U_sm_i, 
                                           K               = K)
  ##
  U_la_i      <- pred.tmp[['x_reconst']]
  Y_cent_la_i <- pred.tmp[['y_reconst']]
  ##  check-plot:
  # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(workGrid), type="l"); abline(v=range(workGrid))
  
  ## ####################################################################
  ## Do we need to predict a further missing *upper* fragment?
  ## ####################################################################
  if(max(U_la_i) !=  max(workGrid)){
    ## Find the first domain value at which the there is no missing cov-value anymore with 
    ## respect to the largest domain values, i.e., 'cov_la_mat[,grid.len]' (could be done outside of the loop!)
    min.loc         <- which.min(cov_la_mat[,grid.len])
    ##
    ## Is the predicted function already (sufficiently far) beyond 'workGrid[min.loc]'?
    ## ('Sufficiently far' >= 'grid.len * fraction' grid points)
    check <- max(U_la_i) > (workGrid[min.loc]) & 
      length(workGrid[workGrid >= workGrid[min.loc] & workGrid <= max(U_la_i)]) >= (grid.len * fraction)
    stop.pred <- FALSE  # initialization
    if(check == FALSE){
      go    <- TRUE
      count <- 1
      ## Predict further parts of the upper missing fragment:
      while(go){ # print("In while loop for upper pred.")
        ## Select new small intervall: 
        ## [(max(U_la_i) - c(max(U_la_i)-min(U_la_i))*fraction),  max(U_la_i)]
        U_slct_sm2_vec  <- workGrid <= max(U_la_i) &  workGrid >= (max(U_la_i) - c(max(U_la_i)-min(U_la_i))*fraction)
        U_sm2_i         <- workGrid[U_slct_sm2_vec]
        ##
        ## Is the new small interval large enough?
        if(length(U_sm2_i) >= grid.len * fraction){
          ## If yes, start the prediction:
          Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
          Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
          ##
          pred2.tmp    <- reconst_fun(cov_la_mat      = cov_la_mat, 
                                                   workGrid        = workGrid, 
                                                   Y_cent_sm_i     = Y_cent_sm2_i, 
                                                   U_sm_i          = U_sm2_i, 
                                                   K               = K)
          U_la2_i      <- pred2.tmp[['x_reconst']]
          Y_cent_la2_i <- pred2.tmp[['y_reconst']]
          ##
          ## Isolate the newly predicted fragment:  
          U_la2_fragm_i      <- U_la2_i[U_la2_i >= max(U_la_i)]
          Y_cent_la2_fragm_i <- Y_cent_la2_i[U_la2_i >= max(U_la_i)]
          ## check-plot:
          # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(workGrid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
          # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, lty=2)
          # abline(h=c(Y_cent_la_i[length(Y_cent_la_i)], Y_cent_la2_fragm_i[1]), col=c("red", "black"))
          ##
          ## Additive alignment:
          align.tmp          <- Y_cent_la_i[length(Y_cent_la_i)] - Y_cent_la2_fragm_i[1]
          Y_cent_la2_fragm_i <- Y_cent_la2_fragm_i + align.tmp
          ## Final joining:
          U_la_i             <- c(U_la_i,       U_la2_fragm_i[-1])
          Y_cent_la_i        <- c(Y_cent_la_i,  Y_cent_la2_fragm_i[-1])
          ## check-plot:
          # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(workGrid), type="l")
          # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
          
          ## Check again:
          ## Is the predicted function now (sufficiently far) beyond 'workGrid[min.loc]'?
          check <- max(U_la_i) > (workGrid[min.loc]) & 
            length(workGrid[workGrid >= workGrid[min.loc] & workGrid <= max(U_la_i)]) >= (grid.len * fraction)
          ## 
          go    <- ifelse(check, yes = FALSE, no = TRUE)
          count <- count + 1
        }else{
          ## If the small interval is too small:
          if(messages){
            warning(cat("The observed interval is too small.\n", 
                        "Prediction process of the upper fragment is stopped.\n", sep = ""))
          }
          go        <- FALSE
          stop.pred <- TRUE
        }
        if(count == max_rep){
          if(messages){
            warning(cat("Maximal repetitions reached.\n", 
                        "Prediction process of the upper fragment is stopped.\n", sep = ""))
          }
          go        <- FALSE
          stop.pred <- TRUE
        }
      }
    }
    ##
    ## Is the predicted function already (sufficiently far) beyond 'workGrid[min.loc]'?
    ## ('Sufficiently far' >= '('grid.len * fraction' grid points)  
    if(check == TRUE & stop.pred==FALSE){ 
      ## If yes, we can achieve a complete prediction of the missing upper fragments! 
      ##
      U_slct_sm2_vec  <- workGrid >= workGrid[min.loc] & workGrid <= max(U_la_i)
      U_sm2_i         <- workGrid[U_slct_sm2_vec]
      Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
      Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
      ##
      pred2.tmp    <- reconst_fun(cov_la_mat      = cov_la_mat, 
                                  workGrid        = workGrid, 
                                  Y_cent_sm_i     = Y_cent_sm2_i, 
                                  U_sm_i          = U_sm2_i, 
                                  K               = K)
      ##
      U_la2_i      <- pred2.tmp[['x_reconst']]
      Y_cent_la2_i <- pred2.tmp[['y_reconst']]
      ##
      ## Isolate the newly predicted fragment:  
      U_la2_fragm_i      <- U_la2_i[U_la2_i >= max(U_la_i)]
      Y_cent_la2_fragm_i <- Y_cent_la2_i[U_la2_i >= max(U_la_i)]
      ## check-plot:
      # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(workGrid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
      # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, lty=2)
      # abline(h=c(Y_cent_la_i[length(Y_cent_la_i)], Y_cent_la2_fragm_i[1]), col=c("red", "black"))
      ##
      ## Additive alignment:
      align.tmp          <- Y_cent_la_i[length(Y_cent_la_i)] - Y_cent_la2_fragm_i[1]
      Y_cent_la2_fragm_i <- Y_cent_la2_fragm_i + align.tmp
      ## Final joining:
      U_la_i             <- c(U_la_i,       U_la2_fragm_i[-1])
      Y_cent_la_i        <- c(Y_cent_la_i,  Y_cent_la2_fragm_i[-1])
      ## check-plot:
      # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(workGrid), type="l")
      # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
    }
  }
  ## ################################################################################
  ## Do we need to predict a further missing *lower* fragement?
  ## ################################################################################
  if(min(U_la_i) !=  min(workGrid)){
    ## Find the last domain value at which the there is still no missing cov-value with 
    ## respect to the smallest domain value, i.e., 'cov_la_mat[,1]' (could be done outside of the loop!)
    max.loc         <- which.max(cov_la_mat[,1])
    ##
    ## Is the predicted function already (sufficiently far) below 'workGrid[max.loc]'?
    ## ('Sufficiently far' >= 'grid.len * fraction' grid points)
    check <- min(U_la_i) < (workGrid[max.loc]) & 
      length(workGrid[workGrid <= workGrid[max.loc] & workGrid >= min(U_la_i)]) >= (grid.len * fraction)
    stop.pred <- FALSE
    if(check == FALSE){
      go    <- TRUE
      count <- 1
      ## Predict further parts of the lower missing fragment:
      while(go){ # print("In while loop for lower pred.")
        ## Select new small intervall: 
        ## [min(U_la_i), (min(U_la_i) + c(max(U_la_i)-min(U_la_i))*fraction)]
        U_slct_sm2_vec  <- workGrid >= min(U_la_i) &  workGrid <= (min(U_la_i) + c(max(U_la_i)-min(U_la_i))*fraction)
        U_sm2_i         <- workGrid[U_slct_sm2_vec]
        ##
        ## Is the new small interval large enough?
        if(length(U_sm2_i) >= grid.len * fraction){ 
          ## If yes, start the prediction:
          Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
          Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
          ##
          pred2.tmp    <- reconst_fun(cov_la_mat      = cov_la_mat, 
                                                   workGrid        = workGrid, 
                                                   Y_cent_sm_i     = Y_cent_sm2_i, 
                                                   U_sm_i          = U_sm2_i, 
                                                   K               = K)
          U_la2_i      <- pred2.tmp[['x_reconst']]
          Y_cent_la2_i <- pred2.tmp[['y_reconst']]
          ##
          ## Isolate the newly predicted fragment:  
          U_la2_fragm_i      <- U_la2_i[U_la2_i <= min(U_la_i)]
          Y_cent_la2_fragm_i <- Y_cent_la2_i[U_la2_i <= min(U_la_i)]
          ## check-plot:
          # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(workGrid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
          # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, lty=2)
          # abline(h=c(Y_cent_la_i[1], Y_cent_la2_fragm_i[length(Y_cent_la2_fragm_i)]), col=c("red", "black"))
          ##
          ## Additive alignment:
          align.tmp          <- Y_cent_la_i[1] - Y_cent_la2_fragm_i[length(Y_cent_la2_fragm_i)]
          Y_cent_la2_fragm_i <- Y_cent_la2_fragm_i + align.tmp
          ## Final joining:
          U_la_i             <- c(U_la2_fragm_i[-length(U_la2_fragm_i)],                U_la_i)
          Y_cent_la_i        <- c(Y_cent_la2_fragm_i[-length(Y_cent_la2_fragm_i)], Y_cent_la_i)
          ## check-plot:
          # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(workGrid), type="l")
          # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
          
          ## Check again:
          ## Is the predicted function now (sufficiently far) below 'workGrid[min.loc]'?
          check <- min(U_la_i) < (workGrid[max.loc]) & 
            length(workGrid[workGrid <= workGrid[max.loc] & workGrid >= min(U_la_i)]) >= (grid.len * fraction)
          ## 
          go    <- ifelse(check, yes = FALSE, no = TRUE)
          count <- count + 1
        }else{
          ## If the small interval is too small:
          if(messages){
            warning(cat("The small interval is too small for t=",t,".\n", 
                        "Prediction process of the lower fragment is stopped.\n", sep = ""))
          }
          go        <- FALSE
          stop.pred <- TRUE
        }
        if(count == max_rep){
          if(messages){
            warning(cat("Maximal repetitions reached for t=",t,".\n", 
                        "Prediction process of the lower fragment is stopped.\n", sep = ""))
          }
          go        <- FALSE
          stop.pred <- TRUE
        }
      }
    }
    ##
    ## Is the predicted function already (sufficiently far) below 'workGrid[min.loc]'?
    ## ('Sufficiently far' >= '('grid.len * fraction' grid points)  
    if(check == TRUE & stop.pred==FALSE){ 
      ## If yes, we can achieve a complete prediction of the missing lower fragment! 
      ##
      U_slct_sm2_vec  <- workGrid <= workGrid[max.loc] & workGrid >= min(U_la_i)
      U_sm2_i         <- workGrid[U_slct_sm2_vec]
      Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
      Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
      ##
      pred2.tmp    <- reconst_fun(cov_la_mat      = cov_la_mat, 
                                               workGrid        = workGrid, 
                                               Y_cent_sm_i     = Y_cent_sm2_i, 
                                               U_sm_i          = U_sm2_i, 
                                               K               = K)
      U_la2_i      <- pred2.tmp[['x_reconst']]
      Y_cent_la2_i <- pred2.tmp[['y_reconst']]
      ##
      ## Isolate the newly predicted fragment:  
      U_la2_fragm_i      <- U_la2_i[U_la2_i <= min(U_la_i)]
      Y_cent_la2_fragm_i <- Y_cent_la2_i[U_la2_i <= min(U_la_i)]
      ## check-plot:
      # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(workGrid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
      # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, lty=2)
      # abline(h=c(Y_cent_la_i[1], Y_cent_la2_fragm_i[length(Y_cent_la2_fragm_i)]), col=c("red", "black"))
      ##
      ## Additive alignment:
      align.tmp          <- Y_cent_la_i[1] - Y_cent_la2_fragm_i[length(Y_cent_la2_fragm_i)]
      Y_cent_la2_fragm_i <- Y_cent_la2_fragm_i + align.tmp
      ## Final joining:
      U_la_i             <- c(U_la2_fragm_i[-length(U_la2_fragm_i)], U_la_i)
      Y_cent_la_i        <- c(Y_cent_la2_fragm_i[-length(Y_cent_la2_fragm_i)], Y_cent_la_i)
      ## check-plot:
      # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(workGrid), type="l")
      # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
    }
  }
  return(list("y_reconst"=Y_cent_la_i,
              "x_reconst"=U_la_i))
}

## ##########################
## Error catching function
## ##########################
## My error checker function
is.error <- function(x){
  bool.result <- inherits(x, "try-error")
  bool.result <- unname(bool.result)
  return(bool.result)
}
## check:
#result   <- try(log("a"), silent=TRUE) # produces and error
#is.error(result)


