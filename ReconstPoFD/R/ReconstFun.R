#' Reconstruction Function
#'
#' This function allows you to reconstruct the missing parts of a function given the observed parts.
#' @param cov_la_mat  Discretized covariance function over [a,b]x[a,b]
#' @param domain_grid Equidistant discretization grid in [a,b]
#' @param Y_cent_sm_i Centered function values of the ith function: Y[{ij}]-hat(mu)(U[{ij}]), j=1,...,m",
#' @param U_sm_i      Discretization points of the ith function: U[{i1}],...,U[{im}]
#' @param K           Truncation parameter
#' @param pre.smooth  If pre.smooth==TRUE:  Pre-smoothing of the 'observed' part.  (Reconstruction operator: L^*). If pre.smooth==FALSE (default): FPCA-estimation of the 'observed' part (Reconstruction operator: L)
#' @export reconst_fun

reconst_fun <- function(
  cov_la_mat,     # Discretized covariance function over [a,b]x[a,b]
  domain_grid,    # Equidistant discretization grid in [a,b]
  Y_cent_sm_i,    # Centered function values of the ith function: "Y_{ij}-\hat\mu(U_{ij}), j=1,...,m",
  U_sm_i,         # Discretization points of the ith function: U_{i1},...,U_{im}
  K,              # Truncation parameter
  pre.smooth=FALSE
  ## pre.smooth==TRUE:  Pre-smoothing of the 'observed' part.  (Reconstruction operator: L^\ast)
  ## pre.smooth==FALSE: FPCA-estimation of the 'observed' part (Reconstruction operator: L)
){

  ## Extracting the [A_i,B_i]^2 part from the large cov-matrix:
  sm_compl_gridloc        <- domain_grid>=min(U_sm_i, na.rm = TRUE) & domain_grid<=max(U_sm_i, na.rm = TRUE)
  cov_sm_compl_mat        <- cov_la_mat[sm_compl_gridloc, sm_compl_gridloc]
  grid_sm_compl_vec       <- domain_grid[sm_compl_gridloc]
  ##
  if(pre.smooth==TRUE){
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
  if(length(eval_sm_compl)<K){
    warning("Less non-negative eigenvalues than K.")
    K  <- length(eval_sm_compl)
  }
  ## Standardize direction of eigenfunctions
  for(k in 1:K){
    evec_sm_compl[,k] <- evec_sm_compl[,k]*sign(stats::cov(evec_sm_compl[,k], 1:length(evec_sm_compl[,k])))
  }
  ## #################################################
  ## 'Extrapolated/Reconstructive' eigenfunctions
  ela_reconst    <- matrix(NA, nrow=length(domain_grid), ncol=K)
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
  xi_sm_i         <- unname(c(stats::lm(c(na.omit(Y_cent_sm_i)) ~ -1 + evec_sm_compl_at_sm_i[,1:K,drop=FALSE])$coefficients))
  ## Refitting (only of effect in the more noisy first run)
  Y_cent_i_fit    <- c(evec_sm_compl[,1:K,drop=FALSE] %*% xi_sm_i)
  xi_sm_i         <- unname(c(stats::lm(Y_cent_i_fit ~ -1 + evec_sm_compl[,1:K,drop=FALSE])$coefficients))
  ##
  ## Recovering: ###########################
  reconst_ls_vec    <- rep(0, length(domain_grid))
  for(k in 1:K){
    reconst_ls_vec    <- c(reconst_ls_vec + c((xi_sm_i[k]/eval_sm_compl[k]) * ela_reconst[,k]))
  }
  y_reconst_vec <- reconst_ls_vec[!is.na(reconst_ls_vec)]
  x_reconst_vec <- domain_grid[!is.na(reconst_ls_vec)]
  if(pre.smooth==TRUE){
    ## Aliment with observed part 'Y_cent_sm_compl_fit_i':
    sm.gr.loc                       <- c(1:length(domain_grid))[sm_compl_gridloc]
    ## lower marginal point
    y_reconst_vec[1:min(sm.gr.loc)] <- y_reconst_vec[1:min(sm.gr.loc)] + Y_cent_sm_compl_fit_i[1] - y_reconst_vec[min(sm.gr.loc)]
    ## upper marginal point
    y_reconst_vec[max(sm.gr.loc):length(domain_grid)] <- y_reconst_vec[max(sm.gr.loc):length(domain_grid)] + Y_cent_sm_compl_fit_i[length(Y_cent_sm_compl_fit_i)] - y_reconst_vec[max(sm.gr.loc)]
    ## estimated ('observed') part
    y_reconst_vec[sm.gr.loc]                          <- Y_cent_sm_compl_fit_i
  }
  ## ######################
  return(list("y_reconst"  = c(y_reconst_vec),
              "x_reconst"  = c(x_reconst_vec)))
  ## ######################
}
