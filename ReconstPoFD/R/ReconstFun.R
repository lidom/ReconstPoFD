#' Reconstruction Function
#'
#' This function allows you to reconstruct the missing parts of a function given the observed parts.
#' @param cov_la_mat  Discretized covariance function over \eqn{[a,b]\times[a,b]}{[a,b]x[a,b]}
#' @param domain_grid Equidistant discretization grid in \eqn{[a,b]}{[a,b]}
#' @param Y_cent_sm_i Centered function values of the ith function: \eqn{Y_{ij}-\hat\mu(U_{ij}), j=1,\dots,m}{Y_{ij}-\hat\mu(U_{ij}), j=1,...,m},
#' @param U_sm_i      Discretization points of the ith function: \eqn{U_{i1},\dots,U_{im}}{U_{i1},...,U_{im}}
#' @param K           Truncation parameter
#' @param pre_smooth  If pre_smooth==TRUE:  Pre-smoothing of the 'observed' part.  (Reconstruction operator: \eqn{L^*}{L*}). If pre_smooth==FALSE (default): FPCA-estimation of the 'observed' part (Reconstruction operator: \eqn{L}{L})
#' @export reconst_fun
reconst_fun <- function(
  cov_la_mat,     
  domain_grid,    
  Y_cent_sm_i,    
  U_sm_i,         
  K,              
  pre_smooth=FALSE
  ## pre_smooth==TRUE:  Pre-smoothing of the 'observed' part.  (Reconstruction operator: L^\ast)
  ## pre_smooth==FALSE: FPCA-estimation of the 'observed' part (Reconstruction operator: L)
){
  
  ## Extracting the [A_i,B_i]^2 part from the large cov-matrix:
  sm_compl_gridloc        <- domain_grid>=min(U_sm_i, na.rm = TRUE) & domain_grid<=max(U_sm_i, na.rm = TRUE)
  cov_sm_compl_mat        <- cov_la_mat[sm_compl_gridloc, sm_compl_gridloc]
  grid_sm_compl_vec       <- domain_grid[sm_compl_gridloc]
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
  if(pre_smooth==TRUE){
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



#' Iterative Reconstruction Algoritm 
#'
#' This function iteratively applies the function reconst_fun() in order to reconstruct the missing parts of a function given the observed parts. 
#' The iterative procedure allows to reconstruct functions when their covariance function cannot be estimated over the total domain. However, the covariance function must be estimated over a sufficiently large part of the domain.  
#' 
#' @param cov_la_mat  Discretized covariance function over \eqn{[a,b]\times[a,b]}{[a,b]x[a,b]}
#' @param domain_grid Equidistant discretization grid in \eqn{[a,b]}{[a,b]}
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
                             domain_grid, 
                             Y_cent_sm_i, 
                             U_sm_i, 
                             K, 
                             fraction   = 0.15,
                             max_rep    = 5,
                             pre_smooth = FALSE,
                             messages   = FALSE
){
  ## ###########################################################
  # t             <- c(82,198, 158)[3]
  # cov_la_mat    <- CV.mat.RE
  # domain_grid   <- x.grid.vec.RE
  # Y_cent_sm_i   <- c(na.omit(Y.s.pred.mat[,t]))
  # U_sm_i        <- x.grid.vec.RE[!is.na(Y.s.pred.mat[,t])]
  # K             <- 2
  # fraction      <- 0.2
  # max_rep       <- 5
  ## ###########################################################
  
  ## extract data:
  grid.len <- length(domain_grid)
  
  
  
  ## Checks #############################################################
  if(grid.len * fraction < K+1){stop("'fraction' too small.")}
  
  
  ## ####################################################################
  ## First Run
  ##
  pred.tmp     <- ReconstPoFD::reconst_fun(cov_la_mat      = cov_la_mat, 
                                           domain_grid     = domain_grid, 
                                           Y_cent_sm_i     = Y_cent_sm_i, 
                                           U_sm_i          = U_sm_i, 
                                           K               = K)
  ##
  U_la_i      <- pred.tmp[['x_reconst']]
  Y_cent_la_i <- pred.tmp[['y_reconst']]
  ##  check-plot:
  # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(domain_grid), type="l"); abline(v=range(domain_grid))
  
  ## ####################################################################
  ## Do we need to predict a further missing *upper* fragment?
  ## ####################################################################
  if(max(U_la_i) !=  max(domain_grid)){
    ## Find the first domain value at which the there is no missing cov-value anymore with 
    ## respect to the largest domain values, i.e., 'cov_la_mat[,grid.len]' (could be done outside of the loop!)
    min.loc         <- which.min(cov_la_mat[,grid.len])
    ##
    ## Is the predicted function already (sufficiently far) beyond 'domain_grid[min.loc]'?
    ## ('Sufficiently far' >= 'grid.len * fraction' grid points)
    check <- max(U_la_i) > (domain_grid[min.loc]) & 
      length(domain_grid[domain_grid >= domain_grid[min.loc] & domain_grid <= max(U_la_i)]) >= (grid.len * fraction)
    stop.pred <- FALSE  # initialization
    if(check == FALSE){
      go    <- TRUE
      count <- 1
      ## Predict further parts of the upper missing fragment:
      while(go){ # print("In while loop for upper pred.")
        ## Select new small intervall: 
        ## [(max(U_la_i) - c(max(U_la_i)-min(U_la_i))*fraction),  max(U_la_i)]
        U_slct_sm2_vec  <- domain_grid <= max(U_la_i) &  domain_grid >= (max(U_la_i) - c(max(U_la_i)-min(U_la_i))*fraction)
        U_sm2_i         <- domain_grid[U_slct_sm2_vec]
        ##
        ## Is the new small interval large enough?
        if(length(U_sm2_i) >= grid.len * fraction){
          ## If yes, start the prediction:
          Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
          Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
          ##
          pred2.tmp    <- ReconstPoFD::reconst_fun(cov_la_mat      = cov_la_mat, 
                                                   domain_grid     = domain_grid, 
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
          # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(domain_grid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
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
          # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(domain_grid), type="l")
          # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
          
          ## Check again:
          ## Is the predicted function now (sufficiently far) beyond 'domain_grid[min.loc]'?
          check <- max(U_la_i) > (domain_grid[min.loc]) & 
            length(domain_grid[domain_grid >= domain_grid[min.loc] & domain_grid <= max(U_la_i)]) >= (grid.len * fraction)
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
    ## Is the predicted function already (sufficiently far) beyond 'domain_grid[min.loc]'?
    ## ('Sufficiently far' >= '('grid.len * fraction' grid points)  
    if(check == TRUE & stop.pred==FALSE){ 
      ## If yes, we can achieve a complete prediction of the missing upper fragments! 
      ##
      U_slct_sm2_vec  <- domain_grid >= domain_grid[min.loc] & domain_grid <= max(U_la_i)
      U_sm2_i         <- domain_grid[U_slct_sm2_vec]
      Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
      Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
      ##
      pred2.tmp    <- ReconstPoFD::reconst_fun(cov_la_mat      = cov_la_mat, 
                                               domain_grid     = domain_grid, 
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
      # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(domain_grid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
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
      # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(domain_grid), type="l")
      # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
    }
  }
  ## ################################################################################
  ## Do we need to predict a further missing *lower* fragement?
  ## ################################################################################
  if(min(U_la_i) !=  min(domain_grid)){
    ## Find the last domain value at which the there is still no missing cov-value with 
    ## respect to the smallest domain value, i.e., 'cov_la_mat[,1]' (could be done outside of the loop!)
    max.loc         <- which.max(cov_la_mat[,1])
    ##
    ## Is the predicted function already (sufficiently far) below 'domain_grid[max.loc]'?
    ## ('Sufficiently far' >= 'grid.len * fraction' grid points)
    check <- min(U_la_i) < (domain_grid[max.loc]) & 
      length(domain_grid[domain_grid <= domain_grid[max.loc] & domain_grid >= min(U_la_i)]) >= (grid.len * fraction)
    stop.pred <- FALSE
    if(check == FALSE){
      go    <- TRUE
      count <- 1
      ## Predict further parts of the lower missing fragment:
      while(go){ # print("In while loop for lower pred.")
        ## Select new small intervall: 
        ## [min(U_la_i), (min(U_la_i) + c(max(U_la_i)-min(U_la_i))*fraction)]
        U_slct_sm2_vec  <- domain_grid >= min(U_la_i) &  domain_grid <= (min(U_la_i) + c(max(U_la_i)-min(U_la_i))*fraction)
        U_sm2_i         <- domain_grid[U_slct_sm2_vec]
        ##
        ## Is the new small interval large enough?
        if(length(U_sm2_i) >= grid.len * fraction){ 
          ## If yes, start the prediction:
          Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
          Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
          ##
          pred2.tmp    <- ReconstPoFD::reconst_fun(cov_la_mat      = cov_la_mat, 
                                                   domain_grid     = domain_grid, 
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
          # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(domain_grid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
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
          # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(domain_grid), type="l")
          # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
          
          ## Check again:
          ## Is the predicted function now (sufficiently far) below 'domain_grid[min.loc]'?
          check <- min(U_la_i) < (domain_grid[max.loc]) & 
            length(domain_grid[domain_grid <= domain_grid[max.loc] & domain_grid >= min(U_la_i)]) >= (grid.len * fraction)
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
    ## Is the predicted function already (sufficiently far) below 'domain_grid[min.loc]'?
    ## ('Sufficiently far' >= '('grid.len * fraction' grid points)  
    if(check == TRUE & stop.pred==FALSE){ 
      ## If yes, we can achieve a complete prediction of the missing lower fragment! 
      ##
      U_slct_sm2_vec  <- domain_grid <= domain_grid[max.loc] & domain_grid >= min(U_la_i)
      U_sm2_i         <- domain_grid[U_slct_sm2_vec]
      Y_slct_sm2_vec  <- U_la_i >= min(U_sm2_i) & U_la_i <= max(U_sm2_i)
      Y_cent_sm2_i    <- Y_cent_la_i[Y_slct_sm2_vec]  
      ##
      pred2.tmp    <- ReconstPoFD::reconst_fun(cov_la_mat      = cov_la_mat, 
                                               domain_grid     = domain_grid, 
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
      # plot(y=Y_cent_la_i, x=U_la_i, xlim=range(domain_grid), ylim=range(Y_cent_la_i,Y_cent_la2_fragm_i), type="l")
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
      # plot(y=Y_cent_la_i,         x=U_la_i, xlim=range(domain_grid), type="l")
      # lines(y=Y_cent_la2_fragm_i, x=U_la2_fragm_i, col="red")
    }
  }
  return(list("y_reconst"=Y_cent_la_i,
              "x_reconst"=U_la_i))
}
