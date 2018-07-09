




DGP        <- "DGP2"
SimDat     <- simuldata(n = 200, m=50, a = 0, b = 1, nRegGrid = 51, DGP=DGP)
#persp(z=var(t(SimDat$Y_true_mat)))
Y_list     <- SimDat[['Y_list']]
U_list     <- SimDat[['U_list']]

par(mfrow=c(1,2))
plot(y=SimDat$Y_true_mat[,1],x=SimDat$U_true_mat[,1], col="red", type="o")

# ##
# test   <- my.fpca(Ly=Y_list, Lu=U_list, CEscores = ifelse(DGP=="DGP1", TRUE, FALSE), reconst_fcts=1)
# ##
# Yhat   <-  t(as.matrix(test$muO[[1]])) + t(test$scoresO[[1]])   %*% t(test$efunctionsO[[1]])
# YhatCE <-  t(as.matrix(test$muO[[1]])) + t(test$CEscoresO[[1]]) %*% t(test$efunctionsO[[1]])
# lines(y=c(Yhat),   x=test$argvalsO[[1]], type="l", lwd=4)
# lines(y=c(YhatCE), x=test$argvalsO[[1]], lty=2)
# ##
# for(k in c(2:6)){
# test2 <- reconstKneipLiebl_fun(mu=test$mu, cov = test$cov, argvals = test$argvals, argvalsO = test$argvalsO[[1]], 
#                       scoresO = test$scoresO[[1]], efunctionsO = test$efunctionsO[[1]], evaluesO = test$evaluesO[[1]], 
#                       fragmO =  c(na.omit(SimDat$Y_mat[,1])), K=k)
# lines( y=c(test2$y_reconst),   x=test2$x_reconst, type="o",pch=as.character(k),lwd=1.2, col=gray(.75))
# }
##
# CEtest2 <- reconstKneipLiebl_fun(mu=test$mu, cov = test$cov, argvals = test$argvals, argvalsO = test$argvalsO[[1]], 
#                                scoresO = test$CEscoresO[[1]], efunctionsO = test$efunctionsO[[1]], evaluesO = test$evaluesO[[1]], K=3)
# lines( y=c(CEtest2$y_reconst), x=CEtest2$x_reconst, lty=4)

for(k in c(2:6)){
  #Ly=Y_list; Lu=U_list; reconst_fcts=1; CEscores = ifelse(DGP=="DGP1", TRUE, FALSE); center = TRUE; K=3
test3 <- reconstructKneipLiebl(Ly=Y_list, Lu=U_list, reconst_fcts=1, method=2, CEscores = ifelse(DGP=="DGP1", TRUE, FALSE), K=k)
lines( y=c(test3$Y_reconst_list[[1]]),   x=test3$U_reconst_list[[1]], type="o",pch=as.character(k),lwd=1.2, col=gray(.75))
}
lines(y=SimDat$Y_mat[,1],x=SimDat$U_mat[,1], col="black", lwd=4)


plot(y=SimDat$Y_true_mat[,1],x=SimDat$U_true_mat[,1], col="red", type="o")
##
result        <- reconstructKraus(X_mat = Y_mat, reconst_fcts = 1)
Y_reconst_mat <- result[['X_reconst_mat']]
lines( y=c(Y_reconst_mat),   x=seq(0,1,len=length(c(Y_reconst_mat))), col=gray(.75),lwd=2, type="o")
lines(y=SimDat$Y_mat[,1],x=SimDat$U_mat[,1], col="black", lwd=4)


reconstructKneipLiebl <- function(Ly,
                                  Lu, 
                                  reconst_fcts = NULL,
                                  method       = c(1,2)[2],
                                  CEscores     = TRUE, 
                                  center       = TRUE,
                                  K            = NULL
){
  ##
  n  <- length(Ly)
  if(is.null(reconst_fcts)){ reconst_fcts <- 1:n }
  ##
  if(method == 1){
    ## method 1: Assumes fully observed fragments and aligns the reconstructed parts to the observed fragment 
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center)
    ##
    Y_reconst_list <- vector("list", length(reconst_fcts))
    U_reconst_list <- vector("list", length(reconst_fcts))
    K_vec          <- rep(NA, length(reconst_fcts))
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( any(fpca_obj$obs_argvalsO[[i]] != fpca_obj$argvalsO[[i]]) ){ stop("The fragment must be fully observed (obs_argvalsO == argvalsO).") }
      ##
      fragmO     <-  c(na.omit(fpca_obj$Y[i,]))
      ##
      ## hier AIC (o.ä.)
      K_vec[i]   <- K
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          cov         = fpca_obj$cov, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efunctionsO = fpca_obj$efunctionsO[[i]], 
                                          evaluesO    = fpca_obj$evaluesO[[i]], 
                                          fragmO      = fragmO, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 2){
    ## method 2: Assumes fully observed fragments, but does not aligns the reconstructed parts to the observed fragment 
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center)
    ##
    Y_reconst_list <- vector("list", length(reconst_fcts))
    U_reconst_list <- vector("list", length(reconst_fcts))
    K_vec          <- rep(NA, length(reconst_fcts))
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      ## hier AIC (o.ä.)
      K_vec[i]   <- K
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          cov         = fpca_obj$cov, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efunctionsO = fpca_obj$efunctionsO[[i]], 
                                          evaluesO    = fpca_obj$evaluesO[[i]], 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  ##
  return(list(
    "Y_reconst_list"  = Y_reconst_list,
    "U_reconst_list"  = U_reconst_list,
    "K_median"        = stats::median(K_vec)
  ))
}




reconstKneipLiebl_fun <- function(mu, cov, argvals, argvalsO, scoresO, efunctionsO, evaluesO, fragmO=NULL, K=NULL){
  ##
  K             <- min(length(evaluesO), K)
  efun_reconst  <- matrix(NA, nrow=length(argvals), ncol=K)
  locO          <- match(argvalsO, argvals)
  ##
  for(k in seq_len(K)){
    if(evaluesO[k] > 0){
      efun_reconst[,k] <- apply(X   = cov[locO,,drop=FALSE], MARGIN = 2,
                                FUN = function(x){pracma::trapz(x=argvalsO, efunctionsO[,k] * x)})
      efun_reconst[,k] <- efun_reconst[,k] / evaluesO[k]
    }else{
      efun_reconst[,k] <- 0
    }
  }
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(scoresO[1:K]) %*% t(efun_reconst) ))
  ##
  ## Align to observed or estimated fragement
  if( !is.null(fragmO) ){
    reconstr[1:min(locO)]               <- reconstr[1:min(locO)]               + fragmO[1]              - reconstr[min(locO)]
    reconstr[max(locO):length(argvals)] <- reconstr[max(locO):length(argvals)] + fragmO[length(fragmO)] - reconstr[max(locO)]
    reconstr[locO]                      <- fragmO
  }
  ##
  ## ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  ## ######################
}


my.fpca <- function(Ly, Lu, reconst_fcts, CEscores = TRUE, center = TRUE, pve = NULL) {
  
  n      <- length(Ly)
  id_vec <- NULL
  for(i in 1:n){id_vec <- c(id_vec, rep(i, length(Ly[[i]])))}
  ydata  <-  data.frame(".id"    = id_vec, 
                        ".index" = unname(unlist(Lu)), 
                        ".value" = unname(unlist(Ly)))
  ##
  if(is.null(pve)){pve <- 1}
  ##
  nbasis         = 10
  maxbins        = 1000
  makePD         = FALSE
  cov.est.method = 2
  integration    = "trapezoidal"
  useSymm        = FALSE
  ##
  # Y (nobs x length(argvals)) with NA's if no data at a argvalue
  Y        <- refund:::irreg2mat(ydata, binning = TRUE, maxbins = maxbins) 
  argvals  <- as.numeric(colnames(Y))
  # functions to be reconstructed
  Y.pred   <- Y[reconst_fcts,,drop=FALSE]
  # argvals of observed fragments to be reconstructed
  argvalsO <- vector("list", length = length(reconst_fcts))
  for(i in seq_len(length(reconst_fcts))){
    minmax        <- range(as.numeric(names(c(na.omit(Y.pred[i,])))))
    argvalsO[[i]] <- argvals[argvals >= minmax[1] & argvals <= minmax[2]]
  }
  ##
  D      <- NCOL(Y)      # nobs per function
  I      <- NROW(Y)      # number of functions
  I.pred <- NROW(Y.pred) # number of functions to be reconstruced
  ##
  d.vec <- rep(argvals, each = I)
  id    <- rep(1:I, rep(D, I))
  
  if (center) {
    ## mean
    gam0    <-  mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu      <-  mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    Y.tilde <- Y - matrix(mu, I, D, byrow = TRUE)
  } else {
    Y.tilde <- Y
    mu      <- rep(0, D)
  }
  ## plot(x=argvals,y=mu, type="l")

  
  if (cov.est.method == 2) {
    # smooth raw covariance estimate
    # 1. pointwise (at argvalues) sample covariance matrix (=naive cov-estimator)
    # 2. smooth this matrix
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
      cov.sum[  obs.points, obs.points] <- cov.sum[  obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
    }
    G.0       <- ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0   <- diag(G.0)
    diag(G.0) <- NA
    if (!useSymm) {
      row.vec <- rep(argvals, each = D)
      col.vec <- rep(argvals, D)
      npc.0   <- matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
                                                   weights = as.vector(cov.count)), newdata = data.frame(row.vec = row.vec,
                                                                                                         col.vec = col.vec)), D, D)
      npc.0 = (npc.0 + t(npc.0))/2
    } else {
      use          <- upper.tri(G.0, diag = TRUE)
      use[2, 1]    <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
      usecov.count <- cov.count
      usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
      usecov.count <- as.vector(usecov.count)[use]
      use          <- as.vector(use)
      vG.0         <- as.vector(G.0)[use]
      row.vec      <- rep(argvals, each = D)[use]
      col.vec      <- rep(argvals, times = D)[use]
      mCov         <- mgcv::gam(vG.0 ~ te(row.vec, col.vec, k = nbasis), weights = usecov.count)
      npc.0        <- matrix(NA, D, D)
      spred        <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
      tpred        <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
      # Estimated covariance function:
      smVCov       <- mgcv::predict.gam(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
      npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
      npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
      # slct <- seq.int(1,length(argvals),len=25)
      # persp(z=npc.0[slct,slct],x=argvals[slct],y=argvals[slct])
    }

    if (makePD) {
      npc.0 <- {
        tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
        as.matrix(tmp$mat)
      }
    }
    
    # Nnumerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
    w          <- refund:::quadWeights(argvals, method = integration)
    Wsqrt      <- diag(sqrt(w))
    Winvsqrt   <- diag(1/(sqrt(w)))
    V          <- Wsqrt %*% npc.0 %*% Wsqrt
    evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues    <- replace(evalues, which(evalues <= 0), 0)
    npc        <- min(which(cumsum(evalues)/sum(evalues) >= pve))
    efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    # Estimated covariance function ('cov' will be returned, 'cov.hat' will be replaced)
    cov        <- cov.hat  <-  efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
    # Numerical integration for estimation of sigma2
    T.len      <- argvals[D] - argvals[1]  # total interval length
    T1.min     <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
    T1.max     <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
    DIAG       <- (diag.G0 - diag(cov.hat))[T1.min:T1.max]  # function values
    w2         <- refund:::quadWeights(argvals[T1.min:T1.max], method = integration)
    sigma2     <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)
    ##
  }
  
  ## computations for observed fragments
  muO          <- vector("list", length(reconst_fcts))
  scoresO      <- vector("list", length(reconst_fcts))
  CEscoresO    <- vector("list", length(reconst_fcts))
  evaluesO     <- vector("list", length(reconst_fcts))
  efunctionsO  <- vector("list", length(reconst_fcts))
  obs_argvalsO <- vector("list", length(reconst_fcts))
  ##
  for(i in seq_len(length(reconst_fcts))){# i <- 1
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    w                 <- refund:::quadWeights(argvalsO[[i]], method = integration)
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    locO              <- match(argvalsO[[i]],argvals)
    # CovOO
    V                 <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalues           <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues           <- replace(evalues, which(evalues <= 0), 0)
    npc               <- min(which(cumsum(evalues)/sum(evalues) >= pve))
    efunctionsO[[i]]  <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evaluesO[[i]]     <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    D.inv             <- diag(1/evaluesO[[i]], nrow = npc, ncol = npc)
    Z                 <- efunctionsO[[i]]
    Y.cent            <- Y.pred[i,,drop=FALSE] - matrix(mu, 1, D)
    obs_locO          <- match(names(c(na.omit((Y.pred[i,])))), as.character(argvalsO[[i]]))
    obs_argvalsO[[i]] <- argvalsO[[i]][obs_locO] 
    ## CEScores (i.e., PACE-Scores)
    if(CEscores){
      if(sigma2 == 0){warning("Measurement error estimated to be zero; CEscores cannot be estimated.")}
      if(length(obs_locO) < npc){warning("There are fewer observed points than PCs; CEscores cannot be estimated.")}
      Zcur           <- Z[obs_locO,,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv)
      CEscoresO[[i]] <- ZtZ_sD.inv %*% t(Zcur) %*% c(na.omit(Y.cent[i,]))
    } else {
      CEscoresO[[i]] <- NA
    }
    ## Classical scores (intergration)
    scoresO[[i]] <- apply(X      = efunctionsO[[i]][obs_locO,], 
                          MARGIN = 2, 
                          FUN    = function(ef){pracma::trapz(y=ef*c(na.omit(Y.cent[i,])),x=argvalsO[[i]][obs_locO])})
    ##
    muO[[i]]     <- mu[locO]
    ##
    # Yhat    <- t(as.matrix(muO[[i]])) + t(scoresO[[i]]) %*% t(efunctionsO[[i]])
    # Yhat_CE <- t(as.matrix(muO[[i]])) + t(CEscoresO[[i]]) %*% t(efunctionsO[[i]])
    # plot( y=c(Yhat),    x=argvalsO[[i]], type="l")
    # lines(y=c(Yhat_CE), x=argvalsO[[i]])
  }
  ## Return results
  ret.objects <- c("Y", "mu", "muO", "cov", "argvals", "argvalsO", "obs_argvalsO", "CEscoresO", "scoresO", "efunctions", "efunctionsO",  "evalues", "evaluesO")
  ret         <- lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret)  <- ret.objects
  return(ret)
}