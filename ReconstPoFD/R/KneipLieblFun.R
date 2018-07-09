reconstructKneipLiebl <- function(Ly,
                                  Lu, 
                                  K            = NULL,
                                  reconst_fcts = NULL){
  ##
  n  <- length(Ly)
  m  <- length(Ly[[1]])
  if(is.null(reconst_fcts)){
    reconst_fcts <- 1:n
  }
  ##
  
}




reconstKneipLiebl_fun <- function(cov, argvals, argvalsO, scoresO, efunctionsO, evaluesO, K){
  ##
  K_max_O       <- length(evaluesO)
  efun_reconst  <- matrix(NA, nrow=length(argvals), ncol=K_max_O)
  locO          <- match(argvalsO, argvals)
  ##
  for(k in seq_len(K_max_O)){
    if(evaluesO[k] > 0){
      efun_reconst[,k] <- apply(X      = cov[locO,,drop=FALSE],
                                MARGIN = 2,
                                FUN    = function(x){pracma::trapz(x=argvalsO, efunctionsO[,k] * x)})
      ##
      efun_reconst[,k] <- efun_reconst[,k] / evaluesO[k]
    }else{
      efun_reconst[,k] <- 0
    }
  }
  
  tmp      <- matrix(rep(scoresO[1:K], each=nrow(efun_reconst)), nrow = nrow(efun_reconst), ncol = ncol(K)) * efun_reconst[,1:K,drop=FALSE]
  reconstr <- rowSums(tmp)
  ## ######################
  return(list("y_reconst"  = c(reconstr),
              "x_reconst"  = c(workGrid)))
  ## ######################
}


SimDat     <- simuldata(n = 100, a = 0, b = 1, nRegGrid = 15, DGP="DGP1")

persp(z=var(t(SimDat$Y_true_mat)),x=seq(0,1,len=15),y=seq(0,1,len=15))

Y_list     <- SimDat[['Y_list']]
U_list     <- SimDat[['U_list']]


Ly=Y_list; Lu=U_list; reconst_fcts=1
test <- my.fpca(Ly=Y_list, Lu=U_list, reconst_fcts=1)


Yhat <-  t(as.matrix(test$muO[[1]])) + t(test$scoresO[[1]]) %*% t(test$efunctionsO[[1]])
plot(y=c(Yhat), x=test$argvalsO[[1]],type="l")
lines(y=SimDat$Y_true_mat[,1],x=SimDat$U_true_mat[,1], col="red")


my.fpca <- function(Ly, Lu, reconst_fcts, center = TRUE, pve = 0.99, npc = NULL) {
  
  n      <- length(Ly)
  id_vec <- NULL
  for(i in 1:n){id_vec <- c(id_vec, rep(i, length(Ly[[i]])))}
  ydata  <-  data.frame(".id"    = id_vec, ".index" = unname(unlist(Lu)), ".value" = unname(unlist(Ly)))
  
  #
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
    npc        <- ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
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
  muO         <- vector("list", length(reconst_fcts))
  scoresO     <- vector("list", length(reconst_fcts))
  evaluesO    <- vector("list", length(reconst_fcts))
  efunctionsO <- vector("list", length(reconst_fcts))
  ##
  for(i in seq_len(length(reconst_fcts))){# i <- 1
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    w                <- refund:::quadWeights(argvalsO[[i]], method = integration)
    Wsqrt            <- diag(sqrt(w))
    Winvsqrt         <- diag(1/(sqrt(w)))
    locO             <- match(argvalsO[[i]],argvals)
    # CovOO
    V                <- Wsqrt %*% npc.0[locO,locO] %*% Wsqrt
    evalues          <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues          <- replace(evalues, which(evalues <= 0), 0)
    npc              <- ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
    efunctionsO[[i]] <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evaluesO[[i]]    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    cov.hat          <- efunctionsO[[i]] %*% tcrossprod(diag(evaluesO[[i]], nrow = npc, ncol = npc), efunctionsO[[i]])
    D.inv            <- diag(1/evaluesO[[i]], nrow = npc, ncol = npc)
    Z                <- efunctionsO[[i]]
    Y.cent           <- Y.pred[i,,drop=FALSE] - matrix(mu, 1, D)
    obs_locO         <- match(names(c(na.omit((Y.pred[i,])))), as.character(argvalsO[[i]]))
    ## CEScores (i.e., PACE-Scores)
    if(error){
      if (sigma2 == 0){stop("Measurement error estimated to be zero; CEscores cannot be estimated.")}
      if (length(obs_locO) < npc){stop("There are fewer observed points than PCs; CEscores cannot be estimated.")}
      Zcur           <- Z[obs_locO,,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv)
      CEscoresO[[i]] <- ZtZ_sD.inv %*% t(Zcur) %*% c(na.omit(Y.cent[i,]))
    } else {
      CEscoresO[[i]] <- NA
    }
    ## Classical scores (intergration)
    scoresO[[i]] <- apply(X      = efunctionsO[[i]][obs_locO,], 
                          MARGIN = 2, 
                          FUN    = function(ef){pracma::trapz(y=ef*c(na.omit(Y.cent[i,])),x=argvalsO[[i]][obs_locO])}
    
    muO[[i]]     <- mu[locO]
    ##
    # Yhat    <- t(as.matrix(muO[[i]])) + t(scoresO[[i]]) %*% t(efunctionsO[[i]])
    # Yhat_CE <- t(as.matrix(muO[[i]])) + t(CEscoresO[[i]]) %*% t(efunctionsO[[i]])
    # plot( y=c(Yhat),    x=argvalsO[[i]], type="l")
    # lines(y=c(Yhat_CE), x=argvalsO[[i]])
  }
  ## Return results
  ret.objects <- c("mu", "muO", "cov", "argvals", "argvalsO", "CEscoresO", "scoresO", "efunctions", "efunctionsO",  "evalues", "evaluesO")
  ret         <- lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret)  <- ret.objects
  return(ret)
}