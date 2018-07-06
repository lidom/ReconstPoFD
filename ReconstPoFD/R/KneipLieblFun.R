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





my.fpca <- function(Ly    = Ly, 
                    Lt    = Lu, reconst_fcts) {
  
  n  <- length(Ly)
  m  <- length(Ly[[1]])
  
  ydata = data.frame(".id"     = rep(1:n, each=m),
                     ".index"  = unname(unlist(Lu)),
                     ".values" = unname(unlist(Ly)))
  
  # remaining arguments of the original fpca.sc() function
  Y = NULL; Y.pred = NULL; argvals = NULL; random.int = FALSE;
  nbasis = 10; pve = 0.99; npc = NULL; var = FALSE; simul = FALSE; sim.alpha = 0.95;
  makePD = FALSE; center = TRUE; cov.est.method = 2; integration = "trapezoidal"
  useSymm = FALSE
  #
  
  stopifnot((!is.null(Y) && is.null(ydata)) || (is.null(Y) && !is.null(ydata)))
  
  # if data.frame version of ydata is provided
  sparseOrNongrid <- !is.null(ydata)
  if (sparseOrNongrid) {
    stopifnot(ncol(ydata) == 3)
    stopifnot(c(".id", ".index", ".value") == colnames(ydata))
    stopifnot(is.null(argvals))
    ##
    # Y (nobs x length(argvals)) with NA's if no data at a argvalue
    Y       = refund:::irreg2mat(ydata)# 
    argvals = sort(unique(ydata$.index))
  }
  
  if (is.null(Y.pred)){
    # functions to be reconstructed
    Y.pred = Y[reconst_fcts,,drop=FALSE]
  }
  
  # argvals of observed fragments to be reconstructed
  argvalsO <- vector("list", length = length(reconst_fcts))
  for(i in seq_len(length(reconst_fcts))){
    minmax        <- range(c(stats:na.omit(Y[i,])))
    argvalsO[[i]] <- argvals[argvals >= minmax[1] & argvals <= minmax[2]]
  }
  
  D      = NCOL(Y)     # nobs per function
  I      = NROW(Y)     # number of functions
  I.pred = NROW(Y.pred)# number of functions to be reconstruced
  
  if (is.null(argvals)){
    argvals = seq(0, 1, length = D)
  }
  
  d.vec = rep(argvals, each = I)
  id    = rep(1:I, rep(D, I))
  
  if (center) {
    if (random.int) {
      ri_data <- data.frame(y = as.vector(Y), d.vec = d.vec, id = factor(id))
      gam0 = gamm4::gamm4(y ~ s(d.vec, k = nbasis), random = ~(1 | id), data = ri_data)$gam
      rm(ri_data)
    } else { gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))}
    ## mean
    mu      = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    ## centered functions
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  } else {
    Y.tilde = Y
    mu = rep(0, D)
  }
  
  if (cov.est.method == 2) {
    # smooth raw covariance estimate
    # 1. pointwise (at argvalues) sample covariance matrix (=naive cov-estimator)
    # 2. smooth this matrix
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
      cov.sum[obs.points, obs.points]   = cov.sum[obs.points, obs.points]   + tcrossprod(Y.tilde[i, obs.points])
    }
    G.0       = ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0   = diag(G.0)
    diag(G.0) = NA
    if (!useSymm) {
      row.vec = rep(argvals, each = D)
      col.vec = rep(argvals, D)
      npc.0   = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
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
      ##
      smVCov       <- mgcv::predict.gam(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
      npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
      npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
    }
    
    
    ## estimated covariance function: npc.0
    
    if (makePD) {
      npc.0 <- {
        tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
        as.matrix(tmp$mat)
      }
    }
    
    
    ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman,
    ### Chapter 8)
    w        <- quadWeights(argvals, method = integration)
    Wsqrt    <- diag(sqrt(w))
    Winvsqrt <- diag(1/(sqrt(w)))
    V        <- Wsqrt %*% npc.0 %*% Wsqrt
    evalues  <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    ###
    evalues    = replace(evalues, which(evalues <= 0), 0)
    npc        = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
    efunctions = matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
    evalues    = eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    ## Estimated covariance function ('cov' will be returned, 'cov.hat' will be replaced)
    cov    <-  cov.hat  <-  efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
    ### numerical integration for estimation of sigma2
    T.len  <- argvals[D] - argvals[1]  # total interval length
    T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
    T1.max <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
    DIAG   <- (diag.G0 - diag(cov.hat))[T1.min:T1.max]  # function values
    w2     <- quadWeights(argvals[T1.min:T1.max], method = integration)
    sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)
    ##
  }
  
  ## computations for observed fragments
  scoresO     <- vector("list", length(reconst_fcts))
  evaluesO    <- vector("list", length(reconst_fcts))
  efunctionsO <- vector("list", length(reconst_fcts))
  ##
  for(i in seq_len(length(reconst_fcts))){
    ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman,
    ### Chapter 8)
    w        <- quadWeights(argvalsO[[i]], method = integration)
    Wsqrt    <- diag(sqrt(w))
    Winvsqrt <- diag(1/(sqrt(w)))
    
    locO     <- match(argvalsO[[i]],argvals)
    V        <- Wsqrt %*% npc.0[locO,locO] %*% Wsqrt
    evalues  <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    ###
    evalues    = replace(evalues, which(evalues <= 0), 0)
    npc        = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
    efunctionsO[[i]] = matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
    evaluesO[[i]]    = eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    cov.hat          = efunctionsO[[i]] %*% tcrossprod(diag(evaluesO[[i]], nrow = npc, ncol = npc), efunctionsO[[i]])
    ####
    D.inv          = diag(1/evaluesO[[i]], nrow = npc, ncol = npc)
    Z              = efunctionsO[[i]]
    ##
    Y.tilde        = Y.pred[i,,drop=FALSE] - matrix(mu, 1, D)
    # ##
    # Yhat           = matrix(0, nrow = I.pred, ncol = D)
    # rownames(Yhat) = rownames(Y.pred)
    # colnames(Yhat) = colnames(Y.pred)
    # scores         = matrix(NA, nrow = I.pred, ncol = npc)
    
    # VarMats        = vector("list", I.pred)
    # ##
    # for (i in 1:I.pred){VarMats[[i]] = matrix(NA, nrow = D, ncol = D)}
    # diag.var = matrix(NA, nrow = I.pred, ncol = D)
    # crit.val = rep(0, I.pred)
    
    #for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i, ]))
    if (sigma2 == 0 & length(obs.points) < npc)
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points), ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    ##
    scoresO[[i]] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i, obs.points])
    ##
    #Yhat[i.subj, ] = t(as.matrix(mu)) + scores[i.subj, ] %*% t(efunctions)
    # if (var) {
    #   VarMats[[i.subj]] = sigma2 * Z %*% ZtZ_sD.inv %*% t(Z)
    #   diag.var[i.subj, ] = diag(VarMats[[i.subj]])
    #   if (simul & sigma2 != 0) {
    #     norm.samp = mvrnorm(2500, mu = rep(0, D), Sigma = VarMats[[i.subj]])/matrix(sqrt(diag(VarMats[[i.subj]])), nrow = 2500, ncol = D, byrow = TRUE)
    #     crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)
    #   }
    # }
    #}
    # 
    # ret.objects = c("Yhat", "Y", "scores", "mu", "efunctions", "evalues", "npc",
    #                 "argvals")
    
    ret.objects = c("mu", "cov", "argvals", "argvalsO", "scoresO", "efunctionsO", "evaluesO")
    
    # if (var) {
    #   ret.objects = c(ret.objects, "sigma2", "diag.var", "VarMats")
    #   if (simul)
    #     ret.objects = c(ret.objects, "crit.val")
    # }
    
    ret        = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
    names(ret) = ret.objects
    # class(ret) = "fpca"
    return(ret)
  }