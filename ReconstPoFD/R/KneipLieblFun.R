
B <- 75
KL_m1 <- rep(NA,B)
KL_m2 <- rep(NA,B)
Kraus <- rep(NA,B)
PACE  <- rep(NA,B)
for(r in 1:B){
DGP        <- "DGP2"
while(TRUE){
SimDat   <- simuldata(n = 100, m=50, a = 0, b = 1, nRegGrid = 51, DGP=DGP)
if(any(is.na(SimDat$Y_mat[,1]))){break}
}
#persp(z=var(t(SimDat$Y_true_mat)))
Y_list     <- SimDat[['Y_list']]
U_list     <- SimDat[['U_list']]
###
test3 <- reconstructKneipLiebl(Ly=Y_list, Lu=U_list, reconst_fcts=1, method = 1)
test4 <- reconstructKneipLiebl(Ly=Y_list, Lu=U_list, reconst_fcts=1, method = 2)
###
result        <- reconstructKraus(X_mat = SimDat$Y_mat, reconst_fcts = 1)
Y_reconst_mat <- result[['X_reconst_mat']]
df            <- round(result$df_median,0)
####
result_PACE <- fdapace::FPCA(Ly = Y_list, Lt = U_list)
Y_PACE_mat <- t(fitted(result_PACE))[,1]
K_PACE <- length(result_PACE$lambda)
###
par(mfrow=c(1,3))
plot(y=SimDat$Y_true_mat[,1],x=SimDat$U_true_mat[,1], col="red", type="o",main="Kneip Liebl")
lines( y=c(test3$Y_reconst_list[[1]]),   x=test3$U_reconst_list[[1]], type="o",pch=as.character(test3$K_median),lwd=1.2, col=gray(.65))
lines( y=c(test4$Y_reconst_list[[1]]),   x=test4$U_reconst_list[[1]], type="o",pch=as.character(test4$K_median),lwd=1.2, col=gray(.85))
lines(y=SimDat$Y_mat[,1],x=SimDat$U_mat[,1], col="black", lwd=4)
#####
plot(y=SimDat$Y_true_mat[,1],x=SimDat$U_true_mat[,1], col="red", type="o",main="Kraus")
lines( y=c(Y_reconst_mat),x=seq(0,1,len=length(c(Y_reconst_mat))), col=gray(.75),lwd=2, type="o",pch=as.character(df))
lines(y=SimDat$Y_mat[,1],x=SimDat$U_mat[,1], col="black", lwd=4)
#####
plot(y=SimDat$Y_true_mat[,1],x=SimDat$U_true_mat[,1], col="red", type="o",main="PACE")
lines( y=Y_PACE_mat,x=seq(0,1,len=length(c(Y_PACE_mat))), col=gray(.75),lwd=2, type="o",pch=as.character(K_PACE))
lines(y=SimDat$Y_mat[,1],x=SimDat$U_mat[,1], col="black", lwd=4)
#
Sys.sleep(1.5)

#KneipLiebl
KL_m1[r] <- round(sum(c(test3$Y_reconst_list[[1]]-SimDat$Y_true_mat[,1])^2),2)#m1
KL_m2[r] <- round(sum(c(test4$Y_reconst_list[[1]]-SimDat$Y_true_mat[,1])^2),2)#m2
#Kraus:
Kraus[r] <- round(sum(c(c(Y_reconst_mat)         -SimDat$Y_true_mat[,1])^2),2)
PACE[r]  <- round(sum(c(c(Y_PACE_mat)            -SimDat$Y_true_mat[,1])^2),2)
cat("r=",r,"\n")
}
mean(KL_m1)
mean(KL_m2)
mean(Kraus)
mean(PACE)

reconstructKneipLiebl <- function(Ly,
                                  Lu, 
                                  reconst_fcts = NULL,
                                  method       = c(1,2)[1],
                                  center       = TRUE,
                                  K            = NULL){
  ##
  n  <- length(Ly)
  if(is.null(reconst_fcts)){ reconst_fcts <- 1:n }
  ##
  if(method == 1){
    ## method 1: This method assumes fully observed fragments 'fragmO' and 
    ## aligns the reconstructed parts to the observed fragment 'fragmO'.
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center)
    ##
    Y_reconst_list <- vector("list", length(reconst_fcts))
    U_reconst_list <- vector("list", length(reconst_fcts))
    K_vec          <- rep(NA,        length(reconst_fcts))
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( any(fpca_obj$obs_argvalsO[[i]] != fpca_obj$argvalsO[[i]]) ){ stop("The fragment must be fully observed (obs_argvalsO == argvalsO).") }
      ##
      fragmO     <- c(na.omit(c(fpca_obj$Y[i,])))
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 1)
      }else{K_vec[i] <- K}
      ##
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
    ## method 2: This method assumes fully observed fragments (as in method 1), but 
    ## does not align the reconstructed parts to the observed fragments. 
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
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 1)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          cov         = fpca_obj$cov, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efunctionsO = fpca_obj$efunctionsO[[i]], 
                                          evaluesO    = fpca_obj$evaluesO[[i]], 
                                          fragmO      = NULL, 
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


gcvKneipLiebl <- function(fpca_obj, argvalsO, method){
  ##
  Y          <- fpca_obj$Y # data pre-processed by irreg2mat()
  mu         <- fpca_obj$mu
  cov        <- fpca_obj$cov
  sigma2     <- fpca_obj$sigma2
  argvals    <- fpca_obj$argvals
  a          <- min(argvals)
  b          <- max(argvals)
  ##
  ## a function is considered as complete if: 
  ## 1. its first and its last observation are non-missing
  ## 2. not considered so far
  compl_fcts <- apply(Y, 1, function(x){
    !is.na(utils::head(x, 1)) & !is.na(utils::tail(x, 1)) 
    # & length(c(stats::na.omit(x)))>=floor(.8*length(argvals))
  })
  n_compl <- length(compl_fcts[compl_fcts==TRUE])
  ##
  if(n_compl <= 10){warning("Very few (<=10) complete functions; do not trust the GCV-result.")}
  if(n_compl <= 1){stop("Too few complete functions.")}
  ##
  # functions to be reconstructed
  Y.compl       <- Y[compl_fcts,,drop=FALSE]
  # make artifical fragements (all over argvalsO)
  locO          <- match(argvalsO, argvals)
  locM          <- c(1:length(argvals))[-locO]
  Y.pred        <- Y.compl
  Y.pred[,locM] <- NA
  ##
  # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
  muO               <- mu[locO]
  w                 <- quadWeights(argvalsO, method = "trapezoidal")
  Wsqrt             <- diag(sqrt(w))
  Winvsqrt          <- diag(1/(sqrt(w)))
  # CovOO
  V                 <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
  evalues           <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  evalues           <- replace(evalues, which(evalues <= 0), 0)
  npc               <- length(evalues[evalues>0])
  efunctionsO       <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
  evaluesO          <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  D.inv             <- diag(1/evaluesO, nrow = npc, ncol = npc)
  Z                 <- efunctionsO
  ##
  rss_mat           <- matrix(NA, n_compl, npc)
  ##
  for(i in seq_len(n_compl)){# i <- 1
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, ncol(Y)))
    obs_locO          <- match(names(c(na.omit((Y.pred[i,])))), as.character(argvalsO))
    ## 
    if(method == 3){
      ## CEScores (i.e., PACE-Scores)
      if(sigma2           ==  0){warning("Measurement error estimated to be zero; CEscores cannot be estimated.")}
      if(length(obs_locO) < npc){warning("There are fewer observed points than PCs; CEscores cannot be estimated.")}
      Zcur           <- Z[obs_locO,,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv)
      CEscoresO      <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    } else {
      CEscoresO      <- NA
    }
    ## Classical scores (via intergration)
    if(any(method == c(1,2))){
      scoresO      <- apply(X      = efunctionsO[obs_locO,,drop=FALSE], 
                            MARGIN = 2, 
                            FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=argvalsO[obs_locO])})
    }
    ##
    for(k in 1:npc){
      ##
      if(method == 1){
        result_tmp <- reconstKneipLiebl_fun(mu=mu,cov=cov,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,  efunctionsO=efunctionsO,evaluesO=evaluesO,fragmO=Y.pred[i,locO], K=k)
      }
      if(method == 2){
        result_tmp <- reconstKneipLiebl_fun(mu=mu,cov=cov,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,  efunctionsO=efunctionsO,evaluesO=evaluesO,fragmO=NULL,       K=k)
      }
      if(method == 3){
        result_tmp <- reconstKneipLiebl_fun(mu=mu,cov=cov,argvals=argvals,argvalsO=argvalsO,scoresO=CEscoresO,efunctionsO=efunctionsO,evaluesO=evaluesO,fragmO=PRESMOOTH,  K=k)
      }
      rss_mat[i,k] <- sum((result_tmp[['y_reconst']][locM] - Y.compl[i,locM])^2)
    }
  }
  
  gcv_k_vec <- colSums(rss_mat)/((1-1:npc/n_compl)^2)
  K.gcv     <- which.min(gcv_k_vec)
  ##
  return(K.gcv)
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
  if( !is.null(fragmO) ){
    ## Aligning the reconstructed parts to fragm0
    reconstr[1:min(locO)]               <- reconstr[1:min(locO)]               + fragmO[1]              - reconstr[min(locO)]
    reconstr[max(locO):length(argvals)] <- reconstr[max(locO):length(argvals)] + fragmO[length(fragmO)] - reconstr[max(locO)]
    reconstr[locO]                      <- fragmO
  }
  ##
  ## ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  ## ######################
}


my.fpca <- function(Ly, Lu, reconst_fcts = NULL, CEscores = TRUE, center = TRUE){
  
  n      <- length(Ly)
  id_vec <- NULL
  for(i in 1:n){id_vec <- c(id_vec, rep(i, length(Ly[[i]])))}
  ydata  <-  data.frame(".id"    = id_vec, 
                        ".index" = unname(unlist(Lu)), 
                        ".value" = unname(unlist(Ly)))
  ##
  nbasis         <- 10
  maxbins        <- 1000
  cov.est.method <- 1     # at the moment there is only one cov.est.method
  useSymm        <- FALSE # if true, there will be no smoothing accross the diagonal
  makePD         <- FALSE # if true, the NP-estimate of cov is forced to be positive definite
  ##
  # Y (nobs x length(argvals)) with NA's if no data at a argvalue
  Y        <- irreg2mat(ydata, binning = TRUE, maxbins = maxbins) 
  argvals  <- as.numeric(colnames(Y))
  ##
  # functions to be reconstructed
  if(is.null(reconst_fcts)){reconst_fcts <- 1:n}
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
  
  if(center){
    ## mean
    gam0    <-  mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu      <-  mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    Y.tilde <- Y - matrix(mu, I, D, byrow = TRUE)
  }else{
    Y.tilde <- Y
    mu      <- rep(0, D)
  }
  ## plot(x=argvals,y=mu, type="l")
  ##
  if(cov.est.method == 1) {
    # smooth raw covariance estimate
    # 1. pointwise (at argvalues) sample covariance matrix (=naive cov-estimator)
    # 2. smooth this matrix
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for(i in 1:I){
      obs.points = which(!is.na(Y[i, ]))
      cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
      cov.sum[  obs.points, obs.points] <- cov.sum[  obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
    }
    G.0       <- ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0   <- diag(G.0)
    diag(G.0) <- NA
    if(!useSymm){
      row.vec <- rep(argvals, each = D)
      col.vec <- rep(argvals, D)
      npc.0   <- matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
                                                    weights = as.vector(cov.count)), 
                                          newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
      npc.0 = (npc.0 + t(npc.0))/2
    }else{
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
    if(makePD){
      npc.0 <- {
        tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
        as.matrix(tmp$mat)
      }
    }
    # Nnumerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
    w          <- quadWeights(argvals, method = "trapezoidal")
    Wsqrt      <- diag(sqrt(w))
    Winvsqrt   <- diag(1/(sqrt(w)))
    V          <- Wsqrt %*% npc.0 %*% Wsqrt
    evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues    <- replace(evalues, which(evalues <= 0), 0)
    npc        <- length(evalues[evalues>0])
    efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    # Estimated covariance function
    cov        <- efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
    # Numerical integration for estimation of sigma2
    T.len      <- argvals[D] - argvals[1]  # total interval length
    T1.min     <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
    T1.max     <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
    DIAG       <- (diag.G0 - diag(cov))[T1.min:T1.max]  # function values
    w2         <- quadWeights(argvals[T1.min:T1.max], method = "trapezoidal")
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
    w                 <- quadWeights(argvalsO[[i]], method = "trapezoidal")
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    locO              <- match(argvalsO[[i]],argvals)
    # CovOO
    V                 <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalues           <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues           <- replace(evalues, which(evalues <= 0), 0)
    npc               <- length(evalues[evalues>0])  
    efunctionsO[[i]]  <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evaluesO[[i]]     <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    D.inv             <- diag(1/evaluesO[[i]], nrow = npc, ncol = npc)
    Z                 <- efunctionsO[[i]]
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
    obs_locO          <- match(names(c(na.omit((Y.pred[i,])))), as.character(argvalsO[[i]]))
    obs_argvalsO[[i]] <- argvalsO[[i]][obs_locO]
    ## 
    if(CEscores){
      ## CEScores (i.e., PACE-Scores)
      if(sigma2           ==  0){warning("Measurement error estimated to be zero; CEscores cannot be estimated.")}
      if(length(obs_locO) < npc){warning("There are fewer observed points than PCs; CEscores cannot be estimated.")}
      Zcur           <- Z[obs_locO,,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv)
      CEscoresO[[i]] <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    } else {
      CEscoresO[[i]] <- NA
    }
    ## Classical scores (via intergration)
    scoresO[[i]] <- apply(X      = efunctionsO[[i]][obs_locO,,drop=FALSE], 
                          MARGIN = 2, 
                          FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=obs_argvalsO[[i]])})
    ##
    muO[[i]]     <- mu[locO]
    ##
    # Yhat    <- t(as.matrix(muO[[i]])) + t(scoresO[[i]]) %*% t(efunctionsO[[i]])
    # Yhat_CE <- t(as.matrix(muO[[i]])) + t(CEscoresO[[i]]) %*% t(efunctionsO[[i]])
    # plot( y=c(Yhat),    x=argvalsO[[i]], type="l")
    # lines(y=c(Yhat_CE), x=argvalsO[[i]])
  }
  ## Return results
  ret.objects <- c("Y", "mu", "muO", "cov", "sigma2",
                   "argvals",    "argvalsO", "obs_argvalsO",
                   "CEscoresO",  "scoresO", 
                   "efunctions", "efunctionsO", 
                   "evalues",    "evaluesO")
  ret         <- lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret)  <- ret.objects
  return(ret)
}


irreg2mat <- function(ydata, binning = FALSE, maxbins = 1000){
  ##
  ydata <- ydata[stats::complete.cases(ydata), ]
  nobs  <- length(unique(ydata$.id))
  newid <- as.numeric(as.factor(ydata$.id))
  bins  <- sort(unique(ydata$.index))
  if(binning && (length(bins) > maxbins)){
    binvalues <- seq((1 - 0.001 * sign(bins[1])) * bins[1], 
                     (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)], 
                     l = maxbins + 1)
    bins      <- binvalues
    binvalues <- utils::head(stats::filter(binvalues, c(0.5, 0.5)), -1)
  }else{
    binvalues <- bins
    bins      <- c((1 - 0.001 * sign(bins[1])) * bins[1], bins[-length(bins)], 
                   (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)])
    if(bins[1] == 0           ){bins[1]            <- -0.001}
    if(bins[length(bins)] == 0){bins[length(bins)] <-  0.001}
  }
  newindex         <- cut(ydata$.index, breaks = bins, include.lowest = TRUE)
  Y                <- matrix(NA, nrow = nobs, ncol = nlevels(newindex))
  colnames(Y)      <- binvalues
  attr(Y, "index") <- binvalues
  ## If there are more than one data-point within a bin, 
  ## then only one of these is used (the last one).
  Y[cbind(newid, as.numeric(newindex))] <- ydata$.value
  ##
  return(Y)
}

quadWeights <- function(argvals, method = "trapezoidal"){
  ret <- switch(method, 
    trapezoidal = {D <- length(argvals); 1/2 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 2)], argvals[D] - argvals[D - 1])}, 
    midpoint    = c(0, diff(argvals)), 
    stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule")
  )
  ##
  return(ret)
}