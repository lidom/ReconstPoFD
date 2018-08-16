## R-packages 
library("devtools")
# install necessary packages:
# install_github("lidom/ReconstPoFD/ReconstPoFD")
# install.packages("fdapace")

library("ReconstPoFD")  # contains the function 'reconstuct()'

## Location to store the simulation results:
setwd("/home/dom/ownCloud/Kneip_Liebl_Reconstruction/Simulation_Submission_2")

## #######################################
## Number of target functions 
R         <-   50
## Number of MC-Repetitions per target function
B         <-   100
## #######################################

## #######################################
## Lower-value of the total domain 
a             <-   0             
## Upper-value of the total domain
b             <-   1
## #######################################

set.seed(1234)

maxbins <- 100

(Start.Time <- Sys.time())

for(DGP in c('DGP1','DGP2','DGP3')[1]){
  for(n in c(50, 100)){
    if(DGP=='DGP1'){m_seq <- c(15, 30)}else{m_seq <- NA}
    for(m in m_seq){
      
      ## DGP <- 'DGP1'; n <- 50; m <- 10; B <- 2; R <- 2
      
      ## #######################################################################
      cat(DGP,"n=",n,"m=",m,"\n")
      ## #######################################################################
      
      BiasSq_mat <- matrix(NA, nrow = R, ncol = 6)
      Var_mat    <- matrix(NA, nrow = R, ncol = 6) 
      
      if(DGP=="DGP1"){nRegGrid  <-  11}else{nRegGrid  <-  51}
      
      ##
      for(r in 1:R){ # r <- 1

        ## ##################################################################
        cat("r/R=",r,"/",R,"\n")
        ## ##################################################################
        
        ## Generate partially observed *target* functions to be reconstructed
        while(TRUE){
          SimDat   <- ReconstPoFD::simuldata(n = 1, m = m, a = a, b = b, DGP=DGP, nRegGrid = nRegGrid)
          if(!all(range(c(na.omit(SimDat$U_mat[,1])))==c(0,1))){
            # takes only a partially observed function as a target function
            break
          }
        }
        ##
        Y_target_true_mat  <- SimDat[['Y_true_mat']]
        U_target_true_mat  <- SimDat[['U_true_mat']]
        ##
        Y_target_mat       <- SimDat[['Y_mat']]
        U_target_mat       <- SimDat[['U_mat']]
        Y_target_list      <- SimDat[['Y_list']]
        U_target_list      <- SimDat[['U_list']]
        ##
        A_target_vec       <- SimDat[['A_vec']]
        B_target_vec       <- SimDat[['B_vec']]
        ##
        missings_of_target <- c(U_target_true_mat[,1] < A_target_vec[1] | U_target_true_mat[,1] > B_target_vec[1])
        ##
        Y_EYes_AYes        <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_EYes_ANo         <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_EYes_AYes_CES    <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_EYes_ANo_CES     <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_PACE             <- matrix(NA, nrow = nRegGrid, ncol = B)
        ##
        Y_ENo_CG_AYes      <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_ENo_CG_ANo       <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_Kraus            <- matrix(NA, nrow = nRegGrid, ncol = B)
        #Y_PACE_ENo_CG      <- matrix(NA, nrow = nRegGrid, ncol = B)
        ##
        ## #######################################################################
        for(bb in 1:B){ # bb <- 1
          ##
          while(TRUE){
            SimDat       <- ReconstPoFD::simuldata(n = n-1, m = m, a = a, b = b, DGP=DGP, nRegGrid = nRegGrid)
            NonNA_fcts   <- apply(SimDat[['Y_mat']],2,function(x)!any(is.na(x)))
            n_NonNA_fcts <- length(NonNA_fcts[NonNA_fcts==TRUE])
            ##
            if(n_NonNA_fcts >= floor(n*0.2)){
              # guarantees that the samples contains at least 20% fully observed functions
              break
            }
          }
          ##
          Y_mat  <- cbind(Y_target_mat, SimDat[['Y_mat']])
          U_mat  <- cbind(U_target_mat, SimDat[['U_mat']])
          ##
          Y_list <- c(Y_target_list,    SimDat[['Y_list']])
          U_list <- c(U_target_list,    SimDat[['U_list']])
          ##
          if(DGP=="DGP1"){
            ##
            ## Reconstruction operator for noisy fragments with alignment at pre-smoothed fragment
            result_EYes_AYes <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                                   Lu           = U_list,
                                                                   method       = 'Error>0_AlignYES',
                                                                   reconst_fcts = 1,
                                                                   nRegGrid     = nRegGrid, 
                                                                   maxbins      = maxbins)
            Y_EYes_AYes[,bb] <- matrix(unlist(result_EYes_AYes[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
            ##
            ## Reconstruction operator for noisy fragments without alignment
            result_EYes_ANo <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                                  Lu           = U_list,
                                                                  method       = 'Error>=0_AlignNO',
                                                                  reconst_fcts = 1,
                                                                  nRegGrid     = nRegGrid, 
                                                                  maxbins      = maxbins)
            Y_EYes_ANo[,bb] <- matrix(unlist(result_EYes_ANo[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
            ##
            ## Reconstruction operator for noisy fragments with alignment at pre-smoothed fragment
            result_EYes_AYes_CES <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list,
                                                                       Lu           = U_list,
                                                                       method       = 'Error>0_AlignYES_CEscores',
                                                                       reconst_fcts = 1,
                                                                       nRegGrid     = nRegGrid,
                                                                       maxbins      = maxbins)
            Y_EYes_AYes_CES[,bb] <- matrix(unlist(result_EYes_AYes_CES[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1)
            # ##
            # ## Reconstruction operator for noisy fragments without alignment
            result_EYes_ANo_CES <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list,
                                                                      Lu           = U_list,
                                                                      method       = 'Error>0_AlignNO_CEscores',
                                                                      reconst_fcts = 1,
                                                                      nRegGrid     = nRegGrid,
                                                                      maxbins      = maxbins)
            Y_EYes_ANo_CES[,bb] <- matrix(unlist(result_EYes_ANo_CES[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1)
            ## 
            ##################
            ## PACE of Yao, Mueller, Wang (2005, JASA)
            # result_PACE           <- fdapace::FPCA(Ly = Y_list, Lt = U_list, optns = list(
            #   "dataType"="Sparse", "methodMuCovEst"="smooth", "error"=TRUE, "methodSelectK"="AIC", "nRegGrid"=nRegGrid))
            # Y_PACE[,bb]        <- t(fitted(result_PACE))[,1]
            ##################
            result_PACE <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                              Lu           = U_list,
                                                              method       = 'PACE',
                                                              reconst_fcts = 1,
                                                              nRegGrid     = nRegGrid, 
                                                              maxbins      = maxbins)
            Y_PACE[,bb] <- matrix(unlist(result_PACE[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
          }
          if(any(DGP==c("DGP2", "DGP3"))){
            ##
            ## Reconstruction operator for fully observed fragments with alignment at fully observed fragment
            result_ENo_CG_AYes <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                                     Lu           = U_list,
                                                                     method       = 'Error=0_AlignYES_CommonGrid',
                                                                     reconst_fcts = 1, 
                                                                     nRegGrid     = nRegGrid)
            Y_ENo_CG_AYes[,bb] <- matrix(unlist(result_ENo_CG_AYes[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
            ##
            ## Reconstruction operator for fully observed fragments without alignment
            result_ENo_CG_ANo  <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                                     Lu           = U_list,
                                                                     method       = 'Error>=0_AlignNO',
                                                                     reconst_fcts = 1, 
                                                                     nRegGrid     = nRegGrid)
            Y_ENo_CG_ANo[,bb]  <- matrix(unlist(result_ENo_CG_ANo[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
            ##
            ## Reconstruction Operator of Kraus (2015, JRSSB)
            result_Kraus       <- ReconstPoFD::reconstructKraus(X_mat = Y_mat, reconst_fcts = 1)
            Y_Kraus[,bb]       <- result_Kraus[['X_reconst_mat']]
          }
          ## ##################################################################
          if(bb %% 25 == 0) cat("bb/B=",bb,"/",B,"\n")
          ## ##################################################################
        } ## End of B-loop
        ##
        slct_M  <- missings_of_target
        ##
        plotting <- FALSE
        # plotting <- TRUE
        if(plotting){
          for(i in 1:min(B,10)){
            if(DGP=="DGP1"){
              par(mfrow=c(2,3))
              ##
              plot(Y_EYes_AYes[,i], type="b", ylim=range(Y_EYes_AYes[,i],Y_target_true_mat[slct_M,1]),main="Y_EYes_AYes")
              lines(Y_target_true_mat[,1]); points(y=Y_EYes_AYes[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_EYes_ANo[,i], type="b", ylim=range(Y_EYes_ANo[,i],Y_target_true_mat[slct_M,1]),main="Y_EYes_ANo")
              lines(Y_target_true_mat[,1]); points(y=Y_EYes_ANo[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              # plot(Y_EYes_AYes_CES[,i], type="b", ylim=range(Y_EYes_AYes_CES[,i],Y_target_true_mat[slct_M,1]),main="Y_EYes_AYes_CES")
              # lines(Y_target_true_mat[,1]); points(y=Y_EYes_AYes_CES[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              # ##
              # plot(Y_EYes_ANo_CES[,i], type="b", ylim=range(Y_EYes_ANo_CES[,i],Y_target_true_mat[slct_M,1]),main="Y_EYes_ANo_CES")
              # lines(Y_target_true_mat[,1]); points(y=Y_EYes_ANo_CES[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_PACE[,i], type="b", ylim=range(Y_PACE[,i],Y_target_true_mat[slct_M,1]),main="Y_PACE")
              lines(Y_target_true_mat[,1]); points(y=Y_PACE[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
            }
            ##
            if(any(DGP==c('DGP2','DGP3'))){
              par(mfrow=c(1,3))
              plot(Y_ENo_CG_AYes[,i], type="b", ylim=range(Y_ENo_CG_AYes[,i],Y_target_true_mat[slct_M,1]),main="Y_ENo_CG_AYes")
              lines(Y_target_true_mat[,1]); points(y=Y_ENo_CG_AYes[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_ENo_CG_ANo[,i], type="b", ylim=range(Y_ENo_CG_ANo[,i],Y_target_true_mat[slct_M,1]),main="Y_ENo_CG_ANo")
              lines(Y_target_true_mat[,1]); points(y=Y_ENo_CG_ANo[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_Kraus[,i], type="b", ylim=range(Y_Kraus[,i],Y_target_true_mat[slct_M,1]),main="Kraus")
              lines(Y_target_true_mat[,1]); points(y=Y_Kraus[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
            }
            par(mfrow=c(1,1))
            Sys.sleep(1)
          }
        }
        ##
        cut <- 0 # (if cut=0, no winsorization)
        ##
        Y_EYes_AYes_Mean     <- apply(Y_EYes_AYes[    slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        Y_EYes_ANo_Mean      <- apply(Y_EYes_ANo[     slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        Y_EYes_AYes_CES_Mean <- apply(Y_EYes_AYes_CES[slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        Y_EYes_ANo_CES_Mean  <- apply(Y_EYes_ANo_CES[ slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        Y_PACE_Mean          <- apply(Y_PACE[         slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        ##
        Y_ENo_CG_AYes_Mean   <- apply(Y_ENo_CG_AYes[  slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        Y_ENo_CG_ANo_Mean    <- apply(Y_ENo_CG_ANo[   slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        Y_Kraus_Mean         <- apply(Y_Kraus[        slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
        ##
        ## ## \int_{Missing} (E(X(t))-X(t))^2 dt ('Integrated squared bias')
        Y_EYes_AYes_Int_BiasSq     <- sum( (Y_EYes_AYes_Mean     - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        Y_EYes_ANo_Int_BiasSq      <- sum( (Y_EYes_ANo_Mean      - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        Y_EYes_AYes_CES_Int_BiasSq <- sum( (Y_EYes_AYes_CES_Mean - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        Y_EYes_ANo_CES_Int_BiasSq  <- sum( (Y_EYes_ANo_CES_Mean  - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        Y_PACE_Int_BiasSq          <- sum( (Y_PACE_Mean          - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        ##
        Y_ENo_CG_AYes_Int_BiasSq   <- sum( (Y_ENo_CG_AYes_Mean   - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        Y_ENo_CG_ANo_Int_BiasSq    <- sum( (Y_ENo_CG_ANo_Mean    - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        Y_Kraus_Int_BiasSq         <- sum( (Y_Kraus_Mean         - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        
        ## \int_{Missing} Var(X(t)) dt ('Integrated variance')
        Y_EYes_AYes_Int_Var     <- sum(apply(Y_EYes_AYes[    slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        Y_EYes_ANo_Int_Var      <- sum(apply(Y_EYes_ANo[     slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        Y_EYes_AYes_CES_Int_Var <- sum(apply(Y_EYes_AYes_CES[slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        Y_EYes_ANo_CES_Int_Var  <- sum(apply(Y_EYes_ANo_CES[ slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        Y_PACE_Int_Var          <- sum(apply(Y_PACE[         slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        ##
        Y_ENo_CG_AYes_Int_Var   <- sum(apply(Y_ENo_CG_AYes[  slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        Y_ENo_CG_ANo_Int_Var    <- sum(apply(Y_ENo_CG_ANo[   slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        Y_Kraus_Int_Var         <- sum(apply(Y_Kraus[        slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        ##
        BiasSq_mat[r,]          <- c(Y_EYes_AYes_Int_BiasSq,   Y_EYes_ANo_Int_BiasSq,   #Y_EYes_AYes_CES_Int_BiasSq, Y_EYes_ANo_CES_Int_BiasSq, 
                                     Y_PACE_Int_BiasSq,
                                     Y_ENo_CG_AYes_Int_BiasSq, Y_ENo_CG_ANo_Int_BiasSq, Y_Kraus_Int_BiasSq)
        Var_mat[r,]             <- c(Y_EYes_AYes_Int_Var,      Y_EYes_ANo_Int_Var,      #Y_EYes_AYes_CES_Int_Var,    Y_EYes_ANo_CES_Int_Var,    
                                     Y_PACE_Int_Var,
                                     Y_ENo_CG_AYes_Int_Var,    Y_ENo_CG_ANo_Int_Var,    Y_Kraus_Int_Var)
      }# end of r-loop
      ## 
      ## Tacking the mean over all target functions
      BiasSq_vec              <- colMeans(BiasSq_mat)
      Var_vec                 <- colMeans(Var_mat)
      names(BiasSq_vec)       <- c("EYes_AYes", "EYes_ANo", "EYes_AYes_CES", "EYes_ANo_CES", 
                                   "PACE", "ENo_CG_AYes", "ENo_CG_ANo", "Kraus")
      names(Var_vec)          <- c("EYes_AYes", "EYes_ANo", "EYes_AYes_CES", "EYes_ANo_CES", 
                                   "PACE", "ENo_CG_AYes", "ENo_CG_ANo", "Kraus")
      ##
      ## c(na.omit(BiasSq_vec))
      ## c(na.omit(Var_vec))
      ## sort(c(c(na.omit(BiasSq_vec)) + c(na.omit(Var_vec))))
      ##
      ## Save results:
      if(DGP=='DGP1'){
        save(BiasSq_vec, Var_vec, file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
      }
      if(any(DGP==c('DGP2','DGP3'))){
        save(BiasSq_vec, Var_vec, file = paste0(DGP,"_n",n,"_simResults.RData"))
      }
    }
  }
}
##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------


DGP <- c('DGP1','DGP2','DGP3')[1]
m   <- c(15,  30)[2] 
n   <- c(50, 100)[2] 

## Load results:
if(any(DGP==c('DGP1'))){
  load(file = paste0(DGP,"_n",n,"_m",m,"_simResults.RData"))
}
if(any(DGP==c('DGP2','DGP3'))){
  load(file = paste0(DGP,"_n",n,"_simResults.RData"))
}
##
col_slct         <- !is.na(BiasSq_vec)
BiasSq_vec       <- BiasSq_vec[col_slct]
Var_vec          <- Var_vec[col_slct]
MSE_vec          <- BiasSq_vec + Var_vec 
##
MSE_norm_vec     <- c(MSE_vec)    / max(MSE_vec)
BiasSq_norm_vec  <- c(BiasSq_vec) / max(MSE_vec)
Var_norm_vec     <- c(Var_vec)    / max(MSE_vec)
##
par(mfrow=c(1,3))
barplot(MSE_norm_vec, main="",    names.arg = names(MSE_norm_vec), ylim = c(0,1))
mtext(text = paste0("MSRE (", round(max(MSE_vec),2),")"), side = 3, line = 1)
barplot(BiasSq_norm_vec, main="", names.arg = names(MSE_norm_vec), ylim = c(0,1))
mtext(text = "Squared Bias", side = 3, line = 1)
barplot(Var_norm_vec, main="",    names.arg = names(MSE_norm_vec), ylim = c(0,1))
mtext(text = "Variance", side = 3, line = 1)
par(mfrow=c(1,1))   


# Interpretation: 
# 1. High variance in Kraus is due to (i) connection-points and (ii) non-smooth covariance estimate
# 2. CE-Scores (PACE-procedure) are extremely bad when estimated without a noise component as the noise component 
# has a ridge-type regularization-effect. 



