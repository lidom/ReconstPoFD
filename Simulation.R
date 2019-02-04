## ################################################
## Installing R-package 'ReconstPoFD' from GitHub
## ################################################
# install.packages("devtools")
# library("devtools")
# install_github("lidom/ReconstPoFD/ReconstPoFD")

## Load R-package 'ReconstPoFD'
library("ReconstPoFD") 

## R-packages used for constructing the LaTeX-outputs:
library("tidyverse")
library("xtable")

## Location to store the simulation results:
your_path <- "your/location/here"

## #######################################
## Number of target functions 
R            <-  50
## Number of MC-Repetitions per target function
B            <- 100
## #######################################
## Lower-value of the total domain 
a            <-  0             
## Upper-value of the total domain
b            <-  1
## #######################################
## Grid-size for evaluating the non-parametric 
## estimates of the mean and the covariance 
## For DGP1 (nRegGrid_vec[1]) we use a smaller 
## value to speed up the computations
nRegGrid_vec <- c(21, 51) 
## Number of binnings
## Smaller values speed up the computations
maxbins      <- 100   
## #######################################
## Sample sized
mm           <- c(15,  30) # number of points per function
nn           <- c(50, 100) # number of functions
## #######################################


## #######################################
## Simulate targets to be reconstructed 
## in each simulation-step
## #######################################

Y_target_true_mat_LIST  <- vector("list", c(4+4+2+2))
U_target_true_mat_LIST  <- vector("list", c(4+4+2+2))
Y_target_mat_LIST       <- vector("list", c(4+4+2+2))
U_target_mat_LIST       <- vector("list", c(4+4+2+2))
Y_target_list_LIST      <- vector("list", c(4+4+2+2))
U_target_list_LIST      <- vector("list", c(4+4+2+2))
A_target_vec_LIST       <- vector("list", c(4+4+2+2))
B_target_vec_LIST       <- vector("list", c(4+4+2+2))
##
list_names <- c("DGP1_n50_m15","DGP1_n50_m30","DGP1_n100_m15","DGP1_n100_m30", # DGP1
                "DGP2_n50_m15","DGP2_n50_m30","DGP2_n100_m15","DGP2_n100_m30", # DGP2
                "DGP3_n50","DGP3_n100",                                        # DGP3
                "DGP4_n50","DGP4_n100")                                        # DGP4
##
names(Y_target_true_mat_LIST)  <- list_names
names(U_target_true_mat_LIST)  <- list_names
names(Y_target_mat_LIST)       <- list_names
names(U_target_mat_LIST)       <- list_names
names(Y_target_list_LIST)      <- list_names
names(A_target_vec_LIST)       <- list_names
names(B_target_vec_LIST)       <- list_names


for(DGP in c('DGP1','DGP2','DGP3','DGP4')){# 
  ##
  set.seed(1)
  ##
  if(any(DGP==c("DGP1","DGP2"))){nRegGrid  <- nRegGrid_vec[1]}
  if(any(DGP==c("DGP3","DGP4"))){nRegGrid  <- nRegGrid_vec[2]}
  ##
  while(TRUE){
    # Simulate data
    SimDat   <- ReconstPoFD::simuldata(n = R, m = mm[2], a = a, b = b, DGP=DGP, nRegGrid = nRegGrid)
    if(!all(range(c(na.omit(SimDat$U_mat[,1])))==c(a,b))){
      # takes only a non-degenerated partially observed function as a target function
      break
    }
  }
  ##
  if(any(DGP==c("DGP1","DGP2"))){
    for(n in nn){
      for(m in mm){ # n <- nn[1]; m <- mm[1]
        if(m == mm[1]){slct <- sort(sample(1:mm[2], size = mm[1], replace = F))} # this selection guarantees that DGP1 and DGP2 have the same target functions, but only sampled at a different number of points m
        if(m == mm[2]){slct <- 1:mm[2]}
        Y_target_true_mat_LIST[[paste0(DGP,"_n",n,"_m",m)]]  <- SimDat[['Y_true_mat']]
        U_target_true_mat_LIST[[paste0(DGP,"_n",n,"_m",m)]]  <- SimDat[['U_true_mat']]
        Y_target_mat_LIST[[paste0(DGP,"_n",n,"_m",m)]]       <- SimDat[['Y_mat']][slct,]
        U_target_mat_LIST[[paste0(DGP,"_n",n,"_m",m)]]       <- SimDat[['U_mat']][slct,]
        for(r in 1:R){ # r <- 1
          Y_target_list_LIST[[paste0(DGP,"_n",n,"_m",m)]][[r]]  <- SimDat[['Y_list']][[r]][slct] # 'Y_mat' and 'Y_list' contain the same data!
          U_target_list_LIST[[paste0(DGP,"_n",n,"_m",m)]][[r]]  <- SimDat[['U_list']][[r]][slct] # 'U_mat' and 'U_list' contain the same data!
        }
        ##
        A_target_vec_LIST[[paste0(DGP,"_n",n,"_m",m)]]       <- SimDat[['A_vec']]
        B_target_vec_LIST[[paste0(DGP,"_n",n,"_m",m)]]       <- SimDat[['B_vec']]
      }
    }
  }
  if(any(DGP==c("DGP3","DGP4"))){
    for(n in nn){
      Y_target_true_mat_LIST[[paste0(DGP,"_n",n)]]  <- SimDat[['Y_true_mat']]
      U_target_true_mat_LIST[[paste0(DGP,"_n",n)]]  <- SimDat[['U_true_mat']]
      Y_target_mat_LIST[[paste0(DGP,"_n",n)]]       <- SimDat[['Y_mat']]
      U_target_mat_LIST[[paste0(DGP,"_n",n)]]       <- SimDat[['U_mat']]
      for(r in 1:R){ # r <- 1
        Y_target_list_LIST[[paste0(DGP,"_n",n)]][[r]]  <- SimDat[['Y_list']][[r]]
        U_target_list_LIST[[paste0(DGP,"_n",n)]][[r]]  <- SimDat[['U_list']][[r]]
      }
      ##
      A_target_vec_LIST[[paste0(DGP,"_n",n)]]       <- SimDat[['A_vec']]
      B_target_vec_LIST[[paste0(DGP,"_n",n)]]       <- SimDat[['B_vec']]
    }
  }
}




## ################################################################
##  Start Simulation ##############################################
## ################################################################
(Start.Time <- Sys.time())

for(DGP in c('DGP1','DGP2','DGP3','DGP4')){
  
  set.seed(2)
  
  for(n in nn){
    ##
    if(any(DGP==c("DGP1","DGP2"))){m_seq <- mm}
    if(any(DGP==c("DGP3","DGP4"))){m_seq <- NA}
    ##
    for(m in m_seq){
      
      ## DGP <- 'DGP1'; n <- nn[1]; m <- mm[1]; B <- 2; R <- 2
      
      ## #######################################################################
      cat(DGP,"n=",n,"m=",m,"\n")
      ## #######################################################################
      
      BiasSq_mat <- matrix(NA, nrow = R, ncol = 9)
      Var_mat    <- matrix(NA, nrow = R, ncol = 9) 
      
      if(any(DGP==c("DGP1","DGP2"))){  
        nRegGrid  <-  nRegGrid_vec[1]
        slct_name <- paste0(DGP,"_n",n,"_m",m)
        }
      if(any(DGP==c("DGP3","DGP4"))){
        nRegGrid  <-  nRegGrid_vec[2]
        slct_name <- paste0(DGP,"_n",n)
      }
      
      ##
      for(r in 1:R){ # r <- 1

        ## ##################################################################
        cat("r/R=",r,"/",R,"\n")
        ## ##################################################################
        
        ## Extract relevant targets
        Y_target_true_mat  <- Y_target_true_mat_LIST[[slct_name]][,r,drop=FALSE]
        U_target_true_mat  <- U_target_true_mat_LIST[[slct_name]][,r,drop=FALSE]
        Y_target_mat       <- Y_target_mat_LIST[[slct_name]][,r,drop=FALSE]
        U_target_mat       <- U_target_mat_LIST[[slct_name]][,r,drop=FALSE]
        Y_target_list      <- list("1"=Y_target_list_LIST[[slct_name]][[r]])
        U_target_list      <- list("1"=U_target_list_LIST[[slct_name]][[r]])
        missings_of_target <- c(U_target_true_mat[,1] < A_target_vec_LIST[[slct_name]][r] | 
                                U_target_true_mat[,1] > B_target_vec_LIST[[slct_name]][r])
        
        ## Declare variables for saving the results
        Y_EYes_AYes        <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_EYes_ANo         <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_EYes_AYes_CES    <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_EYes_ANo_CES     <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_PACE             <- matrix(NA, nrow = nRegGrid, ncol = B)
        ##
        Y_ENo_CG_AYes      <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_ENo_CG_ANo       <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_Kraus            <- matrix(NA, nrow = nRegGrid, ncol = B)
        Y_PACE_ENo_CG      <- matrix(NA, nrow = nRegGrid, ncol = B)
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
          if(any(DGP==c("DGP1", "DGP2"))){
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
            ##
            ## Reconstruction operator for noisy fragments without alignment
            result_EYes_ANo_CES <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list,
                                                                      Lu           = U_list,
                                                                      method       = 'Error>0_AlignNO_CEscores',
                                                                      reconst_fcts = 1,
                                                                      nRegGrid     = nRegGrid,
                                                                      maxbins      = maxbins)
            Y_EYes_ANo_CES[,bb] <- matrix(unlist(result_EYes_ANo_CES[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1)
            ## 
            ## PACE of Yao, Mueller, Wang (2005, JASA)
            result_PACE <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                              Lu           = U_list,
                                                              method       = 'PACE',
                                                              reconst_fcts = 1,
                                                              nRegGrid     = nRegGrid, 
                                                              maxbins      = maxbins)
            Y_PACE[,bb] <- matrix(unlist(result_PACE[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
          }
          if(any(DGP==c("DGP3", "DGP4"))){
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
            ##
            ## PACE for error==0
            result_PACE_E0  <- ReconstPoFD::reconstructKneipLiebl(Ly           = Y_list, 
                                                                  Lu           = U_list,
                                                                  method       = 'Error=0_PACE',
                                                                  reconst_fcts = 1, 
                                                                  nRegGrid     = nRegGrid)
            Y_PACE_ENo_CG[,bb]  <- matrix(unlist(result_PACE_E0[['Y_reconst_list']]), nrow = nRegGrid, ncol = 1) 
          }
          ## ##################################################################
          if(bb %% 50 == 0) cat("bb/B=",bb,"/",B,"\n")
          ## ##################################################################
        } ## End of B-loop
        ##
        slct_M  <- missings_of_target
        ##
        plotting <- FALSE
        # plotting <- TRUE
        if(plotting){
          for(i in 1:min(B,10)){
            if(any(DGP==c("DGP1", "DGP2"))){
              par(mfrow=c(2,3))
              ##
              plot(Y_EYes_AYes[,i], type="b", ylim=range(Y_EYes_AYes[,i],Y_target_true_mat[slct_M,1]),main="AYes")
              lines(Y_target_true_mat[,1]); points(y=Y_EYes_AYes[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_EYes_ANo[,i], type="b", ylim=range(Y_EYes_ANo[,i],Y_target_true_mat[slct_M,1]),main="ANo")
              lines(Y_target_true_mat[,1]); points(y=Y_EYes_ANo[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_EYes_AYes_CES[,i], type="b", ylim=range(Y_EYes_AYes_CES[,i],Y_target_true_mat[slct_M,1]),main="AYesCE")
              lines(Y_target_true_mat[,1]); points(y=Y_EYes_AYes_CES[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_EYes_ANo_CES[,i], type="b", ylim=range(Y_EYes_ANo_CES[,i],Y_target_true_mat[slct_M,1]),main="ANoCE")
              lines(Y_target_true_mat[,1]); points(y=Y_EYes_ANo_CES[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_PACE[,i], type="b", ylim=range(Y_PACE[,i],Y_target_true_mat[slct_M,1]),main="PACE")
              lines(Y_target_true_mat[,1]); points(y=Y_PACE[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
            }
            ##
            if(any(DGP==c('DGP3','DGP4'))){
              par(mfrow=c(2,2))
              plot(Y_ENo_CG_AYes[,i], type="b", ylim=range(Y_ENo_CG_AYes[,i],Y_target_true_mat[slct_M,1]),main="AYes")
              lines(Y_target_true_mat[,1]); points(y=Y_ENo_CG_AYes[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_ENo_CG_ANo[,i], type="b", ylim=range(Y_ENo_CG_ANo[,i],Y_target_true_mat[slct_M,1]),main="ANo")
              lines(Y_target_true_mat[,1]); points(y=Y_ENo_CG_ANo[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_Kraus[,i], type="b", ylim=range(Y_Kraus[,i],Y_target_true_mat[slct_M,1]),main="KRAUS")
              lines(Y_target_true_mat[,1]); points(y=Y_Kraus[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
              ##
              plot(Y_PACE_ENo_CG[,i], type="b", ylim=range(Y_PACE_ENo_CG[,i],Y_target_true_mat[slct_M,1]),main="PACE")
              lines(Y_target_true_mat[,1]); points(y=Y_PACE_ENo_CG[slct_M,i], x=c(1:nRegGrid)[slct_M], col="red")
            }
            par(mfrow=c(1,1))
            Sys.sleep(1)
          }
        }
        ##
        cut <- 0.005 # (if cut=0, no winsorization)
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
        Y_PACE_ENo_Mean      <- apply(Y_PACE_ENo_CG[  slct_M, ], 1, function(x){mean(ReconstPoFD:::winsorize_x(x,cut=cut))})
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
        Y_PACE_ENo_Int_BiasSq      <- sum( (Y_PACE_ENo_Mean      - Y_target_true_mat[slct_M, 1])^2 ) * (b-a)/nRegGrid
        
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
        Y_PACE_ENo_Int_Var      <- sum(apply(Y_PACE_ENo_CG[  slct_M, ], 1, function(x){var(ReconstPoFD:::winsorize_x(x,cut=cut))})) * (b-a)/nRegGrid
        ##
        BiasSq_mat[r,]          <- c(Y_EYes_AYes_Int_BiasSq,     Y_EYes_ANo_Int_BiasSq,   
                                     Y_EYes_AYes_CES_Int_BiasSq, Y_EYes_ANo_CES_Int_BiasSq, 
                                     Y_PACE_Int_BiasSq,          
                                     ##
                                     Y_ENo_CG_AYes_Int_BiasSq,   Y_ENo_CG_ANo_Int_BiasSq, Y_Kraus_Int_BiasSq,  Y_PACE_ENo_Int_BiasSq)
        ###
        Var_mat[r,]             <- c(Y_EYes_AYes_Int_Var,        Y_EYes_ANo_Int_Var,      
                                     Y_EYes_AYes_CES_Int_Var,    Y_EYes_ANo_CES_Int_Var,    
                                     Y_PACE_Int_Var,             
                                     ##
                                     Y_ENo_CG_AYes_Int_Var,      Y_ENo_CG_ANo_Int_Var,    Y_Kraus_Int_Var,     Y_PACE_ENo_Int_Var)
      }# end of r-loop
      ## 
      ## Tacking the mean over all target functions
      BiasSq_vec              <- colMeans(BiasSq_mat)
      Var_vec                 <- colMeans(Var_mat)
      ##
      name_vec                <- c("EYes_AYes",     "EYes_ANo", 
                                   "EYes_AYes_CES", "EYes_ANo_CES", 
                                   "PACE",          
                                   ##
                                   "ENo_CG_AYes",   "ENo_CG_ANo", "Kraus", "PACE_E0")
      names(BiasSq_vec) <- name_vec
      names(Var_vec)    <- name_vec
      ##
      ## Save results:
      if(any(DGP==c("DGP1", "DGP2"))){
        save(BiasSq_vec, Var_vec, file = paste0(your_path,DGP,"_n",n,"_m",m,"_simResults.RData"))
      }
      if(any(DGP==c('DGP3','DGP4'))){
        save(BiasSq_vec, Var_vec, file = paste0(your_path,DGP,"_n",n,"_simResults.RData"))
      }
    }
  }
}
##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------


## ######################################################
## Build LaTeX tables for the simulation results
## ######################################################

## --------------------------------
## DGP1 
## --------------------------------
DGP1_df_0 <- NULL

for(n in c(50, 100)){
  for(m in c(15,  30)){
    ## n <- 50; m <- 15
    load(file = paste0(your_path,"DGP1_n",n,"_m",m,"_simResults.RData"))
    ##
    col_slct         <- !is.na(BiasSq_vec)
    BiasSq_vec       <- BiasSq_vec[col_slct]
    Var_vec          <- Var_vec[col_slct]
    MSE_vec          <- BiasSq_vec + Var_vec 
    ##
    DGP1_df_0 <- rbind.data.frame(DGP1_df_0,
                                  tibble("DGP"   = as.factor("DGP1"),
                                         "Estim" = as.factor(c("EYes_AYes",     "EYes_ANo", 
                                                               "EYes_AYes_CES", "EYes_ANo_CES", 
                                                               "PACE")),
                                         "n"     = as.integer(n),
                                         "m"     = as.integer(m),
                                         "MSE"   = MSE_vec,
                                         "Bias2" = BiasSq_vec, 
                                         "Var"   = Var_vec)
    )
  }
}

DGP1_df <- DGP1_df_0 %>% 
  group_by(n, m) %>% 
  arrange(n, m, MSE) %>% 
  mutate("MSEr" = MSE / min(MSE)) %>% 
  ungroup() %>% 
  mutate(Estim = fct_recode(Estim, 
                            AYes   = "EYes_AYes", 
                            ANo    = "EYes_ANo", 
                            AYesCE = "EYes_AYes_CES",
                            ANoCE  = "EYes_ANo_CES",
                            PACE   = "PACE")) %>%
  select(DGP, n, m, Estim, MSEr, MSE, Bias2, Var)

print(DGP1_df, n = 20)

DGP1_df_LaTeX <- xtable(DGP1_df, digits=c(0,0,0,0,0,2,3,3,3))
print(DGP1_df_LaTeX, include.rownames = FALSE)


## --------------------------------
## DGP2 
## --------------------------------
DGP2_df_0 <- NULL

for(n in c(50, 100)){
  for(m in c(15,  30)){
    ## n <- 50; m <- 15
    load(file = paste0(your_path,"DGP2_n",n,"_m",m,"_simResults.RData"))
    ##
    col_slct         <- !is.na(BiasSq_vec)
    BiasSq_vec       <- BiasSq_vec[col_slct]
    Var_vec          <- Var_vec[col_slct]
    MSE_vec          <- BiasSq_vec + Var_vec 
    ##
    DGP2_df_0 <- rbind.data.frame(DGP2_df_0,
                                  tibble("DGP"   = as.factor("DGP2"),
                                         "Estim" = as.factor(c("EYes_AYes",     "EYes_ANo", 
                                                               "EYes_AYes_CES", "EYes_ANo_CES", 
                                                               "PACE")),
                                         "n"     = as.integer(n),
                                         "m"     = as.integer(m),
                                         "MSE"   = MSE_vec,
                                         "Bias2" = BiasSq_vec, 
                                         "Var"   = Var_vec)
    )
  }
}

DGP2_df <- DGP2_df_0 %>% 
  group_by(n, m) %>% 
  arrange(n, m, MSE) %>% 
  mutate("MSEr" = MSE / min(MSE)) %>% 
  ungroup() %>% 
  mutate(Estim = fct_recode(Estim, 
                            AYes   = "EYes_AYes", 
                            ANo    = "EYes_ANo", 
                            AYesCE = "EYes_AYes_CES",
                            ANoCE  = "EYes_ANo_CES",
                            PACE   = "PACE")) %>%
  select(DGP, n, m, Estim, MSEr, MSE, Bias2, Var)

print(DGP2_df, n = 20)

DGP2_df_LaTeX <- xtable(DGP2_df, digits=c(0,0,0,0,0,2,3,3,3))
print(DGP2_df_LaTeX, include.rownames = FALSE)



## --------------------------------
## DGP3 and DGP4
## --------------------------------
DGP3n4_df_0 <- NULL

for(DGP in c('DGP3', 'DGP4')){ 
  for(n in c(50, 100)){
    ##
    load(file = paste0(your_path,DGP,"_n",n,"_simResults.RData"))
    ##
    col_slct         <- !is.na(BiasSq_vec)
    BiasSq_vec       <- BiasSq_vec[col_slct]
    Var_vec          <- Var_vec[col_slct]
    MSE_vec          <- BiasSq_vec + Var_vec 
    ##
    DGP3n4_df_0 <- rbind.data.frame(DGP3n4_df_0,
                                    tibble("DGP"   = as.factor(DGP),
                                           "Estim" = as.factor(c("ENo_CG_AYes", "ENo_CG_ANo", "Kraus", "PACE_E0")),
                                           "n"     = as.integer(n),
                                           "MSE"   = MSE_vec,
                                           "Bias2" = BiasSq_vec, 
                                           "Var"   = Var_vec)
    )
  }
}

DGP3n4_df <- DGP3n4_df_0 %>% 
  mutate(Estim = fct_recode(Estim, 
                            AYes  = "ENo_CG_AYes", 
                            ANo   = "ENo_CG_ANo", 
                            KRAUS = "Kraus", 
                            PACE  = "PACE_E0")) %>% 
  group_by(DGP, n) %>% 
  arrange(DGP, n, MSE) %>% 
  mutate("MSEr" = MSE / min(MSE)) %>% 
  ungroup() %>% 
  select(DGP, n, Estim, MSEr, MSE, Bias2, Var)
DGP3n4_df

DGP2n3_df_LaTeX <- xtable(DGP2n3_df, digits=c(0,0,0,0,2,3,3,3))
print(DGP2n3_df_LaTeX, include.rownames = FALSE)


