
# Number of categories
NK <- 3 

# Number of covariates for the treatment and reference variables in total
NBeta_X <- ncol(X_valueTrue) 


NRES <- ncol(UnObserved1) # Number of pairs measured ONLY by fallible device 
NOBS <- ncol(Observed1) + NRES # Number of pairs


################################################################################
# Only change the values argument in the below OpenMx matrices.

# Coefficient for the covariates
Beta_X         <- mxMatrix( type="Full", nrow= NBeta_X, ncol=2,
                            free=c(TRUE, FALSE, FALSE,TRUE),
                            values = matrix(c(1,0,0,1),nrow=NBeta_X, ncol=2)*0.5,
                            name="Beta_X")

# Mean of the underlying latent variables for the treatment and reference variable
MiuVec <- mxMatrix(type="Full", nrow=1, ncol=2,
                   free=c(FALSE,TRUE), # The mean of reference variable is fixed at 0
                   values=c(0, 0.5), # The mean of reference and treatment latent variables
                   lbound =c(NA,-2), # The lower bound
                   ubound =c(NA,2), # The upper bound
                   labels = c(NA, "MiuT"), name="MiuVec")

# Matrix transforming the thresholds parameters to the real thresholds for the underlying latent variables 
A         <- mxMatrix( type="Lower", nrow=NK+1, ncol=NK+1,
                       free=FALSE, values=c(1,0,0,0,
                                            0,1,1,0,
                                            0,0,1,0,
                                            0,0,0,1), name="A")
# thresholds parameters before transformation
Y         <- mxMatrix( type="Full", nrow=NK+1, ncol=2, # Each column for reference and treatment variables, respectively
                       free=c(FALSE,TRUE,TRUE,FALSE,
                              FALSE,TRUE,TRUE,FALSE),
                       labels = c(NA,"TH1", "TH2", NA,
                                  NA,"TH1", "TH2",NA), # Matrix arrange by column, equivalently transpose this matrix
                       lbound = c(NA,-10,0.000001,NA,
                                  NA,-10,0.000001,NA),
                       ubound = c(NA,10,20,NA,
                                  NA,10,20,NA),
                       values=c(-1000,-0.25, 0.25, 1000,
                                -1000,-0.25, 0.25, 1000), name="Y")

# Covariance matrix for the underlying latent variables
S         <- mxMatrix( type="Symm", nrow=2, ncol=2,
                       free=c(FALSE,TRUE,TRUE,TRUE),
                       values=c(1, 0.2 , 0.2, 1),
                       lbound = c(NA,NA,NA,0.000001),
                       labels = c(NA,"COV21", "COV21", "COV22"),
                       name="S")

thres <- mxAlgebra(expression =   (A%*%Y), name="thres")