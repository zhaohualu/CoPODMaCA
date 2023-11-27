# What'are needed for Data Analysis
# qtf: a K*K times K*K matrix, 
#          row is the classification by the fallible device, 
#          column by the true device
#
# qff: a K*K vector, counts of pairs classified by the fallible device only
#
# Observed1: a K*K times NOBS-NRES matrix
#         each column corresponds to a pair classified by the true device
#         with only one 1, all other zero. 
#         1 represents the category classified by the true device
#
# UnObserved1: a K*K times NRES matrix
#         each column with only one 1, all other zero.
#         each column corresponds to a pair classified ONLY by the fallible device
#         1 represents the category classified by the fallible device
#
# X_valueTrue: a NRES times NBeta_X matrix
#         each row is the covariate for one pair classified by the true device
#
# X_valueFallible: a NOBS-NRES times NBeta_X matrix
#         each row is the covariate for one pair classified ONLY by the fallible device
#

require(OpenMx)

############################################################################

Q         <- mxMatrix( type="Full", nrow=length(qff), ncol=length(qff),
                       free=FALSE,
                       values = qtf, name="Q")


F1        <- mxMatrix( type="Full", nrow=1, ncol=length(qff),
                       free=FALSE,
                       values = qff, name="F1")


one_vec        <- mxMatrix( type="Unit", nrow=length(qff), ncol=1,
                            free=FALSE, name="one_vec")


Falible_counts_acrosstrue <- mxAlgebra(expression =   Q %*% one_vec, name = "Falible_counts_acrosstrue")

NP <- mxAlgebra(expression =(Falible_counts_acrosstrue + t(F1))/Falible_counts_acrosstrue, name = "NP")

transition <- mxAlgebra(expression =(NP %*% t(one_vec)) * Q * (one_vec %*% t(1/(t(Q)%*%NP))), name = "transition")

XcovTrue         <- mxMatrix( type="Full", nrow=NRES, ncol=NBeta_X,
                              free=FALSE, 
                              values = X_valueTrue, name="XcovTrue")

XcovFallible         <- mxMatrix( type="Full", nrow=NOBS-NRES, ncol=NBeta_X,
                                  free=FALSE, 
                                  values = X_valueFallible, name="XcovFallible")

OneNobsRowVec <- mxMatrix( type="Unit", nrow=1, ncol=NOBS,
                           free=FALSE, 
                           name="OneNobsRowVec")
OneNobsRowVecTrue <- mxMatrix( type="Unit", nrow=1, ncol=NRES,
                               free=FALSE, 
                               name="OneNobsRowVecTrue")
OneNobsRowVecFallible <- mxMatrix( type="Unit", nrow=1, ncol=NOBS-NRES,
                                   free=FALSE, 
                                   name="OneNobsRowVecFallible")

I         <- mxMatrix( type="Unit", nrow=length(qff), ncol=1,
                       free=FALSE, 
                       name="I")


K         <- mxMatrix( type="Full", nrow=1, ncol=1,
                       free=FALSE,
                       values = c(9),
                       name="K")

U         <- mxMatrix( type="Full", nrow=1, ncol=1,
                       free=FALSE,
                       values = NOBS,
                       name="U")

currentAbscissaTrue <- mxMatrix(nrow=2, ncol=1, labels=c("abscissa1True","abscissa2True"), free=c(TRUE,TRUE), name="currentAbscissaTrue")
MeanStrTrue <- mxAlgebra(expression = t(MiuVec)%*%OneNobsRowVecTrue+ t(XcovTrue%*%Beta_X), name="abscissaTrue",dimnames=list(c('abscissa1True','abscissa2True'), NULL))
stuffTrue <- mxAlgebra(omxAllInt(S,cvectorize(currentAbscissaTrue),thres), name="stuffTrue")

currentAbscissaFallible <- mxMatrix(nrow=2, ncol=1, labels=c("abscissa1Fallible","abscissa2Fallible"), free=c(TRUE,TRUE), name="currentAbscissaFallible")
MeanStrFallible <- mxAlgebra(expression = t(MiuVec)%*%OneNobsRowVecFallible+ t(XcovFallible%*%Beta_X), name="abscissaFallible",dimnames=list(c('abscissa1Fallible','abscissa2Fallible'), NULL))
stuffFallible <- mxAlgebra(omxAllInt(S,cvectorize(currentAbscissaFallible),thres), name="stuffFallible")

R1 <- mxAlgebra(mxEvaluateOnGrid(stuffTrue, abscissaTrue), name="R1")
R2Fall <- mxAlgebra(mxEvaluateOnGrid(stuffFallible, abscissaFallible), name="R2Fall")
R2 <- mxAlgebra(expression = transition%*%R2Fall,name="R2")

O1 <- mxMatrix(values = UnObserved1,nrow = 9,ncol = NRES, name = "O1")
O2 <- mxMatrix(values = Observed1,nrow = 9,ncol = NOBS-NRES, name = "O2")

ds_cov_mle <- mxAlgebra(expression = -(sum(log(R1)*O1) + sum(log(R2)*O2)), name="log_Like")

funAl <- mxFitFunctionAlgebra("log_Like")

OneLocusModel <- mxModel("ds_cov_mle", one_vec , Falible_counts_acrosstrue, NP, transition, XcovTrue, XcovFallible, Beta_X, MiuVec, OneNobsRowVecTrue, OneNobsRowVecFallible, Q, F1, I, K, U, A, Y, S,MeanStrTrue, MeanStrFallible, currentAbscissaTrue,currentAbscissaFallible, stuffTrue,stuffFallible, R1, R2Fall, R2, O1, O2, thres, ds_cov_mle, funAl)

OneLocusFit <- mxRun(OneLocusModel)



