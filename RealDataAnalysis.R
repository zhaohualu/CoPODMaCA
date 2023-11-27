library(MASS)
require(mvtnorm) 
library(OpenMx)


source("DataPrepration.R")


save(file="DataInput.RData",list=c("qtf","qff","Observed1","UnObserved1","X_valueTrue","X_valueFallible"))
rm(list=ls())
load("DataInput.RData")


source("SetParameters.R")


tryCatch({
  source("DataAnalysis.R")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})



# Print results: 

print("Estimates:")
OneLocusFit$output$estimate[1:7]
print("SEs:")
OneLocusFit$output$standardErrors

print("Estimated Correlation Matrix:")
cov2cor(matrix(c(1,0.248,0.248,0.205),2))

print("P-values:")
pnorm(abs(OneLocusFit$output$estimate[1:7]/OneLocusFit$output$standardErrors),lower.tail = F)*2

print("Confidence Intervals:")
cbind(OneLocusFit$output$estimate[1:7]-1.96*OneLocusFit$output$standardErrors,
  OneLocusFit$output$estimate[1:7]+1.96*OneLocusFit$output$standardErrors)

# Save the results
save.image(paste(c("DS_COV_",".RData"),collapse=""))

sessionInfo()

