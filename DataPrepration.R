# Prepare data needed in the DataAnalysis.R
#
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


dat <- read.csv("./data/exam_all.csv",header = T,stringsAsFactors = F,encoding = "UTF-8")



flocate=function(b){
  in.d=c(0,0);
  in.d[1]=which(colSums(b[3:4]==set)==2);
  in.d[2]=which(colSums(b[1:2]==set)==2);
  return(in.d)
}

set=matrix(0,2,9);set[1,1:3]=0;set[1,4:6]=1;set[1,7:9]=2;set[2,]=rep(c(0,1,2),3)

SubIdxFallible <- dat[,29]>=200
SubIdxTrue <- dat[,63]>=200
SubIdxTrue[is.na(SubIdxTrue)]=T
SubIdx <- SubIdxFallible&SubIdxTrue
dat <- dat[SubIdx,]
dim(dat)

NOBS <- nrow(dat)
NRES <- sum(!is.na(dat[,52+11]))

print("Total # of pairs:")
print(NOBS)
print("# of pairs measured by the true device:")
print(NRES)


colnames(dat)[c(29,34,63,68)]
head(dat[,c(29,34,63,68)])

NBeta_X <- 2

outable=matrix(0,NOBS,4);


unique(dat[,34])
unique(dat[,29])

Idx <- dat[,34]==0;outable[Idx,1]=0;outable[Idx,2]=0;
Idx <- dat[,34]==1;outable[Idx,1]=0;outable[Idx,2]=1;
Idx <- dat[,34]==2;outable[Idx,1]=0;outable[Idx,2]=2;
Idx <- dat[,34]==10;outable[Idx,1]=1;outable[Idx,2]=0;
Idx <- dat[,34]==11;outable[Idx,1]=1;outable[Idx,2]=1;
Idx <- dat[,34]==12;outable[Idx,1]=1;outable[Idx,2]=2;
Idx <- dat[,34]==20;outable[Idx,1]=2;outable[Idx,2]=0;
Idx <- dat[,34]==21;outable[Idx,1]=2;outable[Idx,2]=1;
Idx <- dat[,34]==22;outable[Idx,1]=2;outable[Idx,2]=2;

unique(dat[,68])

Idx <- dat[,68]=="NANA";outable[Idx,3]=-1;outable[Idx,4]=-1;
Idx <- dat[,68]=="00";outable[Idx,3]=0;outable[Idx,4]=0;
Idx <- dat[,68]=="01";outable[Idx,3]=0;outable[Idx,4]=1;
Idx <- dat[,68]=="02";outable[Idx,3]=0;outable[Idx,4]=2;
Idx <- dat[,68]=="10";outable[Idx,3]=1;outable[Idx,4]=0;
Idx <- dat[,68]=="11";outable[Idx,3]=1;outable[Idx,4]=1;
Idx <- dat[,68]=="12";outable[Idx,3]=1;outable[Idx,4]=2;
Idx <- dat[,68]=="20";outable[Idx,3]=2;outable[Idx,4]=0;
Idx <- dat[,68]=="21";outable[Idx,3]=2;outable[Idx,4]=1;
Idx <- dat[,68]=="22";outable[Idx,3]=2;outable[Idx,4]=2;



tfclean=outable[outable[,3]!=-1,];
tfall=outable[outable[,3]==-1,1:2];

qff=rep(0,9);

qtf=matrix(0,9,9);

tfre=rep(0,9);

for(iii in 1:NRES) {
  tind=flocate(tfclean[iii,]);

  qtf[tind[2],tind[1]]=qtf[tind[1],tind[2]]+1;	
}

for(ii in 1:(NOBS-NRES)){
  fin=which(colSums(tfall[ii,]==set)==2);
  qff[fin]=qff[fin]+1;
}



Observed <- matrix(rep(0, (9*NOBS)),nrow = 9,ncol = NOBS)
UnObserved <- matrix(rep(0, (9*NOBS)),nrow = 9,ncol = NOBS)
vec_both_true <- matrix(0,nrow = NOBS, ncol = 1)
vec_only_fal <- matrix(0,nrow = NOBS, ncol = 1)

R_O1 <- matrix(rep(0, (9*NRES)),nrow = 9,ncol = NRES)
R_O2 <- matrix(rep(0, (9*(NOBS - NRES))),nrow = 9,ncol = (NOBS - NRES))

vec_only_fal[(outable[,3]==-1), ] <- 1
vec_both_true[(outable[,3]!=-1), ] <-1

f <- outable[,1:2]
UnObserved_both_obs <- outable[,3:4]
for (ob_num in 1:NOBS) {
  Observed[which(colSums(set==f[ob_num,])==2),ob_num]=1
  UnObserved[which(colSums(set==UnObserved_both_obs[ob_num,])==2),ob_num]=1
  
}

Observed1 <- Observed[,as.logical(vec_only_fal)]
UnObserved1 <- UnObserved[,as.logical(vec_both_true)]




X_value <- matrix(NA,nrow=NOBS,ncol=NBeta_X)


dat[!is.na(dat[,29]),29] <- scale(dat[!is.na(dat[,29]),29])
dat[!is.na(dat[,63]),63] <- scale(dat[!is.na(dat[,63]),63])


X_value[vec_only_fal==1,] <- dat[vec_only_fal==1,29]
X_value[vec_both_true==1,] <- dat[vec_both_true==1,63]

X_valueTrue <- X_value[as.logical(vec_both_true),]
X_valueFallible <- X_value[as.logical(vec_only_fal),]





