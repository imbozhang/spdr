################################################################################
# R Code for folded dimension reduction                                        
# Last updated: 2011/07/18                                                     
################################################################################

rm(list=ls())

################################################################################
# Preparing                                                                    
################################################################################

#-------------------------------------------------------------------------------
# if (method1==1) dr
# if (method2==1) sir
# if (method3==1) kdr
#-------------------------------------------------------------------------------

method1=1
method2=0
method3=1
pca=0

#-------------------------------------------------------------------------------
# Loading Packages
#-------------------------------------------------------------------------------

library(R.matlab)
library(R.utils)
library(Rlab)
library(mnormt)

#-------------------------------------------------------------------------------
# Call matlab from R
#-------------------------------------------------------------------------------

if (method3==1)
{
matlab <- Matlab()
isOpen <- open(matlab)
}                                                                                
################################################################################
# Simulation Model                                                             
################################################################################

#-------------------------------------------------------------------------------
# Model: Y= 0.1*exp(a'Xb+2) +0.5*e                                                       
#-------------------------------------------------------------------------------

gen.data1<-function(n,plef,prig,dlef,drig)
{
   alph.true <- matrix(c(1, -1, 0, 0, rep(0, plef-4)),plef,dlef)
   beta.true <- matrix(c(-1, 1, 0, 0, rep(0, prig-4)),prig,drig)
   b.true <- kronecker(beta.true,alph.true,)

   x=array(0,c(plef,prig,n))
   vecx = matrix(0,n,plef*prig)
   y = rep(0,n)

   for(i in 1:n){
   vecx[i,] <- rnorm(plef*prig)
   x[,,i] <- matrix(vecx[i,],plef,prig)
   y[i]=0.1*(t(alph.true) %*% x[,,i] %*% beta.true)^2+0.5*rnorm(1);
   }
   return(list(vecx=vecx,x=x,y=y,b.true=b.true,
   alph.true=alph.true,beta.true=beta.true))
}

#-------------------------------------------------------------------------------
# Model: Y=1 if 0.1(a'Xb)^2 > 1; Y=0 o.w.                                                       
#-------------------------------------------------------------------------------

gen.data2 <- function(n,plef,prig,dlef,drig)
{
   alph.true <- matrix(c(1, -1, 1, 0, rep(0, plef-4)),plef,dlef)
   beta.true <- matrix(c(-1, 1, -1, 0, rep(0, prig-4)),prig,drig)
   b.true <- kronecker(beta.true,alph.true)

   x=array(0,c(plef,prig,n))
   vecx = matrix(0,n,plef*prig)
   y = rep(0,n)

   for(i in 1:n){
   vecx[i,] <- rnorm(plef*prig)
   x[,,i] <- matrix(vecx[i,],plef,prig)
   crit <- 0.1*(t(alph.true) %*% x[,,i] %*% beta.true)^2+0.5
   if (crit>1)
   y[i]=1
   else
   y[i]=0
   }
   return(list(vecx=vecx,x=x,y=y,b.true=b.true,
   alph.true=alph.true,beta.true=beta.true))
}


#-------------------------------------------------------------------------------
# Model                                                       
#-------------------------------------------------------------------------------

gen.data3 <- function(n,plef,prig,dlef,drig)
{
   alph.true <- matrix(c(1, 1, 1, 0, rep(0, plef-4)),plef,dlef)
   beta.true <- matrix(c(-1, 1, 1, 0, rep(0, prig-4)),prig,drig)
   b.true <- kronecker(beta.true,alph.true)

   x=array(0,c(plef,prig,n))
   vecx = matrix(0,n,plef*prig)
   y = rep(0,n)

   for(i in 1:n){
   vecx[i,] <- rnorm(plef*prig)
   x[,,i] <- matrix(vecx[i,],plef,prig)
   crit <- t(alph.true) %*% x[,,i] %*% beta.true
   if (crit>0)
   y[i]=1
   else
   y[i]=0
   }
   return(list(vecx=vecx,x=x,y=y,b.true=b.true,
   alph.true=alph.true,beta.true=beta.true))
}

#-------------------------------------------------------------------------------
# Model: Example1 in Bing Li's paper                                                       
#-------------------------------------------------------------------------------

gen.data4 <- function(n,plef,prig,dlef,drig)
{
   alph.true <- matrix(c(1, 0, 0, 0, rep(0, plef-4),
                         0, 1, 0, 0, rep(0, plef-4)),plef,dlef)
   beta.true <- matrix(c(1, 0, 0, 0, rep(0, plef-4),
                         0, 1, 0, 0, rep(0, plef-4)),plef,dlef)
   b.true <- kronecker(beta.true,alph.true)
   # e1 <- matrix(c(1, 0, 0, 0, rep(0, plef-4)),plef,1)
   # e2 <- matrix(c(0, 1, 0, 0, rep(0, prig-4)),prig,1)
   # b.true <- cbind((kronecker(e1,e1),kronecker(e1,e2),
   # kronecker(e2,e1),kronecker(e2,e2))
   
   sigma2 = 0.1
   tao2 = 1.5
   mu = 1

   MU0 = matrix(0,plef*prig,1)
   SIGMA0 = diag(1,plef*prig)
   SIGMA0[2,2] = sigma2
   SIGMA0[prig+1,prig+1] =sigma2
   MU1 = matrix(0,plef*prig,1)
   MU1[1,1] = mu
   MU1[prig+2,1] =mu 
   SIGMA1 = diag(1,plef*prig)
   SIGMA1[2,2] = tao2
   SIGMA1[prig+1,prig+1] = tao2

   x=array(0,c(plef,prig,n))
   vecx = matrix(0,n,plef*prig)
   y = rep(0,n)

   for(i in 1:n){
   y[i]=rbern(1, 0.5)
   if (y[i]==0)
   vecx[i,] <- rmnorm(1,MU0,SIGMA0)
   else 
   vecx[i,] <- rmnorm(1,MU1,SIGMA1)
   x[,,i] <- matrix(vecx[i,],plef,prig)
   }
   return(list(vecx=vecx,x=x,y=y,b.true=b.true,
   alph.true=alph.true,beta.true=beta.true))
}

################################################################################
# Criterion                                                                    
################################################################################

#-------------------------------------------------------------------------------
# LZC distance                                                                           
#-------------------------------------------------------------------------------
dis = function(v1,v2){
p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
d <- sum((p1-p2)*(p1-p2))
      return(d)
}

#-------------------------------------------------------------------------------
# vector correlation and trace correlation between two spaces                                                                          
#-------------------------------------------------------------------------------
# normalize a vector
norm<-function(v)  
{ 
   sumv2<-sum(v^2)
   if(sumv2 == 0) sumv2<-1
   v/sqrt(sumv2)
}

# Gram-Schmidt orthonormalization
orthnorm<-function(X)
{
   X<-as.matrix(X)
   n<-nrow(X)
   p<-ncol(X)

   W<-NULL
   if(p > 1) {
      W<-cbind(W, X[,1])
      for(k in 2:p) {
         gw<-rep(0, n)
         for(i in 1:(k-1)) {
            gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
            gw<-gw + gki * W[,i]
         }
         W<-cbind(W, X[,k] - gw)
      }
   } else {
      W<-cbind(W, X[,1])
   }

   W<-apply(W, 2, norm)
   W
}

# vector correlation and trace correlation between two spaces 
eval.space<-function(A, B, orthnm=TRUE) 
{
   if(!is.matrix(A)) A<-as.matrix(A)
   if(!is.matrix(B)) B<-as.matrix(B)
   if(orthnm) { 
      A<-orthnorm(A)
      B<-orthnorm(B) 
   }

   mat<-t(B) %*% A %*% t(A) %*% B
   d<-eigen(mat)$values
   d<-(d+abs(d))/2
   q<-sqrt(prod(d))
   r<-sqrt(mean(d))
   ans<-list(q=q, r=r)
   return(ans)
}

#-------------------------------------------------------------------------------
# angle between two spaces
#-------------------------------------------------------------------------------
angles<-function(B1, B2)
{
   if(!is.matrix(B1)) B1<-as.matrix(B1)
   if(!is.matrix(B2)) B2<-as.matrix(B2)

   if(ncol(B1) >= ncol(B2)) {
      B<-B1; B.hat<-B2
   } else {
      B<-B2; B.hat<-B1
   }

   P1<-B %*% solve(t(B) %*% B) %*% t(B)
   if(ncol(B.hat) == 1) {
      nume<-as.vector(t(B.hat) %*% P1 %*% B.hat)
      deno<-as.vector(t(B.hat) %*% B.hat)
      ratio<-nume / deno
   } else {
      BtB<-t(B.hat) %*% B.hat
      ei<-eigen(BtB)
      BtB2<-ei$vectors %*% diag(1/sqrt(ei$values)) %*% t(ei$vectors)
      M<-BtB2 %*% t(B.hat) %*% P1 %*% B.hat %*% BtB2
      ratio<-abs(eigen(M)$values[nrow(M)])
   }
   ans<-acos(sqrt(ratio))/pi * 180
   if(ans > 90) ans<-180 - ans
   return(ans)
}

################################################################################
# Supporting function                                                          
################################################################################

#-------------------------------------------------------------------------------
# Subroutine 1:power of a matrix (ignore small evalues)
#-------------------------------------------------------------------------------
matpower = function(a,alpha){
if(ncol(a)==1) return(a^alpha) else
      small <- .00000000001
p1<-nrow(a)
eva<-eigen(a)$values
eve<-eigen(a)$vectors
eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
index<-(1:p1)[eva>small]
evai<-eva
evai[index]<-(eva[index])^(alpha)
ai<-eve%*%diag(evai)%*%t(eve)
return(ai)
}

#-------------------------------------------------------------------------------
# Subroutine 2:three powers of a matrix 
#-------------------------------------------------------------------------------
threematpower = function(a){
p1<-nrow(a)
eva<-eigen(a)$values
eve<-eigen(a)$vectors
eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
evai<- eva^(-1)
sigi <- eve%*%diag(evai)%*%t(eve)
evai <- eva^(-1/2)
signrt <- eve%*%diag(evai)%*%t(eve)
evai <- eva^(1/2)
sigrt <- eve%*%diag(evai)%*%t(eve)
ans<-list(sigi=sigi, signrt=signrt, sigrt=sigrt)
   return(ans)
}

#-------------------------------------------------------------------------------
# Subroutine 3:three powers of a matrix (Moore-Penrose-type) 
#-------------------------------------------------------------------------------
threematpower1 = function(a){
small = .00000000001
p1 = nrow(a)
evall = eigen(a)
eva = evall$values
eve = evall$vectors
eve = eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
index = (1:p1)[eva>small]
evai = eva
evai[index] = eva[index]^(-1)
sigi = eve%*%diag(evai)%*%t(eve)
evai[index] = eva[index]^(-1/2)
signrt = eve%*%diag(evai)%*%t(eve)
evai[index] = eva[index]^(1/2)
sigrt = eve%*%diag(evai)%*%t(eve)
ans<-list(sigi=sigi, signrt=signrt, sigrt=sigrt)
   return(ans)
}

#-------------------------------------------------------------------------------
# Subroutine 4:negative square root
#-------------------------------------------------------------------------------
nsrt<-function(a){
    evec<-eigen(a)$vectors
    eval<-eigen(a)$values
    b<-evec%*%diag(1/sqrt(eval))%*%t(evec)
    return(b)
}

#-------------------------------------------------------------------------------
# Subroutine 5:permutation matrix
#-------------------------------------------------------------------------------
perm = function(plef,prig,dlef,drig){
m = plef*prig*dlef*drig
gamma1 = matrix(0,m,m)
gamma2 = matrix(0,m,m)
counter = 0
for(i in 1:plef){
for(j in 1:dlef){
for(k in 1:prig){
for(ell in 1:drig){
counter = counter + 1
beta = matrix(0,prig,drig)
alph = matrix(0,plef,dlef)
alph[i,j] = 1
beta[k,ell] = 1
gamma1[,counter] = c(kronecker(beta,alph))
gamma2[,counter] = kronecker(c(beta),c(alph))
}}}}
return(gamma1%*%t(gamma2))
}

#-------------------------------------------------------------------------------
# Function for Folded SIR
#-------------------------------------------------------------------------------

foldsir = function(alph,beta,sig,evecx){
plef = nrow(alph)
prig = nrow(beta)
dlef = ncol(alph)
drig = ncol(beta)
nsli = nrow(evecx)
for(i1 in 1:10){
#print(c("iteration",i1))
##############   update c   #######################
mata = kronecker(beta,alph)
c = matrix(0,nsli,dlef*drig)
c[1,] = matpower(t(mata)%*%sig%*%mata,-1)%*%t(mata)%*%evecx[1,]
c[2,] = matpower(t(mata)%*%sig%*%mata,-1)%*%t(mata)%*%evecx[2,]
##############   update beta ######################
mf = array(0,c(dlef,drig,nsli))
for(i in 1:nsli){
mf[,,i] = matrix(c[i,],dlef,drig)}
iprig = diag(1,prig)
term1 = array(0,c(prig*plef,prig*drig,nsli))
for(i in 1:nsli){
term1[,,i] = kronecker(iprig,alph%*%matrix(mf[,,i],dlef,drig))}
evv = matrix(0,prig*drig,prig*drig)
evu = matrix(0,prig*drig,1)
for(i in 1:nsli){
evv = evv + prob[i]*t(term1[,,i])%*%sig%*%term1[,,i]  
evu = evu + prob[i]*t(term1[,,i])%*%evecx[i,]
} 
vecbetat = matpower(evv,-1)%*%evu
beta = t(matrix(vecbetat,drig,))
##############   update alpha #####################
term = array(0,c(prig*plef,plef*dlef,nsli))
iplef = diag(1,plef)
for(i in 1:nsli){
term[,,i] = kronecker(beta%*%t(matrix(mf[,,i],dlef,drig)),iplef)}
evv = matrix(0,plef*dlef,plef*dlef)
evu = matrix(0,plef*dlef,1)
for(i in 1:nsli){
evv = evv + prob[i]*t(term[,,i])%*%sig%*%term[,,i]
evu = evu + prob[i]*t(term[,,i])%*%evecx[i,]}
vecalph = matpower(evv,-1)%*%evu
alph = matrix(vecalph,plef,)
##############   standardize   #####################
#print(alpha)
alph = eigen(alph%*%t(alph))$vectors[,1:dlef]
beta = eigen(beta%*%t(beta))$vectors[,1:drig]
}
alph1 = matrix(c(alph),plef,dlef)
beta1 = matrix(c(beta),prig,drig)
ans<-list(alph1=alph1, beta1=beta1, c=c)
   return(ans)
}

################################################################################
# Simulation 												                                           
################################################################################

#-------------------------------------------------------------------------------
# Parameter Setting 
#-------------------------------------------------------------------------------

M=500                          # number of iterations of data sets
S=10                           # number of iterations of starting value 
T=200                          # number of iner iterations
epsilon = 0.1                  # for folded-dr only


# data related
n <- 500
plef <- 10
prig <- 10
dlef <- 2
drig <- 2
if (method3==1)
{
# set the same value of these parameters in Matlab
setVariable(matlab, pl=plef, pr=prig, dl=dlef, dr=drig, n=n)
}
set.seed(100)
if (method3==1)
{
evaluate(matlab, "rand('state',1);")
evaluate(matlab, "seed=floor(rand(1)*100);")
evaluate(matlab, "randn('state', seed);")
}
# result matrix
vc.alph = matrix(0,M,3)
vc.beta = matrix(0,M,3)
LZC = matrix(0,M,3)

#-------------------------------------------------------------------------------
# generate data
#-------------------------------------------------------------------------------
for(h in 1:M){
data <- gen.data4(n=n, plef=plef, prig=prig,dlef=dlef,drig=drig)
vecx <- data$vecx
x <- data$x
y <- data$y
b.true <- data$b.true
alph.true <- data$alph.true
beta.true <- data$beta.true
n0 = length((1:n)[y==0])
n1 = length((1:n)[y==1])
count = seq(0,0,length=n)
class = rep(0,n)	
kmat = perm(plef,prig,dlef,drig)
# store the original vecx
vecx0 = vecx 
xt = x                                  
yt = y

# sorting y
total <- cbind(y,vecx0)
newtotal <- total[order(y),]
yt <- newtotal[,1]
vecx <- newtotal[,2:(plef*prig+1)]
for (i in 1:n) {
xt[,,i] <- matrix(vecx[i,],plef,prig)
}  

# centering xt                   
xmu = apply(xt,c(1,2),mean)
xc = array(0,c(dim(xt)[1],dim(xt)[2],n))
for(i in 1:n){
xc[,,i] = xt[,,i]-xmu 
}
xt = xc

if (pca==1){

# PCA using E(XX^T) and E(X^TX)
xxt = matrix(0, nrow=dim(xt)[1], ncol=dim(xt)[1])
for (i in 1:n){
	xxt = xxt +xt[ , ,i] %*%t(xt[ , ,i])
	}
exxt =xxt/n
xtx = matrix(0, nrow=dim(xt)[2], ncol=dim(xt)[2])
for(i in 1:n){
	xtx = xtx+ (t(xt[ , ,i]) %*% xt[ , ,i])
	}
extx = xtx/n
pcaleft = eigen(exxt)$vectors[,(1:plef)]
pcaright = eigen(extx)$vectors[ ,(1:prig)]
xpca = array(0,c(plef, prig, n))
for(i in 1:n) {
	xpca[ , ,i] = t(pcaleft) %*% xt[ , ,i] %*% pcaright
	}
vecx = matrix(0 , nrow=n, ncol=plef*prig)
for(i in 1:n){
	vecx[i, ] <- c(xpca[ , ,i])
	}
	
}

sig = var(vecx) + epsilon*diag(1,plef*prig)
sigi = threematpower(sig)$sigi
signrt = threematpower(sig)$signrt
sigrt = threematpower(sig)$sigrt
vxs = vecx%*%signrt
nsli = 2
ind0 = (1:n)[yt==0]
ind1 = (1:n)[yt==1]
n0 = length(ind0)
n1 = length(ind1)
evecx= rbind(apply(vecx[1:n0,],2,mean),apply(vecx[(n0+1):(n0+n1),],2,mean))
prob = c(n0/n,n1/n)
m = plef*prig*dlef*drig*nsli*nsli

c = array(rnorm(m),c(plef*prig,dlef*drig,nsli,nsli))
eu00 = 2*diag(plef*prig)*(n0^2)-2*n0*(n0-1)*var(vxs[ind0,])
eu00 = eu00/(n0^2)
eu11 = 2*diag(plef*prig)*(n1^2)-2*n1*(n1-1)*var(vxs[ind1,])
eu11 = eu11/(n1^2)
eu01 = 2*diag(plef*prig)*(n0*n1)-
(n1*t(vxs[ind0,])%*%vxs[ind0,]-
apply(vxs[ind0,],2,sum)%*%t(apply(vxs[ind1,],2,sum))-
apply(vxs[ind1,],2,sum)%*%t(apply(vxs[ind0,],2,sum))+
n0*t(vxs[ind1,])%*%vxs[ind1,])
eu01 = eu01/(n0*n1)
eu10 = eu01

# starting value
alph = matrix(rnorm(plef*dlef),plef,dlef)
beta = matrix(rnorm(prig*drig),prig,drig)

# store the intital alph and beta
alph0= alph 
beta0 = beta

################################################################################
# folded-dr Method with fixed scale SIG                                                           
################################################################################
if (method1==1) {
for(i1 in 1:T){
###### step 1: update c_{Y, \tilde Y} ######
mata = sigrt%*%kronecker(beta,alph)
c = array(0,c(plef*prig,dlef*drig,nsli,nsli))
aaa = mata%*%solve(t(mata)%*%mata) 
c[,,1,1]= eu00%*%aaa
c[,,1,2]= eu01%*%aaa
c[,,2,1]=c[,,1,2]
c[,,2,2]= eu11%*%aaa
####### step 2: update beta ###############
tmp = kmat%*%kronecker(diag(prig*drig),c(alph))             
ecc = t(c[,,1,1])%*%c[,,1,1]*(n0^2 - n0)/(n^2-n)+
      2*t(c[,,1,2])%*%c[,,1,2]*(n0*n1)/(n^2-n)+
      t(c[,,2,2])%*%c[,,2,2]*(n1^2 - n1)/(n^2-n) 
evv = t(tmp)%*%kronecker(ecc,sig)%*%tmp
euc = eu00%*%c[,,1,1]*(n0^2 - n0)/(n^2-n)+
      2*eu01%*%c[,,1,2]*(n0*n1)/(n^2-n)+
      eu11%*%c[,,2,2]*(n1^2 - n1)/(n^2-n) 
evu = t(tmp)%*%c(sigrt%*%euc)
vecb = matpower(evv,-1)%*%evu
beta = matrix(vecb,prig,drig)
####### step 3: update alph ###############
tmp = kmat%*%kronecker(c(beta),diag(plef*dlef))             
evv = t(tmp)%*%kronecker(ecc,sig)%*%tmp
evu = t(tmp)%*%c(sigrt%*%euc)
veca = matpower(evv,-1)%*%evu
alph = matrix(veca,plef,dlef)   
alph = as.matrix(eigen(alph%*%t(alph))$vectors[,(1:dlef)])
beta = as.matrix(eigen(beta%*%t(beta))$vectors[,(1:drig)])
}
vc.alph[h,1] = eval.space(alph.true, alph, orthnm=TRUE)$q
vc.beta[h,1] = eval.space(beta.true, beta, orthnm=TRUE)$q
LZC[h,1] = dis(b.true,kronecker(beta,alph)) 
}

################################################################################
# folded-sir Method                                                            
################################################################################

if (method2==1)
{
sirresult <- foldsir(alph0,beta0,sig,evecx)
alph <- sirresult$alph
beta <- sirresult$beta
vc.alph[h,2] = eval.space(alph.true, alph ,orthnm=TRUE)$q
vc.beta[h,2] = eval.space(beta.true, beta ,orthnm=TRUE)$q
LZC[h,2] = dis(b.true,kronecker(beta,alph))
}

################################################################################
# folded-kdr Method with fixed scale parameter                                                          
################################################################################
if (method3==1)
{
mintr <- 1e+9;
setVariable(matlab, vecX=vecx0, Y=y)
evaluate(matlab, "vecX=(vecX-repmat(mean(vecX),n,1))./repmat(std(vecX),n,1);")
evaluate(matlab, "sgx0=2*(dl*dr)/(pl*pr)*MedianDist(vecX)^2;")
evaluate(matlab, "sgy0=2*MedianDist(Y)^2;")

for(j in 1:S){
alph0 = matrix(rnorm(plef*dlef),plef,dlef)
beta0 = matrix(rnorm(prig*drig),prig,drig)                      # starting value
evaluate(matlab, "[L,R,B,tr0]=FoldedKDR(vecX,Y,pl,dl,pr,dr,sgx0,sgy0);")
alph <- getVariable(matlab, "L")$L
beta <- getVariable(matlab, "R")$R
tr0 <- getVariable(matlab, "tr0")$tr0

if (tr0 < mintr)
{
mintr=tr0
minalph=alph
minbeta=beta
}
alph=minalph
beta=minbeta
}

vc.alph[h,3] = eval.space(alph.true, alph, orthnm=TRUE)$q
vc.beta[h,3] = eval.space(beta.true, beta, orthnm=TRUE)$q
LZC[h,3] = dis(b.true,kronecker(beta,alph))
}
print(h)
}     
#vc.median <- matrix(c(median(vc.alph[,1]),median(vc.alph[,2]),median(vc.alph[,3]), 
#median(vc.beta[,1]),median(vc.beta[,2]), median(vc.beta[,3])),3,2)
#colnames(vc.median) <- c("alph","beta")
#rownames(vc.median) <- c("Folded-DR","Folded-SIR","Folded-KDR")

#vc.median.mad <- matrix(c(mad(vc.alph[,1]),mad(vc.alph[,2]),mad(vc.alph[,3]), 
#mad(vc.beta[,1]),mad(vc.beta[,2]), mad(vc.beta[,3])),3,2)
#colnames(vc.median.mad) <- c("alph","beta")
#rownames(vc.median.mad) <- c("Folded-DR","Folded-SIR","Folded-KDR")

vc.mean <- matrix(c(mean(vc.alph[,1]),mean(vc.alph[,2]),mean(vc.alph[,3]), 
mean(vc.beta[,1]),mean(vc.beta[,2]), mean(vc.beta[,3])),3,2)
colnames(vc.mean) <- c("alph","beta")
rownames(vc.mean) <- c("Folded-DR","Folded-SIR","Folded-KDR")

vc.mean.se <- matrix(c(sd(vc.alph[,1]),sd(vc.alph[,2]),sd(vc.alph[,3]), 
sd(vc.beta[,1]),sd(vc.beta[,2]), sd(vc.beta[,3])),3,2)/sqrt(M)
colnames(vc.mean.se) <- c("alph","beta")
rownames(vc.mean.se) <- c("Folded-DR","Folded-SIR","Folded-KDR")

#LZC.median <- matrix(c(median(LZC[,1]),median(LZC[,2]),median(LZC[,3])))
#colnames(LZC.median) <- "b"
#rownames(LZC.median) <- c("Folded-DR","Folded-SIR","Folded-KDR")

#LZC.median.mad <- matrix(c(mad(LZC[,1]),mad(LZC[,2]),mad(LZC[,3])))
#colnames(LZC.median.mad) <- "b"
#rownames(LZC.median.mad) <- c("Folded-DR","Folded-SIR","Folded-KDR")

LZC.mean <- matrix(c(mean(LZC[,1]),mean(LZC[,2]),mean(LZC[,3])))
colnames(LZC.mean) <- "b"
rownames(LZC.mean) <- c("Folded-DR","Folded-SIR","Folded-KDR")

LZC.mean.se <- matrix(c(sd(LZC[,1]),sd(LZC[,2]),sd(LZC[,3])))/sqrt(M)
colnames(LZC.mean.se) <- "b"
rownames(LZC.mean.se) <- c("Folded-DR","Folded-SIR","Folded-KDR")

#vc.median
#vc.median.mad
#vc.mean
#vc.mean.se
#LZC.median
#LZC.median.mad
#LZC.mean
#LZC.mean.se

cbind(vc.mean,LZC.mean)

