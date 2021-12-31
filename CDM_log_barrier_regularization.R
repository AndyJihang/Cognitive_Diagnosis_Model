library(CDM)
ecpe<-data.ecpe$data[,-1]
alpha <- matrix(rbinom(15000, 1, 0.5), nrow=5000)
q <- data.ecpe$q.matrix
tmatrix <- alpha%*%t(q)
attsum <- matrix(1, nrow=1, ncol=28)
for (i in 1:28){
  attsum[i] <- sum(q[i,])
}
  
reattsum <- matrix(rep(attsum, 5000), nrow=5000, byrow=T)
eta <- matrix(1, nrow=5000,ncol=28)
for(i in 1:5000){
  for(j in 1:28){
    if(tmatrix[i,j]==reattsum[i,j]){
      eta[i,j] <- 1
    }else{
      eta[i,j] <- 0
    }
  }
}
p <- matrix(1, nrow=5000,ncol=28)
for(i in 1:5000){
  for(j in 1:28){
    p[i,j] <- (0.95^eta[i,j])*(0.20^(1-eta[i,j]))
  }
}
## generating data ##
sim <- matrix(1,nrow=5000,ncol=28)
for(i in 1:5000){
  for(j in 1:28){
    sim[i,j] <- rbinom(1,1,p[i,j])
  }
}
  
#############################################################################################################
as.binary <-
  function(x){
    ans <- NULL
    while(any(x!=0)){
      ans <- cbind(x%%2,ans)
      x <- floor(x/2)
    }
    ans
  }
  
optim.fun <- function(X,Q,pi,a=0.5){
  N = nrow(X)
  J = ncol(X)
  K = ncol(Q)
  all.a = as.binary(0:(2^K-1)) #2^K x K #
  pi = c(1,exp(pi))
  pi = pi/sum(pi)
  Q = exp(Q)
  Q = Q/(1+Q)
  #print(pi)
  pR.a = exp(tcrossprod(1-all.a,log(1-Q))) #2^K x J #
  tmp = tcrossprod(X,log(pR.a)) + tcrossprod(1-X,log(1-pR.a)) # N x 2^K #
  tmp[,2^K] = ifelse(apply(X==1,1,all), 0, -Inf)
  pX = exp(tmp)%*%pi
  result = sum(log(pX)) - a*sum(log(Q)+log(1-Q))
  -2*result
}
  
  
findQ = function(y,ndim=5,a=0.5,trace=0,maxit=500){
  y = as.matrix(y)
  nitem = ncol(y)
  nsubj = nrow(y)
  init = c(rnorm(nitem*ndim),rep(0,2^ndim-1))
  tmp.out = optim(init, function(p){
    optim.fun(y, matrix(p[1:(nitem*ndim)],nitem,ndim), p[-(1:(nitem*ndim))],a=a)},
    control=list(maxit=maxit,trace=trace))
  final = tmp.out$par
  Q = matrix(final[1:(ndim*nitem)], nitem, ndim)
  Q = exp(Q)
  Q = Q/(1+Q)
  pi = c(1,exp(final[-(1:(ndim*nitem))]))
  pi = pi/sum(pi)
  result = list(Q=Q,pi=pi,deviance=tmp.out$value)
  result
}
  
a <- c(0,0.001,0.05,1,2,3,6,9,11,13)
R_1 <- findQ(y=ecpe,ndim=3,a=a[1],trace=0,maxit=10000)
R_2 <- findQ(y=ecpe,ndim=3,a=a[2],trace=0,maxit=10000)
R_3 <- findQ(y=ecpe,ndim=3,a=a[3],trace=0,maxit=10000)
R_4 <- findQ(y=ecpe,ndim=3,a=a[4],trace=0,maxit=10000)
R_5 <- findQ(y=ecpe,ndim=3,a=a[5],trace=0,maxit=10000)
R_6 <- findQ(y=ecpe,ndim=3,a=a[6],trace=0,maxit=10000)
R_7 <- findQ(y=ecpe,ndim=3,a=a[7],trace=0,maxit=10000)
R_8 <- findQ(y=ecpe,ndim=3,a=a[8],trace=0,maxit=10000)
R_9 <- findQ(y=ecpe,ndim=3,a=a[9],trace=0,maxit=10000)
R_10 <- findQ(y=ecpe,ndim=3,a=a[10],trace=0,maxit=10000)
  
##########################################################################################################
  
fn_perm <- function (n, r, v = 1:n){
  if (r == 1)
    matrix(v, n, 1)
  else if (n == 1)
    matrix(v, 1, r)
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], fn_perm(n - 1, r - 1, v[-i])))
    X
  }
}
  
Q2=q
  
discrepQ2=function(Q.EST){
  -sum(Q2*log(Q.EST)+(1-Q2)*log(1-Q.EST))
}

mindiscrepQ2=function(Q.EST){
  all_perm=fn_perm(ncol(Q.EST),ncol(Q.EST))
  all_result <- c()
  for(i in 1:nrow(all_perm)){
    newQ.EST <- Q.EST[,all_perm[i,]]
    one_result <- discrepQ2(newQ.EST)
    all_result <- c(all_result, one_result)
  }
  mindis <- min(all_result)
  index_min <- which(all_result == min(all_result), arr.ind = TRUE)
  min_perm <- all_perm[index_min,]
  FQ.EST=Q.EST[,min_perm]
  result=list(mindis=mindis,min_perm=min_perm,FQ.EST=FQ.EST)
  result
}
  
#disQ2.EST1=mindiscrepQ2(R_1$Q)
#disQ2.EST1=mindiscrepQ2(R_2$Q)
#disQ2.EST1=mindiscrepQ2(R_3$Q)
disQ2.EST1=mindiscrepQ2(R_4$Q)
#disQ2.EST1=mindiscrepQ2(R_5$Q)
#disQ2.EST1=mindiscrepQ2(R_6$Q)
#disQ2.EST1=mindiscrepQ2(R_7$Q)
#disQ2.EST1=mindiscrepQ2(R_8$Q)
#disQ2.EST1=mindiscrepQ2(R_9$Q)
#disQ2.EST1=mindiscrepQ2(R_10$Q)
#disQ2.EST1=mindiscrepQ2(R$Q)

mindisQ2.EST1=disQ2.EST1$mindis
FQ2.EST1=disQ2.EST1$FQ.EST
count <- function(a,b,q,Q){
 qnew=matrix(,nrow=nrow(q),ncol=ncol(q))
  for (i in 1:nrow(q)){
    for (j in 1:ncol(q)){
      if(q[i,j]<=a){qnew[i,j]=0}
      else {if (q[i,j]>=b){qnew[i,j]=1} else {qnew[i,j]=q[i,j]}}
    }
  }
dif=Q-qnew
difnew=matrix(,nrow=nrow(dif),ncol=ncol(dif))
for (i in 1:nrow(dif)){
    for (j in 1:ncol(dif)){
      if(dif[i,j]==0){difnew[i,j]=1}
      else {difnew[i,j]=0}
    }
  }
  rowsum=apply(difnew,1,sum)
  counts=sum(rowsum)
  total=ncol(Q)*nrow(Q)
prob=counts/total
result=list(cutoff=c(a,b),FinalQ=qnew,Correct_Counts=counts,
            Total=total,Correct_Prob=prob)
result
}
  
#count(0.5,0.5,q=R_1$Q,Q=q)
#count(0.5,0.5,q=R_2$Q,Q=q)
#count(0.5,0.5,q=R_3$Q,Q=q)
#count(0.5,0.5,q=R_4$Q,Q=q)
#count(0.5,0.5,q=R_5$Q,Q=q)
#count(0.5,0.5,q=R_6$Q,Q=q)
#count(0.5,0.5,q=R_7$Q,Q=q)
count(0.5,0.5,q=R_8$Q,Q=q)
#count(0.5,0.5,q=R_9$Q,Q=q)
#count(0.5,0.5,q=R_10$Q,Q=q)
#count(0.5,0.5,q=R$Q,Q=q)
mindiscrepQ2(R_8$Q)$min_perm