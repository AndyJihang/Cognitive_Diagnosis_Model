library(gtools);

prof <- function(K) permutations(n=2,r=K,v=c(0,1), repeats.allowed= T)[,K:1]

interact <- function(K){

  X <- t(prof(K));
  Y <- X+1;
  Z <- matrix(0, 2^K, 2^K);
  
  for(i in 1:(2^K)){
    Z[,i] = as.numeric(as.logical(apply(Y - X[,i],2,prod)));
  } 
  Z
} 

logit <- function(x){
  return(log(x/(1-x)));
}

soft_thre <- function(x, mu) sign(x)*pos(abs(x) - mu)

pos <- function(x){
  x[x<0] = 0;
  x;
}


theta.fun <- function(q.index, J, K, s.vec, g.vec){

  c.vec <- 1-s.vec;
  theta <- matrix(0, J, 2^K);
  theta[,1] <- logit(g.vec);
  for(i in 1:J){
    theta[i, q.index[i]] <- logit(c.vec[i]) - logit(g.vec[i]);
  }
  theta;
  
}

prob.fun <- function(theta, Z){
  temp <- theta %*% t(Z);
  exp(temp)/(1+exp(temp));
}

post.fun <- function(theta, pi, Z, data){ #Given theta and pi, compute the posterior distribution
  prob <- prob.fun(theta, Z);
  
  logpost <- t(t(data %*% log(prob) + (1-data) %*% log(1-prob)) + log(pi));
  
  logpost.max <- apply(logpost,1,max);
  temp = exp(logpost - logpost.max); #To avoid unstable result when take exp;
  
  sumtemp = rowSums(temp);
  
  temp/sumtemp;
} 

log.lik <- function(theta, pi, Z, data){ #Given theta and pi, compute the likelihood function

  prob <- prob.fun(theta, Z);
  temp <- exp(data %*% log(prob) + (1-data) %*% log(1-prob)) %*% pi;
  
  sum(log(rowSums(temp)));
} 

obj <- function(theta, pi, Z, data, lambda){
  N = nrow(data);
  -log.lik(theta, pi, Z, data)/N + sum(lambda %*% abs(theta[,-1]))
}

M.step.pi <- function(post, delta){ #delta is from a dirichlet prior to avoid hat.pi = 0.
  #This function gives the hat.pi in the M step.
  temp = colSums(post)
  (temp + 1) / sum(temp+1); 
}

obj.j <- function(post, theta.j, Z, data.j, lambda.j){
  N = length(data.j);
  temp = theta.j %*% t(Z);
  
  obj = data.j %*% post %*% t(temp)  - log(1+exp(temp)) %*% colSums(post)
  -obj/N + lambda.j * sum(abs(theta.j[-1]));
}

grad.j <- function(post, theta.j, Z, data.j){
  
  N = length(data.j);
  temp = theta.j %*% t(Z);
  
  - (data.j %*% post %*% Z - (colSums(post) * (exp(temp)/(1+exp(temp)))) %*% Z)/N;
}

M.step.theta.j <- function(post, theta.j, Z, data.j, lambda.j, step, totstep){ 
  #The optimization is done separately for each item. 
  #We use a Proximal gradient update.
  
  grad <- grad.j(post, theta.j, Z, data.j); 
  obj0 <- obj.j(post, theta.j, Z, data.j, lambda.j); 
  
  y <- theta.j - step * grad;
  theta.new <- soft_thre(y, lambda.j * step);
  
  obj1 <- obj.j(post, theta.new, Z, data.j, lambda.j); 
  
  z = 1;
  while(obj1 > obj0 & z < totstep){
    step = step * 0.5;
    y <- theta.j - step * grad;
    theta.new[1] <- y[1];
    theta.new[-1] <- soft_thre(y[-1], lambda.j * step);
    
    obj1 <- obj.j(post, theta.new, Z, data.j, lambda.j); 
    z = z+1;
  }
  theta.new;
}

M.step <- function(post, theta, Z, data, lambda, delta, step, totstep){
  N = nrow(data);
  J = ncol(data);
  
  hat.pi <- M.step.pi(post,delta);
  hat.theta <- theta;
  
  for(j in 1:J){
    hat.theta[j,] <- M.step.theta.j(post, theta[j,], Z, data[,j], lambda[j], step, totstep);
  }
  
  obj.value = obj(hat.theta, hat.pi, Z, data, lambda);
  
  list(hat.theta = hat.theta, hat.pi = hat.pi, obj.value = obj.value);
}

EM <- function(theta, pi, Z, data, lambda, delta = 1, step = 1, totstep = 20, tol = 1e-6){
  
  obj0 <- obj(theta, pi, Z, data, lambda);
  post = post.fun(theta, pi, Z, data);
  M.result = M.step(post, theta, Z, data, lambda, delta, step, totstep)
  
  obj1 <- M.result$obj.value;

  
  while(obj0-obj1 > tol){
    #print(obj1)
    obj0 <- obj1;
    theta <- M.result$hat.theta;
    pi <- M.result$hat.pi;
    post = post.fun(theta, pi, Z, data);
    M.result = M.step(post, theta, Z, data, lambda, delta, step, totstep)
    obj1 <- M.result$obj.value;
  }
  M.result;
}

Path <- function(theta, pi, Z, data, lambda, X, h.thre = 0.1, path.tot = 30, lambda.step = 0.01, delta = 1, step = 1, totstep = 20, tol = 1e-6){
  
  for(i in 1:path.tot){
    result <- EM(theta, pi, Z, data, lambda, delta, step, totstep, tol);
    ind.matr <- abs(result$hat.theta[,-1]) > h.thre; #a hard thresholding step;
    ind = rowSums(ind.matr); 
    
    if(sum(ind > 1)==0){
      break;
    }else{
      lambda[ind > 1] = lambda[ind > 1] + lambda.step;
      theta = result$hat.theta; #use warm start
      pi = result$hat.pi; #use warm start
    }
  }
  
  
  qindex <- apply(abs(result$hat.theta[,-1]), 1, which.max)+1;   #If still did not shrink to the sparsity level, we still give an answer based on the maximal coefficient. 
  ind <- (rowSums(abs(result$hat.theta[,-1])) < tol);
  qindex[ind] <- 1;
  Q = X[qindex,];
 
  list(lambda = lambda, Q = Q, num.iter = i, reg.est = result);
}






















