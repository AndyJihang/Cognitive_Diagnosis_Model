library(CDM);

K = 3;
J = 15;
N = 1000;

X = prof(K);
Z = interact(K); 

s.vec <- runif(J, 0.1, 0.3);
g.vec <- runif(J, 0.1, 0.3);

alpha.index <- sample(1:(2^K), N, replace = T); 
alpha <- X[alpha.index,];

q.index <- c(2,2,2,3,3,3,5,5,5,4,4,6,6,7,7); 
Q <- X[q.index,];

theta.true <- theta.fun(q.index, J, K, s.vec, g.vec);
pi.true <- rep(1/(2^K), 2^K)

data <- sim.din(q.matrix = Q, guess = g.vec, slip = s.vec, alpha = alpha)$dat;



save(data,theta.true,pi.true, file = "data.Rdata");

rm(list=ls(all=TRUE));

load("data.Rdata");
source("DINA_LASSO.R");
J = ncol(data);

lambda <- rep(0.01, J); #intial value of lambda
K = 3;

X = prof(K);
Z = interact(K); 

theta0 <- theta.true + runif(length(theta.true), -1, 1); 
pi0 <- rep(1/(2^K), 2^K); 

#Regularized EM for given lambda;
em.result <- EM(theta0, pi0, Z, data, lambda) 

#Running a path: Iteratively increase lambda, until each item has only one nonzero slope parameter. 
result <- Path(theta0, pi0, Z, data, lambda, X)




