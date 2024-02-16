#####################################################################
# Tests for Applications 1 and 2
# Paper: Detecting p-hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied.
#####################################################################
#Install the following packages
# install.packages("NlcOptim")
# install.packages("fdrtool")
# install.packages("pracma")
# install.packages("gdata")
# install.packages("spatstat")
# install.packages("rddensity")

library("fdrtool")
library("pracma")
library("gdata")
library("spatstat")
library("rddensity")
library("ggplot2")
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

set.seed(123)

##### Binomial test #####
# use type = "c" to test on [p_min, p_max] 
# and type = "o" to test on (p_min, p_max)
Binomial <-function(P,p_min, p_max, type){
  # closed interval
  if (type == "c"){
    P = P[P<=p_max & P>=p_min]
  }
  # open interval
  if (type == "o"){
    P = P[P<p_max & P>p_min]
  }
  nn = length(P)
  kk = sum(P>(p_max+p_min)/2)
  return(1 - pbinom(kk-1, nn, 0.5))
}

##### LCM test #####

# Simulate Brownian Bridge (BB) and ||LCM(BB)-BB||
# M -- number of repetitions
SimBB <-function(M){
N = 10000
x = 1:N
x = x/N
BBsup = matrix(0, nrow=M, ncol = 1)
for (m in 1:M){
  eps <- rnorm(n = N, sd = 1, mean = 0)
  eps = eps/sqrt(N)
  W = cumsum(eps)
  B = W - x*W[N]
  C = c(0,x)
  B = c(0,B)
  lcmaj = gcmlcm(C, B, type="lcm")
  f_lcm = approxfun(lcmaj$x.knots, lcmaj$y.knots)
  y = f_lcm(C)
  BBsup[m,1] = max(abs(y - B))
}
return(BBsup)
}

# cdf for LCM test
F_LCMsup = ecdf(SimBB(10000))

#LCM test on [p_min, p_max]; P -- vector of p-values
LCM <-function(P, p_min, p_max){
  P = P[P<=p_max & P>=p_min]
  nn = length(P)
  f <- ecdf(P)
  x = seq(0, 1, length=1000)
  y = f(x*(p_max-p_min)+p_min)
  lcmaj = gcmlcm(x, y, type="lcm")
  ff_lcm = approxfun(lcmaj$x.knots, lcmaj$y.knots)
  z = as.numeric(ff_lcm(x))
  Test_stat = sqrt(nn)*max(abs(y-z)) #test statistics
  return(1 - F_LCMsup(Test_stat))    #p-value
}

##### Fisher's test #####
# P -- vector of p-values; p_min -- left end point of the testing interval,
#p_max -- right end point
Fisher <- function(P, p_min, p_max){
  P = P[P<p_max & P>=p_min]
  nn = length(P)
  statFM = -2*sum(log(1 - (P-p_min)/(p_max-p_min)))
  return(1 - pchisq(statFM, df = 2*nn))
}

##### Discontinuity test #####
#P -- vector of p-values; c -- potential discontinuity point
Discontinuity_test <- function(P, c){
  res = rddensity(P, c = c)
  return(res$test$p_jk)
}

##### Cox&Shi test #####

# Critical value function for 2-sided t-test
cv2 <- function(p){
  cv = qnorm(1 - p/2)
  return(cv)
}

lambda2 <- function(x1, x2, h){
  lambda = pnorm(cv2(x1) - h) - pnorm(cv2(x2) - h)+pnorm(cv2(x1) + h) - pnorm(cv2(x2) + h)
  return(lambda)
}

#Bounds on the proportions (two-sided t-test); 
#[p_min, p_max] interval with J bins
Bound0 <- function(p_min, p_max, J){
  h = seq(0,100,by = 0.001)
  X = linspace(p_min, p_max, J+1)
  B = matrix(0, J, 1)
  for (j in 1:(J)){
    Obj1 = lambda2(X[j], X[j+1], h)
    B[j] = max(Obj1)
  }
  if (p_min==0){
  B[1]=1
  }
  return(B)
}

#Bounds on the first differences of proportions (two-sided t-test); 
#[p_min, p_max] interval with J bins
Bound1 <- function(p_min, p_max, J){
  h = seq(0,100, by = 0.001)
  X = linspace(p_min, p_max, J+1)
  B = matrix(0, J-1, 1)
  for (j in 1:(J-1)){
    Obj1 = lambda2(X[j], X[j+1], h)
    Obj2 = lambda2(X[j+1], X[j+2], h)
    A = Obj2-Obj1
    B[j] = max(abs(A))
  }
  if (p_min==0){
  B[1]=1
  }
  return(B)
}

#Bounds on the second differences of proportions (two-sided t-test); 
#[p_min, p_max] interval with J bins
Bound2 <- function(p_min, p_max, J){
  h = seq(0,100, by = 0.001)
  X = linspace(p_min, p_max, J+1)
  B = matrix(0, J-2, 1)
  for (j in 1:(J-2)){
   Obj1 = lambda2(X[j], X[j+1], h)
   Obj2 = lambda2(X[j+1], X[j+2], h)
   Obj3 = lambda2(X[j+2], X[j+3], h)
   A = Obj3 - 2*Obj2 + Obj1
   B[j] = max(abs(A))
  }
  if (p_min==0){
  B[1]=1
  }
  return(B)
}
  
#######Bounds for One-sided t-tests###########
Bound0_1s <- function(pmax, J){
  incr = pmax/J
  X = linspace(0, pmax, J+1)
  B = matrix(0, J, 1)
  for (j in 1:(J)){
    B[j] = 2*pnorm((qnorm(1 - X[j]) - qnorm(1 - X[j+1]))/2)-1
  }
  B[1]=1
  return(B)
}
Bound1_1s <- function(pmax, J){
  incr = pmax/J
  h = linspace(0, 100, n=100000)
  X = linspace(0, pmax, J+1)
  B = matrix(0, J-1, 1)
  for (j in 1:(J-1)){
    Obj1 = pnorm(qnorm(1 - X[j]) - h) - pnorm(qnorm(1 - X[j+1]) - h)
    Obj2 = pnorm(qnorm(1 - X[j+1]) - h) - pnorm(qnorm(1 - X[j+2]) - h)
    A = Obj2-Obj1
    B[j] = max(abs(A))
  }
  B[1]=1
  return(B)
}
Bound2_1s <- function(pmax, J){
  incr = pmax/J
  h = linspace(0, 100, n=100000)
  X = linspace(0, pmax, J+1)
  B = matrix(0, J-2, 1)
  for (j in 1:(J-2)){
    Obj1 = pnorm(qnorm(1 - X[j]) - h) - pnorm(qnorm(1 - X[j+1]) - h)
    Obj2 = pnorm(qnorm(1 - X[j+1]) - h) - pnorm(qnorm(1 - X[j+2]) - h)
    A1 = Obj2-Obj1
    Obj3 = pnorm(qnorm(1 - X[j+1]) - h) - pnorm(qnorm(1 - X[j+2]) - h)
    Obj4 = pnorm(qnorm(1 - X[j+2]) - h) - pnorm(qnorm(1 - X[j+3]) - h)
    A2 = Obj4-Obj3
    B[j] = max(abs(A2-A1))
  }
  B[1]=1
  return(B)
}

#Cox-Shi test for K-monotonicity and bounds on [p_min, p_max] interval
#Q - vector of p-values; ind - vector of paper ids; J - number of subintervals;
# B=1 to use bounds, B=0 to test without bounds

CoxShi <- function(Q, ind, p_min, p_max, J, K, B){
  B0 = Bound0(p_min, p_max, J)
  B1 = Bound1(p_min, p_max, J)
  B2 = Bound2(p_min, p_max, J)

  P = Q[Q<=p_max & Q>=p_min]
  if (length(ind)>1){
    ind = ind[Q<=p_max & Q>=p_min]
    indu = unique(ind)
  }
  Bnd_adj = length(P)/length(Q)
  N = length(P)
  bin = seq(p_min,p_max,length=J+1)
  Phat = matrix(0, nrow=J-1, ncol = 1)
  
  for (s in 1:(J-1)){
    Phat[s] = sum((P>bin[s])*(P<=bin[s+1]))/N
  }
  Phat[1] = Phat[1]+sum(P==bin[1])/N
  if (B==0){
    B0 = -matrix(1, nrow = J, ncol = 1)
  }
  if (B==1){
    B0=-B0/Bnd_adj
    B1=-B1/Bnd_adj
    B2=-B2/Bnd_adj
    if (p_min==0){
      B0[1]=-1
      B1[1]=-1
      B2[1]=-1
    }
  }
  
  if (length(ind)>1){
    Omega = matrix(0, J-1, J-1)
    for (i in c(indu)){
      X = P[ind==i]
      
      mq = repmat(matrix(X, 1, length(X)), J-1,1)
      a1 = (mq<= bin[2:J])
      a2 = (mq>bin[1:(J-1)])
      mq0 = 1*(mq==0)
      mq0[2:(J-1),] = 0 
      mq = 1*(a1*a2) + mq0
      mq = mq - repmat((Phat), 1,length(X))
      
      Omega = Omega + (mq)%*%matrix(1, length(X), length(X))%*%t(mq)
    }
    Omega = Omega/length(P)
  }
  if (length(ind)==1){
    if (min(Phat)==0){
      Qhat = Phat*N/(N+1) + 1/(J*(N+1))
      Omega = diag(c(Qhat)) - Qhat%*%t(Qhat)
    }
    if (min(Phat)>0){
      Omega = diag(c(Phat)) - Phat%*%t(Phat)
    }
  }
  D = matrix(0, J-1, J)
  for (i in 1:(J-1)){
    for (j in 1:J){
      if (i==j){
        D[i,j] = -1
      }
      if (i+1==j){
        D[i,j] = 1
      }
      
    }
  }
  Dk = -D
  if (K>1){
    d = D
    for (k in 2:K){
      d = D[1:(J-k), 1:(J-k+1)]%*%d
      Dk = rbind(Dk, (-1)^k*d)
    }
  }
  if (B==0){
    Dk = rbind(-diag(J), diag(J), Dk)
  }
  if (B==1){
    Dk = rbind(-diag(J),-Dk, diag(J), Dk)
  }
  eJ = matrix(0, J, 1)
  eJ[J] = 1
  
  F1 = rbind(-diag(J-1), matrix(1, 1,J-1))
  c = matrix(0, (K+1)*(J-K/2), 1)

  if (B==0){
    c = rbind(matrix(-1, J,1), c)
  }
  if (B==1){
    
    if (K==0){
      c=rbind(B0, c)
    }
    if (K==1){
      c=rbind(B0, B1, c)
    }
    if (K==2){
      c=rbind(B0, B1, B2, c)
    }
  }
  A = Dk%*%F1
  b = Dk%*%eJ - c
  
  if (abs(det(Omega))>0){
  myQ = solve(Omega)
  fn <- function(t){
    Obj = N*(t(Phat - t))%*%myQ%*%((Phat - t))
    return(Obj)
  }
  t0 = matrix(1, (J-1), 1)/(J)
  res = fmincon(t0, fn, A = A, b = b)
  t_opt = t(t(res$par))
  print(t_opt)
  T = fn(t_opt)
  Ba = A[which(res$info$lambda$ineqlin>0),]
  JX = qr(Ba)$rank
  if (res$convergence==0){
    return(1 - pchisq(T, df = JX)*(JX>0))
  }
  else {
    return(999)
  }
  }
  else {
    return (888)
  }
}
