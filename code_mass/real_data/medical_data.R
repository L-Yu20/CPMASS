##Medical Data
setwd("E:/code_mass")
library(tictoc)
library(MASS)
library(class)
library(kernlab)
library(caret)
library(e1071)
library(pracma)
library(mcompanion)
library(ManifoldOptim)
library("R.matlab")
library("orthoDr")
library("Rcpp")
library(ecp)
library("wbs")
library("InspectChangepoint")
library("hdbinseg")
library('mvtnorm')
library(moments)
remove(list = ls())

sourceDir <- function(path, trace = TRUE, ...){
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$"))
  {
    if(trace) cat(nm, ":")
    source(file.path(path, nm), ...)
  }
}
sourceDir("fn")
sourceCpp("F4DDS_binary.cpp")
sourceCpp("dF4DDS_binary.cpp")
sourceCpp("F4DDS_pairwise.cpp")
sourceCpp("dF4DDS_pairwise.cpp")

tic()

set.seed(789)

data(ACGH, package = "ecp")
x = ACGH$data
dimx = dim(x)
n = dimx[1]
px = dimx[2]
tau = 0.5
an = floor(n / 10)

###
qm = 10
dt = seq(1, qm)
lambda = matrix(nrow = length(qm), ncol = 1)
lambda10 = matrix(nrow = length(qm), ncol = 1)
alphaj = matrix(nrow = px, ncol = sum(dt))
akk = rep(0, qm + 1)
for (i in 1:qm){
  akk[i + 1] = sum(seq(1, i))
}
qmax = min(px / 2, 10)
cn0 = 0.15 * log(n) * n^{-2 / (4 + qmax)}
stop <- NA
prev_lambda10 <- NA

for (ajk in 1:length(dt)){
  ##bandwidth
  d = dt[ajk]
  h1 = n^{-1 / (4 + d)}
  maxiter = 200
  
  ##set a random initial label
  zs = floor(n / an)
  k0 = rep(1, n)
  for (j in 1:zs){
    k0[((j - 1) * an + 1):(j * an)] = j
  }
  if ((n - zs * an) > 0){
    k0[(zs * an + 1):n] = zs + 1
  }
  y = k0
  
  ##the first dimension reduction
  p <- ncol(x)
  k = max(y)
  if (k == 2){
    Fhandle <- F(F4DDS_binary_cpp, x, y, h1)
  } else if (k > 2){
    Fhandle <- F(F4DDS_pairwise_cpp, x, y, h1)
  }
  beta <- as.matrix(orth(matrix(runif(p * d), p, d)))
  f <- Fhandle(beta)
  for (i in 1:maxiter){
    beta1 <- orth(matrix(runif(p * d), p, d))
    f1 <- Fhandle(beta1)
    if (f1 < f){
      f <- f1
      beta <- beta1
    }
  }
  B0 = beta
  ###find the dimension reduction matrix
  if (k == 2){
    alpha10 = ybi(x, y, h1, B0)
  } else if (k > 2){
    alpha10 = ypa(x, y, h1, B0)
  }
  alphaj[1:px, (akk[ajk] + 1):(akk[ajk + 1])] = alpha10
  ##first dimension reduction
  lambda[ajk] = Fhandle(alpha10)
  lambda10[ajk] = lambda[ajk] / (2 * n - 2 * an)
  
  if (!is.na(prev_lambda10)) {
    if (abs(prev_lambda10 - lambda10[ajk]) < cn0) {
      stop <- ajk
      break
    }
  }
  prev_lambda10 <- lambda10[ajk]
}

##
if (is.na(stop)) {
  stop <- length(dt)
} 
##
if (stop > 2){
  qm = stop
  lambda2 = lambda10[1:(qm - 1)] - lambda10[2:(qm)]
  lplus = lambda2 + cn0
  l3 = lplus[2:(qm - 1)] / lplus[1:(qm - 2)]
  l4 = which(l3 < 0.5)
  aaa = min(l3)
  if (aaa > 0.5){
    hatd = 1
  } else{
    hatd = max(l4) + 1
  }
} else {
  hatd = 1
}

##
alpha100 = alphaj[1:px, (akk[hatd] + 1):(akk[hatd + 1])]
BX1 = x %*% alpha100

set.seed(789)

##ecp
output1 = e.divisive(BX1, R = 499, alpha = 1)
kk1 = output1$estimates

##figure
kk2 = kk1
kk = length(kk1)
kk1 = kk1[2:(kk - 1)]
##
cpnum = kk - 2
plot(BX1[1:2215], main = "The E-Divisive method after MASS",
     xlab = "Index",
     ylab = expression(paste("B"["1"]^"T", "X"["i"])), cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')

toc()
