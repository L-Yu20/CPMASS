##To choose the alpha_n for balanced dataset
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
library(kcpRS)
remove(list = ls())

sourceDir <- function(path, trace = TRUE, ...){
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$"))
  {
    if (trace) cat(nm, ":")
    source(file.path(path, nm), ...)
  }
}
sourceDir("fn")
sourceCpp("F4DDS_binary.cpp")
sourceCpp("dF4DDS_binary.cpp")
sourceCpp("F4DDS_pairwise.cpp")
sourceCpp("dF4DDS_pairwise.cpp")

tic()

times = 100
n = 500
an0 = rep(0, 1)
an0[1] = floor(n / 5)
an0[2] = floor(n / 10)
an0[3] = floor(n / 15)
an0[4] = floor(n / 20)
an0[5] = floor(n / 25)
an0[6] = floor(n / 30)
cmz = matrix(nrow = length(an0), ncol = times)
lan = length(an0)

nuu1 = matrix(nrow = lan, ncol = 1000)
poo1 = array(dim = c(lan, 1000, 50))

##calculate the hatk, MSE and rand index
mk = rep(0, length(an0))  ##hatk
me = rep(0, length(an0))  ##MSE
ri = rep(0, length(an0))  ##rand index

for (kkk in 1:times){
  set.seed(1234 + kkk)
  ##balanced
  n = 500
  px = 50  ##px=50 or 100
  n1 = 125
  n2 = 250
  n3 = 375
  n4 = 500
  c1 = n1
  c2 = n2 - n1
  c3 = n3 - n2
  c4 = n4 - n3
  
  ##data setting
  af = 4
  px0 = floor(sqrt(px))  
  x = matrix(nrow = n, ncol = px)
  mu3 = rep(0, px - px0)
  sigma3 = diag(rep(1, px - px0))
  x[1:n, (px0 + 1):px] = mvrnorm(n = n, mu3, sigma3)
  ##change dimension
  mu1 = rep(0, px0)
  sigma1 = diag(rep(1, px0))
  rho = 7
  sigma2 = matrix(nrow = px0, ncol = px0)
  for (i in 1:px0){
    for (j in 1:px0){
      sigma2[i, j] = rho
    }
  }
  for (i in 1:px0){
    sigma2[i, i] = rho + 0.25
  }
  x[1:n1, 1:px0] = mvrnorm(n = c1, mu1, sigma2)
  x[(n1 + 1):n2, 1:px0] = mvrnorm(n = c2, mu1, sigma1)
  x[(n2 + 1):n3, 1:px0] = mvrnorm(n = c3, mu1, sigma2)
  x[(n3 + 1):n4, 1:px0] = mvrnorm(n = c4, mu1, sigma1)
  
  ##real label
  yreal = rep(1, n)
  yreal[1:n1] = 1
  yreal[(n1 + 1):n2] = 2
  yreal[(n2 + 1):n3] = 3
  yreal[(n3 + 1):n4] = 4
  
  for (ajk0 in 1:length(an0)){
    an = an0[ajk0]
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
      set.seed(1234 + kkk)
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
      if (k == 2){
        alpha10 = ybi(x, y, h1, B0)
      } else if (k > 2){
        alpha10 = ypa(x, y, h1, B0)
      }
      alphaj[1:px, (akk[ajk] + 1):(akk[ajk + 1])] = alpha10
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
    
    if (is.na(stop)) {
      stop <- length(dt)
    } 
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
    
    alpha100 = alphaj[1:px, (akk[hatd] + 1):(akk[hatd + 1])]
    BX1 = x %*% alpha100
    
    set.seed(1234 + kkk)
    
    out1 = e.divisive(BX1, R = 499, alpha = 1)
    kk1 = out1$estimates
    nuu1[ajk0, kkk] = length(kk1)
    poo1[ajk0, kkk, 1:nuu1[ajk0, kkk]] = kk1
  }
}

for (ajk0 in 1:length(an0)){
  kreal = 3
  numberr1 = nuu1[ajk0, ]
  positionn1 = poo1[ajk0, , ]
  for (k in 1:times){
    kk = positionn1[k, ]
    he1 = length(which(kk > 0))
    if (he1 == 2){
      positionn1[k, ] = 0
    } else {
      positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
      positionn1[k, (he1 - 1):he1] = c(0, 0)
    }
  }
  numberr1 = numberr1 - 2
  mk[ajk0] = mean(numberr1[1:times])
  mse1 = 0
  for (j in 1:times){
    a1 = (numberr1[j] - 3)^2
    mse1 = mse1 + a1
  }
  mse1 = mse1 / times
  mse1 = sqrt(mse1)
  me[ajk0] = mse1
  posi = positionn1[1:times, ]
  num = numberr1[1:times]
  ri[ajk0] = randi(yreal, posi, num, n, times)
}

result = cbind(mk, me, ri)
result
write.csv(result, "E:/re.csv")

toc()
