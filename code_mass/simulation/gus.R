## Balanced dataset in Experiment 1
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

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
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

nu1 = matrix(nrow = 1, ncol = 1000)
po1 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nu2 = matrix(nrow = 1, ncol = 1000)
po2 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu1 = matrix(nrow = 1, ncol = 1000)
poo1 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu2 = matrix(nrow = 1, ncol = 1000)
poo2 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu3 = matrix(nrow = 1, ncol = 1000)
poo3 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu4 = matrix(nrow = 1, ncol = 1000)
poo4 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu5 = matrix(nrow = 1, ncol = 1000)
poo5 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu6 = matrix(nrow = 1, ncol = 1000)
poo6 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu7 = matrix(nrow = 1, ncol = 1000)
poo7 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu8 = matrix(nrow = 1, ncol = 1000)
poo8 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
nuu9 = matrix(nrow = 1, ncol = 1000)
poo9 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)

times = 100

for (kkk in 1:times) {
  set.seed(1234 + kkk)
  ##Balanced dataset
  n = 500
  px = 50 ##px=50 or 100
  n1 = 125
  n2 = 250
  n3 = 375
  n4 = 500
  c1 = n1
  c2 = n2 - n1
  c3 = n3 - n2
  c4 = n4 - n3
  
  an = floor(n / 10)
  ##data setting
  af = 4
  px0 = floor(sqrt(px))
  # px0 = 3 * floor(sqrt(px))
  x = matrix(nrow = n, ncol = px)
  mu3 = rep(0, px - px0)
  sigma3 = diag(rep(1, px - px0))
  x[(1):(n), (px0 + 1):px] = mvrnorm(n = n, mu3, sigma3)
  ##change dimension
  mu1 = rep(0, px0)
  # sigma1 = 0.5 * diag(rep(1, px0))
  sigma1 = diag(rep(1, px0))
  rho = 7
  # rho = 2
  sigma2 = matrix(nrow = px0, ncol = px0)
  for (i in 1:px0) {
    for (j in 1:px0) {
      sigma2[i, j] = rho
    }
  }
  for (i in 1:px0) {
    sigma2[i, i] = rho + 0.25
  }
  x[(1):(n1), 1:px0] = mvrnorm(n = c1, mu1, sigma2)
  x[(n1 + 1):(n2), 1:px0] = mvrnorm(n = c2, mu1, sigma1)
  x[(n2 + 1):(n3), 1:px0] = mvrnorm(n = c3, mu1, sigma2)
  x[(n3 + 1):(n4), 1:px0] = mvrnorm(n = c4, mu1, sigma1)
  
  ##real label
  yreal = rep(1, n)
  yreal[1:(n1)] = 1
  yreal[(n1 + 1):(n2)] = 2
  yreal[(n2 + 1):(n3)] = 3
  yreal[(n3 + 1):(n4)] = 4
  
  ###
  qm = 10
  dt = seq(1, qm)
  lambda = matrix(nrow = length(qm), ncol = 1)
  lambda10 = matrix(nrow = length(qm), ncol = 1)
  alphaj = matrix(nrow = px, ncol = sum(dt))
  akk = rep(0, qm + 1)
  for (i in 1:qm) {
    akk[i + 1] = sum(seq(1, i))
  }
  qmax = min(px / 2, 10)
  cn0 = 0.15 * log(n) * n^{-2 / (4 + qmax)}
  stop <- NA
  prev_lambda10 <- NA
  
  for (ajk in 1:length(dt)) {
    set.seed(1234 + kkk)
    ##bandwidth
    d = dt[ajk]
    h1 = n^{-1 / (4 + d)}
    maxiter = 200
    
    ##set a random initial label
    zs = floor(n / an)
    k0 = rep(1, n)
    for (j in 1:zs) {
      k0[((j - 1) * an + 1):(j * an)] = j
    }
    if ((n - zs * an) > 0) {
      k0[(zs * an + 1):n] = zs + 1
    }
    y = k0
    
    ##the first dimension reduction
    p <- ncol(x)
    k = max(y)
    if (k == 2) {
      Fhandle <- F(F4DDS_binary_cpp, x, y, h1)
    } else if (k > 2) {
      Fhandle <- F(F4DDS_pairwise_cpp, x, y, h1)
    }
    beta <- as.matrix(orth(matrix(runif(p * d), p, d)))
    f <- Fhandle(beta)
    for (i in 1:maxiter) {
      beta1 <- orth(matrix(runif(p * d), p, d))
      f1 <- Fhandle(beta1)
      if (f1 < f) {
        f <- f1
        beta <- beta1
      }
    }
    B0 = beta
    ###find the dimension reduction matrix
    if (k == 2) {
      alpha10 = ybi(x, y, h1, B0)
    } else if (k > 2) {
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
  if (stop > 2) {
    qm = stop
    lambda2 = lambda10[1:(qm - 1)] - lambda10[2:(qm)]
    lplus = lambda2 + cn0
    l3 = lplus[2:(qm - 1)] / lplus[1:(qm - 2)]
    l4 = which(l3 < 0.5)
    aaa = min(l3)
    if (aaa > 0.5) {
      hatd = 1
    } else {
      hatd = max(l4) + 1
    }
  } else {
    hatd = 1
  }
  
  ##
  alpha100 = alphaj[1:px, (akk[hatd] + 1):(akk[hatd + 1])]
  BX1 = x %*% alpha100
  
  ##Mahalanobis matrix
  Tx = matrix(nrow = n, ncol = n * px)
  a = matrix(nrow = n * n, ncol = px)
  Mx = matrix(nrow = px, ncol = px)
  for (i in 1:px) {
    tt = x[, i] %*% matrix(rep(1, n), nrow = 1, ncol = n) - matrix(rep(1, n), nrow = n, ncol = 1) %*% t(x[, i])
    Tx = matrix(tt, nrow = 1)
    a[, i] = Tx
  }
  Mx = t(a) %*% a / (n * (n - 1))
  sigma = matrix(rep(0, px * px), nrow = px, ncol = px)
  bn = floor(sqrt(n))
  zzs = floor(n / bn)
  for (j in 1:(zzs - 1)) {
    sigma = sigma + cov(x[((j - 1) * bn + 1):(j * bn), ])
  }
  sigma = sigma + cov(x[(((zs - 1) * bn) + 1):n, ])
  sigma = sigma / (zzs)
  Dx = Mx - 2 * sigma
  spec = eigen(Dx)
  evalue = rep(0, px + 1)
  va = abs(spec$values)
  evalue[1:px] = va
  evector = spec$vectors
  ccn = log(log(n)) * sqrt(px / n) / 2
  evalueplus = evalue + ccn
  lambda1 = evalueplus[2:(px + 1)] / evalueplus[1:px]
  lambda2 = which(lambda1[1:(px - 1)] < 0.5)
  aaa = min(lambda1[1:(px - 1)])
  if (aaa > 0.5) {
    hatq = 1
  } else {
    hatq = max(lambda2)
  }
  B = evector[, 1:(hatq)]
  BXM = t(B) %*% t(x)
  
  set.seed(1234 + kkk)
  
  ##ecpdr
  out1 = e.divisive(BX1, R = 499, alpha = 1)
  kk1 = out1$estimates
  nn1 = length(kk1)
  poo1[kkk, 1:nn1] = kk1
  nuu1[kkk] = nn1
  
  ##ecpm
  out2 = e.divisive(t(BXM), R = 499, alpha = 1)
  kk2 = out2$estimates
  nn2 = length(kk2)
  poo2[kkk, 1:nn2] = kk2
  nuu2[kkk] = nn2
  
  ##ecp
  out3 = e.divisive(x, R = 499, alpha = 1)
  kk3 = out3$estimates
  nn3 = length(kk3)
  poo3[kkk, 1:nn3] = kk3
  nuu3[kkk] = nn3
  
  ##kcpdr
  kk4 = kcpa(BX1, zs, 30)
  nuu4[kkk] = length(kk4)
  if (nuu4[kkk] > 0) {
    poo4[kkk, 1:(nuu4[kkk])] = kk4
  }
  
  ##kcpm
  kk5 = kcpa(t(BXM), zs, 200)
  nuu5[kkk] = length(kk5)
  if (nuu5[kkk] > 0) {
    poo5[kkk, 1:(nuu5[kkk])] = kk5
  }
  
  ##kcp
  kk6 = kcpa(x, zs, 5)
  nuu6[kkk] = length(kk6)
  if (nuu6[kkk] > 0) {
    poo6[kkk, 1:(nuu6[kkk])] = kk6
  }
  
  ##RSdr
  res = kcpRS(data = BX1, RS_fun = runVar, RS_name = "Var", wsize = 25,
              nperm = 1000, Kmax = 10, alpha = .05, varTest = FALSE, ncpu = 1)
  kk7 = res$changePoints
  nuu7[kkk] = length(kk7)
  if (nuu7[kkk] > 0) {
    poo7[kkk, 1:(nuu7[kkk])] = kk7
  }
  
  ##RSm
  res1 = kcpRS(data = t(BXM), RS_fun = runVar, RS_name = "Var", wsize = 25,
               nperm = 1000, Kmax = 10, alpha = .05, varTest = FALSE, ncpu = 1)
  kk8 = res1$changePoints
  nuu8[kkk] = length(kk8)
  if (nuu8[kkk] > 0) {
    poo8[kkk, 1:(nuu8[kkk])] = kk8
  }
  
  ##RS
  res2 = kcpRS(data = x, RS_fun = runVar, RS_name = "Var", wsize = 25,
               nperm = 1000, Kmax = 10, alpha = .05, varTest = FALSE, ncpu = 1)
  kk9 = res2$changePoints
  nuu9[kkk] = length(kk9)
  if (nuu9[kkk] > 0) {
    poo9[kkk, 1:(nuu9[kkk])] = kk9
  }
  
  ##
  plot(BX1[1:n])
  plot(BXM[1:n])
}

mk = rep(0, 9)
me = rep(0, 9)
ri = rep(0, 9)
kreal = 3

##ecpdr
numberr1 = nuu1
positionn1 = poo1
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):(he1)] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[1] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal) ^ 2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[1] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[1] = randi(yreal, posi, num, n, times)

##ecpm
numberr1 = nuu2
positionn1 = poo2
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):(he1)] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[2] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal) ^ 2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[2] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[2] = randi(yreal, posi, num, n, times)

##ecp
numberr1 = nuu3
positionn1 = poo3
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):(he1)] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[3] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal) ^ 2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[3] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[3] = randi(yreal, posi, num, n, times)

## kcpdr
numberr1 = nuu4
positionn1 = poo4
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):(he1)] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[4] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal)^2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[4] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[4] = randi(yreal, posi, num, n, times)

## kcpm
numberr1 = nuu5
positionn1 = poo5
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):(he1)] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[5] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal)^2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[5] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[5] = randi(yreal, posi, num, n, times)

## kcp
numberr1 = nuu6
positionn1 = poo6
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):(he1)] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[6] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal)^2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[6] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[6] = randi(yreal, posi, num, n, times)

## rsdr
numberr1 = nuu7
positionn1 = poo7
mk[7] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal)^2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[7] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[7] = randi(yreal, posi, num, n, times)

## rsm
numberr1 = nuu8
positionn1 = poo8
mk[8] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal)^2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[8] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[8] = randi(yreal, posi, num, n, times)

## rs
numberr1 = nuu9
positionn1 = poo9
mk[9] = mean(numberr1[1:times])
mse1 = 0
for (j in 1:times) {
  a1 = (numberr1[j] - kreal)^2
  mse1 = mse1 + a1
}
mse1 = mse1 / times
mse1 = sqrt(mse1)
me[9] = mse1
posi = positionn1[1:times, ]
num = numberr1[1:times]
ri[9] = randi(yreal, posi, num, n, times)

res = cbind(mk, me, ri)
res
write.csv(res, "E:/re.csv")

toc()
