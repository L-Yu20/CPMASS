randi <- function(yreal, posi, num, n, times) {
  numberr1 = num
  positionn1 = posi
  ## rand index
  index = matrix(nrow = 1, ncol = n)
  index0 = yreal
  randindex = matrix(nrow = 1, ncol = times)
  all0 = n * (n - 1) / 2
  for (j in 1:times) {
    if (numberr1[j] == 1) {
      index[1:positionn1[j, 1]] = 1
      index[(positionn1[j, (numberr1[j])] + 1):n] = ((numberr1[j]) + 1) 
    } else if (numberr1[j] == 0) {
      index[1:n] = 1
    } else {
      index[1:positionn1[j, 1]] = 1
      index[(positionn1[j, (numberr1[j])] + 1):n] = ((numberr1[j]) + 1)
      for (i in 1:((numberr1[j]) - 1)) {
        index[(positionn1[j, i] + 1):positionn1[j, i + 1]] = i + 1
      }
    }
    ## rand index = (tp + tn) / (n * (n - 1) / 2)
    tp = 0
    tn = 0
    for (jj in 2:n) {
      for (hh in 1:(jj - 1)) {
        if (index[jj] == index[hh] & index0[jj] == index0[hh])
          tp = tp + 1 
        else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
      } 
    }
    randindex[j] = (tp + tn) / all0
  }
  ## the mean of the rand index
  randi = mean(randindex)
  
  return(randi)
}

