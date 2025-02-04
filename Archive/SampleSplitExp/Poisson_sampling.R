# Unbiasedness simulation using all possible samples under Poisson
library(utils)

N = 10

Index = 1:N

l <- rep(list(c(F, T)), N)
A_set = apply(expand.grid(l), 1, function(x) Index[x])

pis = runif(N, 0.25, 0.75)
d = 1 / pis

x = rnorm(N, 2, 1)
y = x + rnorm(N) * 0.1
t_y = sum(y)

res = NULL
res2 = NULL
print(length(A_set))
for(k in 1:length(A_set)){
  A = A_set[[k]]
  A_c = Index[!(Index %in% A)]
  prob = prod(pis[A]) * prod(1 - pis[A_c])
  n = length(A)
  
  y_s = y[A]
  x_s = x[A]
  d_s = d[A]
  
  if(length(A) == 0){
    model = NULL
    y_GREG = 0
  }else{
    model = lm(y_s ~ 0 + x_s, weights = d_s)
    y_GREG = sum(predict(model, data.frame(x_s = x))) + 
      sum((y_s - predict(model, data.frame(x_s = x_s))) * d_s)
  }
  
  if(length(A) != 0){
    res_SRB = NULL
    if(length(A) == 1){
      A1_set = as.matrix(A, nr = 1, nc = 1)
    }
    else{
      # A1_set = combn(A, ceiling(length(A) / 2))
      A1_set = combn(A, 1)
    } 
    for(j in 1:ncol(A1_set)){
      A1 = A1_set[,j]
      A2 = A[!(A %in% A1)]
      
      if(F){
        print("new")
        print(A)
        print(A1_set)
        print(A1)
        print(A2)
      } 
      
      # prob2 = prod(pis[A1]) * prod(1 - pis[A2])

      d_s1 = d[A1]
      d_s2 = d[A2]
      
      model1 = lm(y_s ~ 0 + x_s, weights = d_s1, data = data.frame(x_s = x[A1], y_s = y[A1]))
      yhat = predict(model1, data.frame(x_s = x))
      res_SRB = c(res_SRB, (sum(yhat) + sum(y[A1] - yhat[A1]) + sum((y[A2] - yhat[A2]) * d_s2)) )
      
    }
  }else{
    res_SRB = 0
  }

  
  res = rbind(res, c(HT = sum(y_s * d_s),
                     GREG = y_GREG, SRB = mean(res_SRB)))
  res2 = c(res2, prob)
}

# apply(res, 2, var)
# 

colSums(res * res2)

colSums(res^2 * res2) - colSums(res * res2)^2

sum(res2)

t_y


combn(5,1)
# HT estimator
# GREG estimator
# SRB estimator
# Sample split estimator
# Sample split estimator2
