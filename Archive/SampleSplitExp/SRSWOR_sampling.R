# Unbiasedness simulation using all possible samples under SRSWOR
library(utils)

N = 20
n = 4
n1 = 2

Index = 1:N
A_set = combn(1:N, n)

x = rnorm(N, 2, 1)
y = x + rnorm(N) * 0.1
t_y = sum(y)

res = NULL
print(ncol(A_set))
for(k in 1:ncol(A_set)){
  A = A_set[,k]
  # print(A)
  
  y_s = y[A]
  x_s = x[A]
  model = lm(y_s ~ x_s)
  
  res_SRB = NULL
  A1_set = combn(A, n1)
  for(l in 1:ncol(A1_set)){
    A1 = A1_set[,l]
    A2 = A[!(A %in% A1)]
    
    model1 = lm(y_s ~ 0 + x_s, data = data.frame(x_s = x[A1], y_s = y[A1]))
    
    yhat = predict(model1, data.frame(x_s = x))
    
    res_SRB = c(res_SRB, sum(yhat) + sum(y[A1] - yhat[A1]) + sum(y[A2] - yhat[A2]) * (N - n1) / (n - n1))
  }
  
  res = rbind(res, c(HT = sum(y_s * N / n),
              GREG = sum(predict(model, data.frame(x_s = x))),
              SRB = mean(res_SRB)))
}

colMeans(res)
apply(res, 2, var)

apply(res, 2, mean)

t_y
# HT estimator
# GREG estimator
# SRB estimator
# Sample split estimator
# Sample split estimator2
