# Unbiasedness simulation using all possible samples under SRSWOR
library(utils)

N_A = 6
n_A = 3
n_A1 = 2

N_B = 6
n_B = 3
n_B1 = 2

N = N_A + N_B

A_set = combn(1:N_A, n_A)
B_set = combn(1:N_B, n_B)

x = rnorm(N, 2, 1)
y = x + rnorm(N)
t_y = sum(y)

res = NULL
for(k in 1:ncol(A_set)){
  print(k)
  for(k2 in 1:ncol(B_set)){
    A = A_set[,k]
    B = B_set[,k2]
    # print(A)
    
    # x_A = x[1:N_A]
    # x_B = x[(N_A + 1):N_B]
    
    y_sA = y[A]
    y_sB = y[N_A + B]
    
    x_sA = x[A]
    x_sB = x[N_A + B]
    
    y_s = c(y_sA, y_sB)
    x_s = c(x_sA, x_sB)
    model = lm(y_s ~ x_s, weights = c(rep(N_A / n_A, n_A), rep(N_B / n_B, n_B)))
    
    res_SRB = NULL
    A1_set = combn(A, n_A1)
    B1_set = combn(B, n_B1)
    for(l in 1:ncol(A1_set)){
      for(l2 in 1:ncol(B1_set)){
        A1 = A1_set[,l]
        A2 = A[!(A %in% A1)]
        
        B1 = B1_set[,l2]
        B2 = B[!(B %in% B1)]
        
        model1 = lm(y_s ~ 0 + x_s, data = data.frame(x_s = c(x[A1], x[N_A + B1]), y_s = c(y[A1], y[N_A + B1])), 
                    weights = c(rep(N_A / n_A, n_A1), rep(N_B / n_B, n_B1)))
        
        yhat = predict(model1, data.frame(x_s = x))
        
        res_SRB = c(res_SRB, sum(yhat) + sum(y[c(A1, N_A + B1)] - yhat[c(A1, N_A + B1)]) + 
                      sum(y[c(A2)] - yhat[c(A2)]) * (N_A - n_A1) / (n_A - n_A1) + sum(y[c(N_A + B2)] - yhat[c(N_A + B2)]) * (N_B - n_B1) / (n_B - n_B1))
      }
    }
    
    res = rbind(res, c(HT = sum(y_sA * N_A / n_A) + sum(y_sB * N_B / n_B),
                       GREG = sum(predict(model, data.frame(x_s = x))),
                       SRB = mean(res_SRB)))
  }
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
