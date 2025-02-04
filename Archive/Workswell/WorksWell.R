N = 2000
n = 100
m = 50
p = 40
s = 20

set.seed(2)

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)

library(mvtnorm)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

X = rmvnorm(n = N, rep(0, p), ar1_cor(p, 0.2))
X = pnorm(X)
# X = scale(X, T, T)

beta = c(rep(1, s), rep(0, p - s))
e = rnorm(N, 0, 1)
y = X %*% beta + e 
t_y = sum(y)

# pi1 = e^2 / sum(e^2) * n

r = 0.75

ratio = (1 - r) / r
z = rnorm(n = N, 0, sqrt(ratio)) + e
# z = rnorm(n = N, 0, 1)
z_sorted_idx = order(z)
len = round(N / 4)

# pi1 = rep(n / N, N)
# pi1 = rep(n / N, N) + rnorm(N, 0, 0.0025)
# pi1 = rep(n / N, N) + e * n / N * 0.2

# pi1 = (y - min(y)) / (max(y) - min(y)) / 2 + 0.25
# pi1 = (e - min(e)) / (max(e) - min(e)) / 2 + 0.10
# pi1 = pi1 / sum(pi1) * n
# d1 = 1 / pi1

res = NULL
res2 = NULL
res3 = NULL
SIMNUM = 500
X = cbind(1, X)
beta = c(0, beta)
for(simnum in 1:SIMNUM){
  print(simnum)
  # set.seed(simnum)
  
  Index = rep(0, n)
  d1 = rep(0, N)
  n_h = c(15, 20, 30, 35) / 100 * n
  cumn_h = cumsum(n_h)
  for (i in 1:4){
    Idx_z = (len * (i - 1)+ 1): (len * i)
    Idx = z_sorted_idx[Idx_z]
    d1[Idx] = len / n_h[i]
    from = ifelse(i == 1, 0, cumn_h[i-1])
    Index[(from + 1) : cumn_h[i]] = sample(Idx, size = n_h[i], replace = FALSE)
    
  }
  
  
  # Index = sample(1:N, size = n, replace = FALSE, prob = pi1)
  y_s = y[Index]
  X_s = X[Index, ,drop = F]
  d1_s = d1[Index]
  lm_obj = lm(y_s ~  0 + X_s, weights = d1_s)
  beta_hat = lm_obj$coefficients
  #drop(solve(t(X_s) %*% diag(d1_s) %*% X_s, t(X_s) %*% diag(d1_s) %*% y_s))
  
  # beta_hat = unname(lm(y_s ~  t(apply(X_s, 1, function(k) k 
  # - colSums(X_s * d1_s) / N)), weights = d1_s)$coefficients[-1])
  
  y_HT = sum(y_s * d1_s)
  
  y_diff = drop(colSums(X) %*% beta + sum((y_s - drop(X_s %*% beta)) * d1_s))
  
  y_GREG = drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  
  y_model = drop(colSums(X) %*% beta_hat)
  
  Index_sub1 = Index[sample(1:n, size = m, replace = FALSE)]
  Index_sub2 = Index[!(Index %in% Index_sub1)]
  
  y_s1 = y[Index_sub1]
  X_s1 = X[Index_sub1, ,drop = F]
  d1_s1 = d1[Index_sub1]
  lm_obj1 = lm(y_s1 ~  0 +X_s1, weights = d1_s1)
  beta_hat1 = lm_obj1$coefficients
  
  y_s2 = y[Index_sub2]
  X_s2 = X[Index_sub2, ,drop = F]
  d1_s2 = d1[Index_sub2]
  lm_obj2 = lm(y_s2 ~ 0 +  X_s2, weights = d1_s2)
  beta_hat2 = lm_obj2$coefficients
  
  # y_debiased = drop(colSums(X) %*% beta_hat + sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) 
  # + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
  
  y_debiased = drop(colSums(X) %*% (sum(d1_s2) / sum(d1_s) * beta_hat1 + sum(d1_s1) / sum(d1_s) * beta_hat2) + 
                      sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
  
  y_unbiased = drop(colSums(X) %*% beta_hat1 + sum((y_s1 - drop(X_s1 %*% beta_hat1))) 
                    + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * 2))
  
  Omega = matrix(1 - n / (n-1) * (N-1) / N, nr = n, nc = n)
  diag(Omega) = 1 - n / N
  
  
  sigma_HT = sqrt(drop(t(y_s)  %*% Omega %*% y_s * N^2 / n^2))
  
  e_s = y_s - X_s %*% beta
  sigma_diff = sqrt(drop(t(e_s)  %*% Omega %*% e_s * N^2 / n^2))
  
  # sum(sapply(1:n, function(k) sapply(1:n, function(l) 
  #   ifelse(k == l, 1 - n / N, 1 - n / (n - 1) * (N - 1) / N) * e_s[k] * e_s[l] * N / n * N / n  ))) / N^2
  
  e_s = y_s - X_s %*% beta_hat
  sigma_GREG = sqrt(drop(t(e_s)  %*% Omega %*% e_s * N^2 / n^2))
  
  # e_s = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
  # # e_s = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat1, y_s - X_s %*% beta_hat2)
  # # e_s = c(y_s1 - X_s1 %*% beta_hat2, y_s2 - X_s2 %*% beta_hat1)
  # sigma_Debiased = sqrt(drop(t(e_s)  %*% Omega %*% e_s * N^2 / n^2))
  
  e_s1 = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat1, 0)
  e_s2 = ifelse(Index %in% Index_sub2, y_s - X_s %*% beta_hat2, 0)
  
  sigma_Debiased = sqrt(drop((t(e_s1)  %*% Omega %*% e_s1 + t(e_s2)  %*% Omega %*% e_s2) * N^2 / n^2))
  
  that_k = sapply(1:n, function(k)  sum((y_s[-k]) * d1_s[-k]) * (n) / (n-1))
  sigma_HT2 = sqrt(sum((that_k- y_HT)^2) * (n - 1) / n)
  
  that_k = sapply(1:n, function(k)  sum((y_s[-k]) * d1_s[-k]) * (n) / (n-1)  + drop((colSums(X) - colSums(X_s[-k,,drop = F] * d1_s[-k]) * (n) / (n-1)) %*% beta))
  
  sigma_diff2 = sqrt(sum((that_k- y_diff)^2) * (n - 1) / n)
  
  that_k = sapply(1:n, function(k)  sum((y_s[-k]) * d1_s[-k]) * (n) / (n-1)  + drop((colSums(X) - colSums(X_s[-k,,drop = F] * d1_s[-k]) * (n) / (n-1)) %*% drop(solve(t(X_s[-k,,drop = F]) %*% diag(d1_s[-k]) %*% X_s[-k,,drop = F], t(X_s[-k,,drop = F]) %*% diag(d1_s[-k]) %*% y_s[-k]))))
  
  # that_k = sapply(1:n, function(k)  sum((y_s) * d1_s)  + drop((colSums(X) - colSums(X_s * d1_s)) %*% drop(solve(t(X_s[-k,]) %*% diag(d1_s[-k]) %*% X_s[-k,], t(X_s[-k,]) %*% diag(d1_s[-k]) %*% y_s[-k]))))
  
  sigma_GREG2 = sqrt(sum((that_k- y_GREG)^2) * (n - 1) / n)
  
  that_k1 = sapply(1:length(d1_s1), function(k) {betahat_k1 = drop(solve(t(X_s1[-k,,drop = F]) %*% diag(d1_s1[-k]) %*% X_s1[-k,,drop = F], t(X_s1[-k,,drop = F]) %*% diag(d1_s1[-k]) %*% y_s1[-k]))
  drop(colSums(X) %*% (sum(d1_s2) / sum(d1_s) * betahat_k1 + sum(d1_s1) / sum(d1_s) * beta_hat2) + 
         sum((y_s1[-k] - drop(X_s1[-k,,drop = F] %*% beta_hat2)) * d1_s1[-k]) + sum((y_s2 - drop(X_s2 %*% betahat_k1)) * d1_s2))
  })
  
  that_k2 = sapply(1:length(d1_s2), function(k) {betahat_k2 = drop(solve(t(X_s2[-k,,drop = F]) %*% diag(d1_s2[-k]) %*% X_s2[-k,,drop = F], t(X_s2[-k,,drop = F]) %*% diag(d1_s2[-k]) %*% y_s2[-k]))
  drop(colSums(X) %*% (sum(d1_s2) / sum(d1_s) * beta_hat1 + sum(d1_s1) / sum(d1_s) * betahat_k2) + 
         sum((y_s1 - drop(X_s1 %*% betahat_k2)) * d1_s1)) + sum((y_s2[-k] - drop(X_s2[-k,,drop = F] %*% beta_hat1)) * d1_s2[-k])
  })
  
  that_k_null = rep(0, length(that_k))
  that_k_null[Index %in% Index_sub1] = that_k1
  that_k_null[Index %in% Index_sub2] = that_k2
  
  sigma_Debiased2 = sqrt(sum((that_k_null - y_debiased)^2) * (n - 1) / n)
  
  res = rbind(res, c(Diff = y_diff, GREG = y_GREG, Model = y_model, Debiased = y_debiased, Unbiased = y_unbiased, HT = y_HT))
  
  res2 = rbind(res2, c(HT = sigma_HT, Diff = sigma_diff, GREG = sigma_GREG, Debiased = sigma_Debiased, HT2 = sigma_HT2, Diff2 = sigma_diff2, GREG2 = sigma_GREG2, Debiased2 = sigma_Debiased2))
  
  # res3 = rbind(res3, c(sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1), sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2)))
  
  res3 = rbind(res3, c(colMeans(X_s), beta_hat))
  
}

BIAS = colMeans(res - t_y)
SE = apply(res, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE = apply(res - t_y, 2, function(x) sqrt(mean(x^2)))

# colMeans(res2)

cbind(BIAS, SE, RMSE)

BIAS2 = colMeans(res2 - rep(SE[rep(c(6,1,2,4),2)], each = SIMNUM))
SE2 = apply(res2, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE2 = apply(res2 - rep(SE[rep(c(6,1,2,4),2)], each = SIMNUM), 2, function(x) sqrt(mean(x^2)))
MidErr2 = apply(res2 - rep(SE[rep(c(6,1,2,4),2)], each = SIMNUM), 2, median)

cbind(BIAS2, SE2, RMSE2, REl_BIAS = BIAS2 / SE[rep(c(6,1,2,4),2)])
cbind(BIAS2, REl_BIAS = BIAS2 / SE[rep(c(6,1,2,4),2)])

library(ggplot2)
ggplot() + 
  geom_density(data = data.frame(res), aes(x = Diff, fill = "Diff"), alpha = 0.3) +
  geom_density(data = data.frame(res), aes(x = GREG, fill = "GREG"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Debiased, fill = "Debiased"), alpha = 0.3) + 
  geom_vline(xintercept = t_y, col = "red")

