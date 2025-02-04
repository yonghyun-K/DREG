N = 2000
n = 140
m = 70
p = 40
s = 5

set.seed(2)

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
library(xtable)
library(mvtnorm)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

X = rmvnorm(n = N, rep(0, p), ar1_cor(p, 0.2)) + 2
# X = pnorm(X)

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)

# X = X * c(rep(1, N), rep(-1, N))

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
# X = scale(X, T, T)

beta = c(rep(1, s), rep(0, p - s))
e = rnorm(N, 0, 1)
y = X %*% beta + e 
t_y = sum(y)

# pi1 = e^2 / sum(e^2) * n

r = 0
if(r == 0){
  z = rnorm(n = N, 0, 1)
}else{
  z = rnorm(n = N, 0, sqrt((1 - r) / r)) + e
}
# z = rnorm(n = N, 0, sqrt(ratio)) + e
# z = e
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
  set.seed(simnum)
  
  Index = rep(0, n)
  d1 = rep(0, N)
  n_h = c(15, 20, 30, 35) / 100 * n
  # n_h = c(25, 25, 25, 25) / 100 * n
  cumn_h = cumsum(n_h)
  
  Omega = matrix(0, nr = n, nc = n)
  
  for (i in 1:4){
    Idx_z = (len * (i - 1)+ 1): (len * i)
    Idx = z_sorted_idx[Idx_z]
    d1[Idx] = len / n_h[i]
    from = ifelse(i == 1, 0, cumn_h[i-1])
    Index[(from + 1) : cumn_h[i]] = sample(Idx, size = n_h[i], replace = FALSE)
    # len^2 / n_h[i]^2 - len * (len - 1) / n_h[i] / (n_h[i] - 1) 
    # (1 - n_h[i] / (n_h[i]-1) * (len-1) / len)* len^2 / n_h[i]^2
    Omega_tmp = matrix(len^2 / n_h[i]^2 - len * (len - 1) / n_h[i] / (n_h[i] - 1) , nr = n_h[i], nc = n_h[i])
    diag(Omega_tmp) = (1 - n_h[i] / len)  * len^2 / n_h[i]^2
    
    Omega[(from + 1) : cumn_h[i], (from + 1) : cumn_h[i]] = Omega_tmp
    
  }
  
  if(simnum == 1) plot(1 / d1, e, xlab = "Inclusion Probability", ylab = "error")
  
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
  
  library(glmnet)
  if(simnum == 1) cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  
  best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = best_lambda)
  beta_hat_Lasso = as.vector(coef(best_model))[-1]
  
  y_Lasso = drop(colSums(X) %*% beta_hat_Lasso + sum((y_s - drop(X_s %*% beta_hat_Lasso)) * d1_s))
  
  # y_model = drop(colSums(X) %*% beta_hat)
  
  y_debiased_vec = NULL
  y_debiased_Lasso_vec = NULL
  e_s_vec = NULL
  e_s_vec_Lasso = NULL
  for(k in 1:10){
    # for(k in 1){
    set.seed(k)
    SubIndex = sample(1:N, size = round(N / 2), replace = FALSE)
    
    Index_sub1 = Index[Index %in% SubIndex]
    Index_sub2 = Index[!(Index %in% SubIndex)]
    # Index_sub1 = Index[sample(1:n, size = m, replace = FALSE)]
    # Index_sub2 = Index[!(Index %in% Index_sub1)]
    
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
    
    # y_debiased = drop(colSums(X) %*% (sum(d1_s2) / sum(d1_s) * beta_hat1 + sum(d1_s1) / sum(d1_s) * beta_hat2) + 
    #                     sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
    
    t1_hat = colSums(X[SubIndex, ]) %*% beta_hat2 + sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1)
    
    y_debiased = drop(colSums(X[SubIndex, ]) %*% beta_hat2 + colSums(X[-SubIndex, ]) %*% beta_hat1 +
                        sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))

    # y_debiased = drop(colSums(X[SubIndex, ]) %*% beta_hat2 +
    #                     sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1))*2
    
        
    e_s_tmp = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
    e_s_vec = rbind(e_s_vec, e_s_tmp)
    
    if(is.na(y_debiased)) stop()
    
    y_debiased_vec = c(y_debiased_vec, y_debiased)

    if(simnum == 1) cv_model <- cv.glmnet(X_s1, y_s1, weights = d1_s1)

    #find optimal lambda value that minimizes test MSE
    best_lambda1 <- cv_model$lambda.min
    best_lambda1

    best_model1 <- glmnet(X_s1, y_s1, weights = d1_s1, lambda = best_lambda1)
    beta_hat1 = as.vector(coef(best_model1))[-1]

    if(simnum == 1) cv_model <- cv.glmnet(X_s2, y_s2, weights = d1_s2)

    #find optimal lambda value that minimizes test MSE
    best_lambda2 <- cv_model$lambda.min
    best_lambda2

    best_model2 <- glmnet(X_s2, y_s2, weights = d1_s2, lambda = best_lambda2)
    beta_hat2 = as.vector(coef(best_model2))[-1]

    y_debiased_Lasso = drop(colSums(X[SubIndex, ]) %*% beta_hat2 + colSums(X[-SubIndex, ]) %*% beta_hat1 +
                              sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
    
    e_s_tmp_Lasso = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
    e_s_vec_Lasso = rbind(e_s_vec_Lasso, e_s_tmp_Lasso)
    
    y_debiased_Lasso_vec = c(y_debiased_Lasso_vec, y_debiased_Lasso)
  }
  
  e_s_tmp = colMeans(e_s_vec)
  e_s_tmp_Lasso = colMeans(e_s_vec_Lasso)
  
  y_debiased = mean(y_debiased_vec)
  
  y_debiased_Lasso = mean(y_debiased_Lasso_vec)
  
  # y_unbiased = drop(colSums(X) %*% beta_hat1 + sum((y_s1 - drop(X_s1 %*% beta_hat1))) 
  #                   + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * 2))
  
  
  sigma_HT = sqrt(drop(t(y_s)  %*% Omega %*% y_s))
  
  # tmp <- unclass(by(y_s, d1_s, sd))
  # attr(tmp, "call") <- NULL
  # tmp <- tmp[4:1]
  # sqrt(sum(len^2 * (1 - n_h / len) / n_h * tmp^2))
  
  e_s = y_s - X_s %*% beta
  sigma_diff = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
  
  # sum(sapply(1:n, function(k) sapply(1:n, function(l) 
  #   ifelse(k == l, 1 - n / N, 1 - n / (n - 1) * (N - 1) / N) * e_s[k] * e_s[l] * N / n * N / n  ))) / N^2
  
  e_s = y_s - X_s %*% beta_hat
  sigma_GREG = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
  
  e_s = y_s - X_s %*% beta_hat_Lasso
  sigma_Lasso = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
  
  sigma_Debiased = sqrt(drop(t(e_s_tmp)  %*% Omega %*% e_s_tmp))
  # e_s1 = ifelse(Index %in% Index_sub1, y_s - X_s %*% lm_obj2$coefficients, 0)
  # e_s2 = ifelse(Index %in% Index_sub2, y_s - X_s %*% lm_obj1$coefficients, 0)
  # sigma_Debiased = sqrt(drop((t(e_s1)  %*% Omega %*% e_s1 + t(e_s2)  %*% Omega %*% e_s2)))
  
  sigma_Debiased_Lasso = sqrt(drop(t(e_s_tmp_Lasso)  %*% Omega %*% e_s_tmp_Lasso))
  
  y_res = c(HT = y_HT, Diff = y_diff, GREG = y_GREG, Lasso = y_Lasso, 
            Double = y_debiased, DoubleL = y_debiased_Lasso)
  
  res = rbind(res, y_res)
  
  sigma_res = c(HT = sigma_HT, Diff = sigma_diff, GREG = sigma_GREG, Lasso = sigma_Lasso,
                Double = sigma_Debiased, DoubleL = sigma_Debiased_Lasso)
  
  res2 = rbind(res2, sigma_res)
  
  res3 = rbind(res3, ifelse(abs(y_res - t_y) > 1.96 * sigma_res, 0, 1))
}

BIAS = colMeans(res - t_y)
SE = apply(res, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE = apply(res - t_y, 2, function(x) sqrt(mean(x^2)))

# colMeans(res)
tmpdf21 = cbind(BIAS, SE, RMSE)
xtable(cbind(BIAS, SE, RMSE), digits = 3, caption = "Summary of point estimation")

BIAS2 = colMeans(res2 - rep(SE, each = nrow(res2)))
SE2 = apply(res2, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE2 = apply(res2 - rep(SE, each = nrow(res2)), 2, function(x) sqrt(mean(x^2)))

# cbind(BIAS2, REl_BIAS = BIAS2 / SE)
tmpdf22 = cbind(RB = BIAS2 / SE, CR = colMeans(res3))
xtable(cbind(RB = BIAS2 / SE, CR = colMeans(res3)), digits = c(0,4,3), caption = "Summary of variance estimation")

library(ggplot2)
ggplot() + 
  geom_density(data = data.frame(res), aes(x = Diff, fill = "Diff"), alpha = 0.3) +
  geom_density(data = data.frame(res), aes(x = GREG, fill = "GREG"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Lasso, fill = "Lasso"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Debiased, fill = "GREG_Debiased"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Lasso_Debiased, fill = "Lasso_Debiased"), alpha = 0.3) +
  geom_vline(xintercept = t_y, col = "red")


mat = matrix((n - N) / n / (N-1), nr = N, nc = N)
diag(mat) = 0

mat2 = matrix((n - N) / n / (N-1), nr = N, nc = N)
diag(mat2) = (N - n) / n

mat = matrix(0, nr = N, nc = N)
diag(mat) = 0

mat2 = matrix(0, nr = N, nc = N)
diag(mat2) = d1 - 1

# N^(-2) * sum((d1 - 1) * X[,2] * X[,3] * e)
# cov(X[,2], X[,3] * e) / N

e_res = lm(y ~ X)$residuals

tmp1= diag(t(X*e_res) %*% mat %*% X) / N/ (N-1) - diag(cov(X, X*e_res)) / N
tmp2 = diag(t(X*e_res) %*% mat2 %*% X) / N / (N-1)

tmp1= t(X*e_res) %*% mat %*% X / N/ (N-1) - cov(X, X*e_res) / N
tmp2 = t(X*e_res) %*% mat2 %*% X / N / (N-1)
diag(t(X) %*% mat2 %*% (X*e_res) / N / (N-1))
plot(tmp1, tmp2)

sum(diag((t(X) %*% X) %*% t(X*e_res) %*% diag(d1 - 1) %*% X)) / N^2

# plot(diag(cov(X, X*e_res)) / N, diag(t(X*e_res) %*% diag(d1 - 1) %*% X)/ N^2)
abline(a = 0, b= -1)
colSums(X) - colSums(X, weights = d)
