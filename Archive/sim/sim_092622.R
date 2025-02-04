N = 2000
n = 200
BIAS_res = NULL
SE_res = NULL
pvec = 1:18 * 10
set.seed(10)
# ar1_cor <- function(n, rho) {
#   exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
#                     (1:n - 1))
#   rho^exponent
# }
# 
# X_whole = matrix(mvtnorm::rmvnorm(N, sigma = ar1_cor(200, 0.8)), nr = N, nc = 200)

X_whole = matrix(rnorm(N * 500, 200, 1), nr = N, nc = 500)

for(p in pvec){
  # p = 10
  set.seed(1)
  print(p)
s = 5

X = X_whole[,1:p]
beta = c(rep(1, s), rep(0, p - s))
e = rnorm(N, 0, 1)
y = X %*% beta + e
t_y = sum(y)

# pi1 = X[,1]^2 / sum(X[,1]^2) * n
pi1 = rep(n / N, N)
d1 = 1 / pi1

res = NULL
res2 = NULL
res3 = NULL
SIMNUM = 500

library(glmnet)

for(simnum in 1:SIMNUM){
  # set.seed(simnum)
  Index = sample(1:N, size = n, replace = FALSE, prob = pi1)
  y_s = y[Index]
  X_s = X[Index, ]
  d1_s = d1[Index]
  lm_obj = lm(y_s ~  X_s, weights = d1_s)
  beta_hat = lm_obj$coefficients[-1]
  #drop(solve(t(X_s) %*% diag(d1_s) %*% X_s, t(X_s) %*% diag(d1_s) %*% y_s))
  
  # beta_hat = unname(lm(y_s ~  t(apply(X_s, 1, function(k) k 
  # - colSums(X_s * d1_s) / N)), weights = d1_s)$coefficients[-1])
  
  y_HT = sum(y_s * d1_s)
  
  y_diff = drop(colSums(X) %*% beta + sum((y_s - drop(X_s %*% beta)) * d1_s))
  
  y_GREG = drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  
  y_model = drop(colSums(X) %*% beta_hat)

  
  #perform k-fold cross-validation to find optimal lambda value
  if(simnum == 1) cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  
  best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = best_lambda)
  beta_hat2 = as.vector(coef(best_model))[-1]
  
  y_Lasso = drop(colSums(X) %*% beta_hat2 + sum((y_s - drop(X_s %*% beta_hat2)) * d1_s))
  
  # sigma_Debiased = sqrt(drop(t(e_s)  %*% Omega %*% e_s * N^2 / n^2))
  
  res = rbind(res, c(HT = y_HT, Diff = y_diff, GREG = y_GREG, Lasso = y_Lasso))

  
}

BIAS = colMeans(res - t_y)
SE = apply(res, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE = apply(res - t_y, 2, function(x) sqrt(mean(x^2)))

cbind(BIAS, SE, RMSE)

BIAS_res = rbind(BIAS_res, BIAS)
SE_res = rbind(SE_res, SE)
}

# plot(pvec, SE_res[,3], type = "l", col = 1, ylim = c(min(SE_res), max(SE_res)))
# points(pvec, SE_res[,2], type = "l", col = 2)
# points(pvec, SE_res[,4], type = "l", col = 3)
# points(pvec, SE_res[,1], type = "l", col = 4)
# legend("topleft", c("GREG", "Diff", "HT"), fill = c(1,2,3, 4))
# 
# plot(pvec, BIAS_res[,3], type = "l", col = 1, ylim = c(min(BIAS_res), max(BIAS_res)))
# points(pvec, BIAS_res[,2], type = "l", col = 2)
# points(pvec, BIAS_res[,1], type = "l", col = 4)


plot(pvec, SE_res[,3], type = "l", col = 1, ylim = c(min(BIAS_res), max(SE_res)), main = sprintf("n = %g", n), xlab = "p", ylab = " ")
points(pvec, SE_res[,2], type = "l", col = 2)
points(pvec, SE_res[,4], type = "l", col = 3)
points(pvec, SE_res[,1], type = "l", col = 4)
legend("topleft", c("GREG", "Diff", "Lasso", "HT"), col = c(1,2,3,4), lty = 1, title = "SE")

points(pvec, BIAS_res[,3], type = "l", col = 1, lty = 2)
points(pvec, BIAS_res[,2], type = "l", col = 2, lty = 2)
points(pvec, BIAS_res[,4], type = "l", col = 3, lty = 2)
points(pvec, BIAS_res[,1], type = "l", col = 4, lty = 2)
legend("bottomleft", c("GREG", "Diff", "Lasso", "HT"), col = c(1,2,3, 4), lty = 2, title = "BIAS")

