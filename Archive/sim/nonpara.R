N = 2000
n = 100
m = 50
p = 4
s = 3
g = N

set.seed(2)

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)

nonSRS = T

library(mvtnorm)
library(glmnet)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

# X = rmvnorm(n = N, rep(0, p), ar1_cor(p, 0.2))
# X = pnorm(X)

X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
# X = scale(X, T, T)

beta = c(rep(1, s), rep(0, p - s))
muX =  as.vector(X^2 %*% beta)
# muX =  as.vector(X %*% beta)
e = rnorm(N, 0, sd(muX) / 3)
y = muX + e
t_y = sum(y)
plot(y ~ muX)
# pi1 = e^2 / sum(e^2) * n

library(CVST)
library(mgcv)

b1 <- gam(y ~ te(X[,1], X[,2]))
b1$coefficients

plot(b1$fitted.values, y)
# krr <- constructKRRLearner()
# 
# dat <- constructData(X,y)
# X2 <- matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
# dat_tst <- constructData(X2 ,0)
# 
# par(mfrow=c(3,3),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)
# lambdas= 10^(-8:0)
# 
# for(lambda in lambdas){
#   # param <- list(kernel="rbfdot", sigma=0.0005, lambda=lambda)
#   param <- list(kernel="polydot", lambda=lambda, degree = 2, scale = 1, offset = 1)
#   # param <- list(kernel="vanilladot", lambda=lambda)
#   # param <- list(kernel="splinedot", lambda=lambda)
#   
#   krr.model <- krr$learn(dat,param)
#   pred <- krr$predict(krr.model,dat_tst)
#   plot(muX, y, xaxt='n', yaxt='n', main=paste('lambda =',signif(lambda,digits=3)) )
#   # lines(gridx,fgridx,col='red')
#   Index = order(as.vector(X2^2 %*% beta))
#   lines(as.vector(X2^2 %*% beta)[Index],pred[Index],col='blue')
# }

r = 1

if(nonSRS == TRUE){
  if(r == 0){
    z = rnorm(n = N, 0, 1)
  }else{
    ratio = (1 - r) / r
    z = rnorm(n = N, 0, sqrt(ratio)) - e
  }
  # z = -e
  # z = rnorm(n = N, 0, 1)
  z_sorted_idx = order(z)
  len = round(N / 4)
}

if(nonSRS == FALSE){
  pi1 = rep(n / N, N)
  # pi1 = rep(n / N, N) + rnorm(N, 0, 0.0025)
  # pi1 = rep(n / N, N) + e * n / N * 0.2
  
  # pi1 = (y - min(y)) / (max(y) - min(y)) / 2 + 0.25
  pi1 = (e - min(e)) / (max(e) - min(e)) / 2 + 0.10
  pi1 = pi1 / sum(pi1) * n
  d1 = 1 / pi1
  
  Omega = matrix(0, nr = n, nc = n)
  diag(Omega) = (1 - pi1) / pi1^2
}


res = NULL
res2 = NULL
res3 = NULL
SIMNUM = 500
X = cbind(1, X)
beta = c(0, beta)
for(simnum in 1:SIMNUM){
  print(simnum)
  set.seed(simnum)
  
  if(nonSRS == TRUE){
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
  }
  

  
  if(simnum == 1) plot(1 / d1, e, xlab = "Inclusion Probability", ylab = "error")
  
  if(nonSRS == FALSE) Index = sample(1:N, size = n, replace = FALSE, prob = pi1)
  y_s = y[Index]
  X_s = X[Index, ,drop = F]
  d1_s = d1[Index]
  lm_obj = lm(y_s ~  0 + X_s, weights = d1_s)
  beta_hat = lm_obj$coefficients
  #drop(solve(t(X_s) %*% diag(d1_s) %*% X_s, t(X_s) %*% diag(d1_s) %*% y_s))
  
  # beta_hat = unname(lm(y_s ~  t(apply(X_s, 1, function(k) k 
  # - colSums(X_s * d1_s) / N)), weights = d1_s)$coefficients[-1])
  
  y_HT = sum(y_s * d1_s)
  
  e_s = e[Index]
  y_diff = drop(sum(muX) + sum(e_s * d1_s))
  
  y_GREG = drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  
  # lm_obj = lm(y_s ~  0 + X_s, weights = d1_s - 1)
  # beta_hat = lm_obj$coefficients
  # y_Peff = drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  
  # library(glmnet)
  if(simnum == 1) cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  
  best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = best_lambda)
  beta_hat2 = as.vector(coef(best_model))[-1]
  
  y_Lasso = drop(colSums(X) %*% beta_hat2 + sum((y_s - drop(X_s %*% beta_hat2)) * d1_s))
  
  # if(simnum == 1){
  #   test <- gausspr(X_s, y_s, kernel = "rbfdot", scaled = F)
  #   kpar <- test@kernelf@kpar
  # }else{
  #   test <- gausspr(X_s, y_s, kernel = "rbfdot", scaled = F, kpar = kpar)
  # }
  
  if(simnum == 1){
    test <- gausspr(X_s, y_s, kernel = "polydot", scaled = F, kpar = list(degree = 2))
    kpar <- test@kernelf@kpar
  }else{
    test <- gausspr(X_s, y_s, kernel = "polydot", scaled = F, kpar = kpar)
  }
  
  # test <- gausspr(X_s, y_s, kernel = "splinedot", scaled = F)
  m_X <- predict(test,X)
  m_X_s <- predict(test,X_s)
  y_kernel = sum(m_X) + sum((y_s - m_X_s) * d1_s)
  
  # plot(m_X, y)
  # plot(m_X_s, y_s)
  
  # y_model = drop(colSums(X) %*% beta_hat)

  tmp <- sample(factor(rep(1:g, length.out=N), 
                labels=paste0(1:g)))
  split_group = by(1:N, tmp, c)
  
  # split_group = split(1:N, sample(1:g, N, replace=T))
  # split_group = split(1:N, 1:N)
  y_proposed = 0
  for(k in 1:length(split_group)){
    SubIndex = split_group[[k]]
    
    Index_sub1 = Index[Index %in% split_group[[k]]]
    Index_sub2 = Index[!(Index %in% Index_sub1)]
    
    # if(length(Index_sub1) != 0){
    y_s1 = y[Index_sub1]
    X_s1 = X[Index_sub1, ,drop = F]
    d1_s1 = d1[Index_sub1]
    
    y_s2 = y[Index_sub2]
    X_s2 = X[Index_sub2, ,drop = F]
    d1_s2 = d1[Index_sub2]
    lm_obj2 = lm(y_s2 ~ 0 +  X_s2, weights = d1_s2)
    beta_hat2 = lm_obj2$coefficients
    
    y_proposed = y_proposed + drop(colSums(X[SubIndex, , drop = F]) %*% beta_hat2  +
           sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1))
  }
  
  
  y_debiased_vec = NULL
  y_debiased_Lasso_vec = NULL
  y_debiased_Kernel_vec = NULL
  e_s_vec = NULL
  for(k in 1:10){
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
    
    y_debiased_Lasso_vec = c(y_debiased_Lasso_vec, y_debiased_Lasso)
    
    if(simnum == 1){
      test1 <- gausspr(X_s1, y_s1, kernel = "polydot", scaled = F, kpar = list(degree = 2))
      kpar1 <- test1@kernelf@kpar
    }else{
      test1 <- gausspr(X_s1, y_s1, kernel = "polydot", scaled = F, kpar = kpar1)
    }
    
    m_X2 <- predict(test1,X[-SubIndex, ])
    m_X_s2 <- predict(test1,X_s2)
    
    if(simnum == 1){
      test2 <- gausspr(X_s2, y_s2, kernel = "polydot", scaled = F, kpar = list(degree = 2))
      kpar2 <- test2@kernelf@kpar
    }else{
      test2 <- gausspr(X_s2, y_s2, kernel = "polydot", scaled = F, kpar = kpar2)
    }
    
    m_X1 <- predict(test2,X[SubIndex, ])
    m_X_s1 <- predict(test2,X_s1)
    
    y_debiased_Kernel = sum(m_X2) + sum((y_s2 - m_X_s2) * d1_s2) + sum(m_X1) + sum((y_s1 - m_X_s1) * d1_s1)
    y_debiased_Kernel_vec = c(y_debiased_Kernel_vec, y_debiased_Kernel)
  }
  
  e_s_tmp = colMeans(e_s_vec)
  
  y_debiased = mean(y_debiased_vec)
  
  y_debiased_Lasso = mean(y_debiased_Lasso_vec)
  
  y_debiased_Kernel = mean(y_debiased_Kernel_vec)
  
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
  
  sigma_Debiased = sqrt(drop(t(e_s_tmp)  %*% Omega %*% e_s_tmp))
  
  res = rbind(res, c(HT = y_HT, Diff = y_diff, GREG = y_GREG, Lasso = y_Lasso, Kernel = y_kernel, y_proposed = y_proposed,
                     Debiased = y_debiased, Lasso_Debiased = y_debiased_Lasso, Kernel_Debiased = y_debiased_Kernel))
  
  res2 = rbind(res2, c(HT = sigma_HT, Diff = sigma_diff, GREG = sigma_GREG, 
                       Debiased = sigma_Debiased))
  
}

BIAS = colMeans(res - t_y)
SE = apply(res, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE = apply(res - t_y, 2, function(x) sqrt(mean(x^2)))

# colMeans(res)

cbind(BIAS, SE, RMSE)
xtable::xtable(cbind(BIAS, SE, RMSE), caption = "Informative Sampmling")

# BIAS2 = colMeans(res2 - rep(SE[c(1:3, 5)], each = nrow(res2)))
# SE2 = apply(res2, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
# RMSE2 = apply(res2 - rep(SE[c(1:3, 5)], each = nrow(res2)), 2, function(x) sqrt(mean(x^2)))
# 
# cbind(BIAS2, REl_BIAS = BIAS2 / SE[c(1:3, 5)])

library(ggplot2)
ggplot() + 
  geom_density(data = data.frame(res), aes(x = Diff, fill = "Diff"), alpha = 0.3) +
  geom_density(data = data.frame(res), aes(x = GREG, fill = "GREG"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Lasso, fill = "Lasso"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Kernel, fill = "Kernel"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Debiased, fill = "GREG_Debiased"), alpha = 0.3) + 
  geom_density(data = data.frame(res), aes(x = Lasso_Debiased, fill = "Lasso_Debiased"), alpha = 0.3) +
  geom_density(data = data.frame(res), aes(x = Kernel_Debiased, fill = "Kernel_Debiased"), alpha = 0.3) +
  geom_vline(xintercept = t_y, col = "red")

