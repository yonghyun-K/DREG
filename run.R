Index = rep(0, n)
d1 = rep(0, N)
n_h = c(15, 20, 30, 35) / 100 * n
# n_h = c(25, 25, 25, 25) / 100 * n
cumn_h = cumsum(n_h)

Omega = matrix(0, nr = n, nc = n)

for (i in 1:4) {
  Idx_z = (len * (i - 1) + 1):(len * i)
  Idx = z_sorted_idx[Idx_z]
  d1[Idx] = len / n_h[i]
  from = ifelse(i == 1, 0, cumn_h[i - 1])
  Index[(from + 1):cumn_h[i]] = sample(Idx, size = n_h[i], replace = FALSE)
  # len^2 / n_h[i]^2 - len * (len - 1) / n_h[i] / (n_h[i] - 1)
  # (1 - n_h[i] / (n_h[i]-1) * (len-1) / len)* len^2 / n_h[i]^2
  Omega_tmp = matrix(len ^ 2 / n_h[i] ^ 2 - len * (len - 1) / n_h[i] / (n_h[i] - 1) ,
                     nr = n_h[i],
                     nc = n_h[i])
  diag(Omega_tmp) = (1 - n_h[i] / len)  * len ^ 2 / n_h[i] ^ 2
  
  Omega[(from + 1):cumn_h[i], (from + 1):cumn_h[i]] = Omega_tmp
  
}

# if(simnum == 1) plot(1 / d1, e, xlab = "Inclusion Probability", ylab = "error")

# Index = sample(1:N, size = n, replace = FALSE, prob = pi1)
y_s = y[Index]
X_s = X[Index, , drop = F]
d1_s = d1[Index]
if (wls) {
  lm_obj = lm(y_s ~  0 + X_s, weights = d1_s)
} else{
  lm_obj = lm(y_s ~  0 + X_s)
}
beta_hat = lm_obj$coefficients
#drop(solve(t(X_s) %*% diag(d1_s) %*% X_s, t(X_s) %*% diag(d1_s) %*% y_s))

X_s_oracle = X[Index, 1:(1 + s), drop = F]
if (wls) {
  lm_obj_oracle = lm(y_s ~  0 + X_s_oracle, weights = d1_s)
} else{
  lm_obj_oracle = lm(y_s ~  0 + X_s_oracle)
}
beta_hat_oracle = lm_obj_oracle$coefficients

y_HT = sum(y_s * d1_s)

y_diff = drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))

y_GREG = drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))

y_GREG_oracle = drop(colSums(X[, 1:(1 + s), drop = F]) %*% beta_hat_oracle + sum((y_s - drop(X_s_oracle %*% beta_hat_oracle)) * d1_s))

if(simnum == 0){
  if(wls){
    cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
  }else{
    cv_model <- cv.glmnet(X_s, y_s)
  }
}


#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda
if (wls) {
  best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = best_lambda)
} else{
  best_model <- glmnet(X_s, y_s, lambda = best_lambda)
}
beta_hat_Lasso = as.vector(coef(best_model))[-1]
y_Lasso = drop(colSums(X) %*% beta_hat_Lasso + sum((y_s - drop(
  X_s %*% beta_hat_Lasso
)) * d1_s))

# sum(beta_hat_Lasso[-c(1 : (1 + s))] != 0) / (p - s) # FDR

X_s_refit = X_s[, beta_hat_Lasso != 0, drop = F]
if (wls) {
  refit_model = lm(y_s ~  0 + X_s_refit, weights = d1_s)
} else{
  refit_model = lm(y_s ~  0 + X_s_refit)
}
beta_hat_refit = refit_model$coefficients
y_Lasso2 = drop(colSums(X[, beta_hat_Lasso != 0, drop = F]) %*% beta_hat_refit + sum((y_s - drop(
  X_s_refit %*% beta_hat_refit
)) * d1_s))

sigma_HT = sqrt(drop(t(y_s)  %*% Omega %*% y_s))

# tmp <- unclass(by(y_s, d1_s, sd))
# attr(tmp, "call") <- NULL
# tmp <- tmp[4:1]
# sqrt(sum(len^2 * (1 - n_h / len) / n_h * tmp^2))

e_s = y_s - mu[Index]
sigma_diff = sqrt(drop(t(e_s)  %*% Omega %*% e_s))

# sum(sapply(1:n, function(k) sapply(1:n, function(l)
#   ifelse(k == l, 1 - n / N, 1 - n / (n - 1) * (N - 1) / N) * e_s[k] * e_s[l] * N / n * N / n  ))) / N^2

e_s = y_s - X_s %*% beta_hat
sigma_GREG = sqrt(drop(t(e_s)  %*% Omega %*% e_s))

e_s = y_s - X_s_oracle %*% beta_hat_oracle
sigma_GREG_oracle = sqrt(drop(t(e_s)  %*% Omega %*% e_s))

e_s = y_s - X_s %*% beta_hat_Lasso
sigma_Lasso = sqrt(drop(t(e_s)  %*% Omega %*% e_s))

e_s = y_s - X_s[, beta_hat_Lasso != 0, drop = F] %*% beta_hat_refit
sigma_Lasso2 = sqrt(drop(t(e_s)  %*% Omega %*% e_s))

y_res = c(
  HT = y_HT,
  Diff = y_diff,
  GREG = y_GREG,
  GREG_oracle = y_GREG_oracle,
  Lasso = y_Lasso,
  LassoRefit = y_Lasso2
)

sigma_res = c(
  HT = sigma_HT,
  Diff = sigma_diff,
  GREG = sigma_GREG,
  GREG_oracle = sigma_GREG_oracle,
  Lasso = sigma_Lasso,
  LassoRefit = sigma_Lasso2
)

# seq_K = round(c(2^(1:floor(log2(n))), n))
seq_K = round(c(2, n))
# seq_K = round(c(2, N))
# seq_K = 2
for (K in seq_K) {
  # print(K)
  
  y_debiased_vec = NULL
  y_debiased_Lasso_vec = NULL
  y_debiased_Lasso_vec2 = NULL
  e_s_vec = list()
  e_s_vec_Lasso = list()
  e_s_vec_Lasso2 = list()
  # e_s_vec = numeric(n)
  # e_s_vec_Lasso = numeric(n)
  
  shuffled <- sample(1:N)
  
  # Create a grouping factor for partitioning
  groups <- cut(seq_along(shuffled),
                breaks = K,
                labels = FALSE)
  
  # Split the shuffled vector into K partitions
  SubIndex_list = split(shuffled, groups)
  
  for (k in 1:K) {
    # for(k in 1){
    # set.seed(k)
    # SubIndex = sample(1:N, size = round(N / 2), replace = FALSE)
    SubIndex = SubIndex_list[[k]]
    
    Index_sub1 = Index[Index %in% SubIndex]
    Index_sub2 = Index[!(Index %in% SubIndex)]
    # Index_sub1 = Index[sample(1:n, size = m, replace = FALSE)]
    # Index_sub2 = Index[!(Index %in% Index_sub1)]
    
    y_s1 = y[Index_sub1]
    X_s1 = X[Index_sub1, , drop = F]
    d1_s1 = d1[Index_sub1]
    # if(wls){
    #   lm_obj1 = lm(y_s1 ~  0 +X_s1, weights = d1_s1)
    # }else{
    #   lm_obj1 = lm(y_s1 ~  0 +X_s1)
    # }
    # beta_hat1 = lm_obj1$coefficients
    
    y_s2 = y[Index_sub2]
    X_s2 = X[Index_sub2, , drop = F]
    d1_s2 = d1[Index_sub2]
    if (wls) {
      lm_obj2 = lm(y_s2 ~ 0 +  X_s2, weights = d1_s2)
    } else{
      lm_obj2 = lm(y_s2 ~ 0 +  X_s2)
    }
    beta_hat2 = lm_obj2$coefficients
    
    # y_debiased = drop(colSums(X) %*% beta_hat + sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1)
    # + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
    
    # y_debiased = drop(colSums(X) %*% (sum(d1_s2) / sum(d1_s) * beta_hat1 + sum(d1_s1) / sum(d1_s) * beta_hat2) +
    #                     sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
    
    y_debiased = drop(colSums(X[SubIndex, , drop = F]) %*% beta_hat2 + sum((y_s1 - drop(
      X_s1 %*% beta_hat2
    )) * d1_s1))
    
    # y_debiased = drop(colSums(X[SubIndex, ]) %*% beta_hat2 + colSums(X[-SubIndex, ]) %*% beta_hat1 +
    #                     sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
    
    
    # e_s_tmp = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
    e_s_tmp = drop(y_s1 - X_s1 %*% beta_hat2)
    names(e_s_tmp) = Index_sub1
    e_s_vec = append(e_s_vec, list(e_s_tmp))
    # e_s_vec = e_s_vec + e_s_tmp
    
    if (is.na(y_debiased))
      stop()
    
    y_debiased_vec = c(y_debiased_vec, y_debiased)
    
    # if(simnum == 1){
    #   if(wls){
    #     cv_model <- cv.glmnet(X_s1, y_s1, weights = d1_s1)
    #   }else{
    #     cv_model <- cv.glmnet(X_s1, y_s1)
    #   }
    # }
    
    #find optimal lambda value that minimizes test MSE
    # best_lambda1 <- cv_model$lambda.min
    # best_lambda1
    # if(wls){
    #   best_model1 <- glmnet(X_s1, y_s1, weights = d1_s1, lambda = best_lambda1)
    # }else{
    #   best_model1 <- glmnet(X_s1, y_s1, lambda = best_lambda1)
    # }
    # beta_hat1 = as.vector(coef(best_model1))[-1]
    
    if(simnum == 0){
      if(wls){
        cv_model2 <- cv.glmnet(X_s2, y_s2, weights = d1_s2)
      }else{
        cv_model2 <- cv.glmnet(X_s2, y_s2)
      }
    }

    
    #find optimal lambda value that minimizes test MSE
    best_lambda2 <- cv_model2$lambda.min
    best_lambda2
    
    if (wls) {
      best_model2 <- glmnet(X_s2, y_s2, weights = d1_s2, lambda = best_lambda2)
    } else{
      best_model2 <- glmnet(X_s2, y_s2, lambda = best_lambda2)
    }
    beta_hat2 = as.vector(coef(best_model2))[-1]
    
    # y_debiased_Lasso = drop(colSums(X[SubIndex, ]) %*% beta_hat2 + colSums(X[-SubIndex, ]) %*% beta_hat1 +
    #                           sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
    
    y_debiased_Lasso = drop(colSums(X[SubIndex, , drop = F]) %*% beta_hat2 +
                              sum((y_s1 - drop(
                                X_s1 %*% beta_hat2
                              )) * d1_s1))
    
    # e_s_tmp_Lasso = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
    
    e_s_tmp_Lasso = drop(y_s1 - X_s1 %*% beta_hat2)
    names(e_s_tmp_Lasso) = Index_sub1
    e_s_vec_Lasso = append(e_s_vec_Lasso, list(e_s_tmp_Lasso))
    # e_s_vec_Lasso = e_s_vec_Lasso + e_s_tmp_Lasso
    
    y_debiased_Lasso_vec = c(y_debiased_Lasso_vec, y_debiased_Lasso)
    
    X_s_refit2 = X_s2[, beta_hat2 != 0, drop = F]
    if (wls) {
      refit_model = lm(y_s2 ~  0 + X_s_refit2, weights = d1_s2)
    } else{
      refit_model = lm(y_s2 ~  0 + X_s_refit2)
    }
    beta_hat_refit2 = refit_model$coefficients
    y_debiased_Lasso2 = drop(colSums(X[SubIndex, beta_hat2 != 0, drop = F]) %*% beta_hat_refit2 +
                               sum((
                                 y_s1 - drop(X_s1[, beta_hat2 != 0, drop = F] %*% beta_hat_refit2)
                               ) * d1_s1))
    
    e_s_tmp_Lasso2 = drop(y_s1 - X_s1[, beta_hat2 != 0, drop = F] %*% beta_hat_refit2)
    names(e_s_tmp_Lasso2) = Index_sub1
    e_s_vec_Lasso2 = append(e_s_vec_Lasso2, list(e_s_tmp_Lasso2))
    
    y_debiased_Lasso_vec2 = c(y_debiased_Lasso_vec2, y_debiased_Lasso2)
  }
  
  # e_s_tmp = e_s_vec / K
  # e_s_tmp_Lasso = e_s_vec_Lasso / K
  
  y_debiased = sum(y_debiased_vec)
  y_debiased_Lasso = sum(y_debiased_Lasso_vec)
  y_debiased_Lasso2 = sum(y_debiased_Lasso_vec2)
  
  # y_unbiased = drop(colSums(X) %*% beta_hat1 + sum((y_s1 - drop(X_s1 %*% beta_hat1)))
  #                   + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * 2))
  
  # SubIndex_unlist = unlist(SubIndex_list)
  # e_s = unlist(sapply(Index, function(i) e_s_vec[SubIndex_unlist == i]))
  # sigma_Debiased = sqrt(drop(t(e_s_tmp)  %*% Omega %*% e_s_tmp))
  e_s = unlist(e_s_vec)[as.character(Index)]
  sigma_Debiased = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
  
  # cor(y_s - X_s %*% beta_hat, e_s)
  
  # e_s1 = ifelse(Index %in% Index_sub1, y_s - X_s %*% lm_obj2$coefficients, 0)
  # e_s2 = ifelse(Index %in% Index_sub2, y_s - X_s %*% lm_obj1$coefficients, 0)
  # sigma_Debiased = sqrt(drop((t(e_s1)  %*% Omega %*% e_s1 + t(e_s2)  %*% Omega %*% e_s2)))
  
  # sigma_Debiased_Lasso = sqrt(drop(t(e_s_tmp_Lasso)  %*% Omega %*% e_s_tmp_Lasso))
  # e_s = unlist(sapply(Index, function(i) e_s_vec_Lasso[SubIndex_unlist == i]))
  e_s = unlist(e_s_vec_Lasso)[as.character(Index)]
  sigma_Debiased_Lasso = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
  
  e_s = unlist(e_s_vec_Lasso2)[as.character(Index)]
  sigma_Debiased_Lasso2 = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
  
  y_res = c(y_res,
            setNames(c(SS = y_debiased, SSLasso = y_debiased_Lasso, SSLassoRefit = y_debiased_Lasso2), 
                     paste(c("SSGREG", "SSLasso", "SSLassoRefit"), "(", K, ")", sep = "")))
  
  sigma_res = c(sigma_res,
                setNames(c(SS = sigma_Debiased, SSLasso = sigma_Debiased_Lasso, SSLassoRefit = sigma_Debiased_Lasso2), 
                         paste(c("SSGREG", "SSLasso", "SSLassoRefit"), "(", K, ")", sep = "")))
  
}

# sum(beta_hat_Lasso[-c(2 : (1 + s))] != 0)
