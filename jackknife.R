N = 2000
n = 500
p = 40
s = 5
# K = 140
wls = T

set.seed(2)

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
library(xtable)
library(mvtnorm)
library(glmnet)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
# suppressMessages(library(doMC))
suppressMessages(library(doRNG))

# Determine the number of CPU cores ####
cores = min(detectCores() - 3, 101)
print(paste("cores =", cores))

# Record the current timestamp. ####
timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

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
mu = X %*% beta
# mu = cbind(exp(1.25 * sin(X[,1:5])), X[,6:ncol(X)]) %*% beta
y = mu + e 
# var(exp(1.25 * sin(X[,2])))

t_y = sum(y)

cv_model <- cv.glmnet(X, y)

# pi1 = e^2 / sum(e^2) * n

BIAS_res = NULL
SE_res = NULL
RMSE_res = NULL
RB_res = NULL
CR_res = NULL
X = cbind(1, X)
beta = c(0, beta)

# Set.seed for the multi-clusters ####
# cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
cl <- makeCluster(cores)
registerDoParallel(cl)

# seq_K = round(c(2^(1:floor(log2(n))), n))
# seq_K = round(c(2, n))
seq_K = 2
for(K in seq_K){
  registerDoRNG(seed = 11)
  print(K)

r = 0.75 # To be changed
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
  
  # for(simnum in 1:SIMNUM){
  #   # print(simnum)
  #   set.seed(simnum)
  final_res <- foreach(
    simnum = 1:SIMNUM, 
    .packages = c("glmnet"), 
    .errorhandling="pass") %dopar% {  
    
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
    
    # if(simnum == 1) plot(1 / d1, e, xlab = "Inclusion Probability", ylab = "error")
    
    # Index = sample(1:N, size = n, replace = FALSE, prob = pi1)
    y_s = y[Index]
    X_s = X[Index, ,drop = F]
    d1_s = d1[Index]
    if(wls){
      lm_obj = lm(y_s ~  0 + X_s, weights = d1_s)
    }else{
      lm_obj = lm(y_s ~  0 + X_s)
    }
    
    beta_hat = lm_obj$coefficients
    #drop(solve(t(X_s) %*% diag(d1_s) %*% X_s, t(X_s) %*% diag(d1_s) %*% y_s))
    
    # beta_hat = unname(lm(y_s ~  t(apply(X_s, 1, function(k) k 
    # - colSums(X_s * d1_s) / N)), weights = d1_s)$coefficients[-1])
    
    y_HT = sum(y_s * d1_s)
    
    y_diff = drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))
    
    y_GREG = drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
    
    # if(wls){
    #   cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
    # }else{
    #   cv_model <- cv.glmnet(X_s, y_s)
    # }
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    
    if(wls){
      best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = best_lambda)
    }else{
      best_model <- glmnet(X_s, y_s, lambda = best_lambda)
    }
    beta_hat_Lasso = as.vector(coef(best_model))[-1]
    
    y_Lasso = drop(colSums(X) %*% beta_hat_Lasso + sum((y_s - drop(X_s %*% beta_hat_Lasso)) * d1_s))
    
    # y_model = drop(colSums(X) %*% beta_hat)
    
    y_debiased_vec = NULL
    y_debiased_Lasso_vec = NULL
    e_s_vec = list()
    e_s_vec_Lasso = list()
    # e_s_vec = numeric(n)
    # e_s_vec_Lasso = numeric(n)
    
    shuffled <- sample(1:N)
    
    # Create a grouping factor for partitioning
    groups <- cut(seq_along(shuffled), breaks = K, labels = FALSE)
    
    # Split the shuffled vector into K partitions
    SubIndex_list = split(shuffled, groups)
    
    for(k in 1:K){
      # for(k in 1){
      # set.seed(k)
      # SubIndex = sample(1:N, size = round(N / 2), replace = FALSE)
      SubIndex = SubIndex_list[[k]]
      
      Index_sub1 = Index[Index %in% SubIndex]
      Index_sub2 = Index[!(Index %in% SubIndex)]
      # Index_sub1 = Index[sample(1:n, size = m, replace = FALSE)]
      # Index_sub2 = Index[!(Index %in% Index_sub1)]
      
      y_s1 = y[Index_sub1]
      X_s1 = X[Index_sub1, ,drop = F]
      d1_s1 = d1[Index_sub1]
      # if(wls){
      #   lm_obj1 = lm(y_s1 ~  0 +X_s1, weights = d1_s1)
      # }else{
      #   lm_obj1 = lm(y_s1 ~  0 +X_s1)
      # }
      # beta_hat1 = lm_obj1$coefficients
      
      y_s2 = y[Index_sub2]
      X_s2 = X[Index_sub2, ,drop = F]
      d1_s2 = d1[Index_sub2]
      if(wls){
        lm_obj2 = lm(y_s2 ~ 0 +  X_s2, weights = d1_s2)
      }else{
        lm_obj2 = lm(y_s2 ~ 0 +  X_s2)
      }
      beta_hat2 = lm_obj2$coefficients
      
      # y_debiased = drop(colSums(X) %*% beta_hat + sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) 
      # + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
      
      # y_debiased = drop(colSums(X) %*% (sum(d1_s2) / sum(d1_s) * beta_hat1 + sum(d1_s1) / sum(d1_s) * beta_hat2) + 
      #                     sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
      
      y_debiased = drop(colSums(X[SubIndex, ,drop = F]) %*% beta_hat2 + sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1))
      
      # y_debiased = drop(colSums(X[SubIndex, ]) %*% beta_hat2 + colSums(X[-SubIndex, ]) %*% beta_hat1 +
      #                     sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
      
      
      # e_s_tmp = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
      e_s_tmp = drop(y_s1 - X_s1 %*% beta_hat2)
      e_s_vec = append(e_s_vec, list(e_s_tmp))
      # e_s_vec = e_s_vec + e_s_tmp
      
      if(is.na(y_debiased)) stop()
      
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
      
      # if(wls){
      #   cv_model <- cv.glmnet(X_s2, y_s2, weights = d1_s2)
      # }else{
      #   cv_model <- cv.glmnet(X_s2, y_s2)
      # }
      
      #find optimal lambda value that minimizes test MSE
      best_lambda2 <- cv_model$lambda.min
      best_lambda2
      
      if(wls){
        best_model2 <- glmnet(X_s2, y_s2, weights = d1_s2, lambda = best_lambda2)
      }else{
        best_model2 <- glmnet(X_s2, y_s2, lambda = best_lambda2)
      }
      beta_hat2 = as.vector(coef(best_model2))[-1]

      # y_debiased_Lasso = drop(colSums(X[SubIndex, ]) %*% beta_hat2 + colSums(X[-SubIndex, ]) %*% beta_hat1 +
      #                           sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1) + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * d1_s2))
            
      y_debiased_Lasso = drop(colSums(X[SubIndex, , drop = F]) %*% beta_hat2 + 
                                sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1))
      
      # e_s_tmp_Lasso = ifelse(Index %in% Index_sub1, y_s - X_s %*% beta_hat2, y_s - X_s %*% beta_hat1)
      
      e_s_tmp_Lasso = drop(y_s1 - X_s1 %*% beta_hat2)
      e_s_vec_Lasso = append(e_s_vec_Lasso, list(e_s_tmp_Lasso))
      # e_s_vec_Lasso = e_s_vec_Lasso + e_s_tmp_Lasso
      
      y_debiased_Lasso_vec = c(y_debiased_Lasso_vec, y_debiased_Lasso)
    }
    
    # e_s_tmp = e_s_vec / K
    # e_s_tmp_Lasso = e_s_vec_Lasso / K
    
    y_debiased = sum(y_debiased_vec)
    
    y_debiased_Lasso = sum(y_debiased_Lasso_vec)
    
    # y_unbiased = drop(colSums(X) %*% beta_hat1 + sum((y_s1 - drop(X_s1 %*% beta_hat1))) 
    #                   + sum((y_s2 - drop(X_s2 %*% beta_hat1)) * 2))
    
    
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
    
    e_s = y_s - X_s %*% beta_hat_Lasso
    sigma_Lasso = sqrt(drop(t(e_s)  %*% Omega %*% e_s))
    
    # sigma_Debiased = sqrt(drop(t(e_s_tmp)  %*% Omega %*% e_s_tmp))
    sigma_Debiased = sqrt(sum(mapply(function(i, e){t(e) %*% 
        Omega[which(Index %in% i), which(Index %in% i)] %*% e}, 
           SubIndex_list, e_s_vec)))
    
    # e_s1 = ifelse(Index %in% Index_sub1, y_s - X_s %*% lm_obj2$coefficients, 0)
    # e_s2 = ifelse(Index %in% Index_sub2, y_s - X_s %*% lm_obj1$coefficients, 0)
    # sigma_Debiased = sqrt(drop((t(e_s1)  %*% Omega %*% e_s1 + t(e_s2)  %*% Omega %*% e_s2)))
  
    # sigma_Debiased_Lasso = sqrt(drop(t(e_s_tmp_Lasso)  %*% Omega %*% e_s_tmp_Lasso))
    sigma_Debiased_Lasso = sqrt(sum(mapply(function(i, e){t(e) %*% 
        Omega[which(Index %in% i), which(Index %in% i)] %*% e}, 
        SubIndex_list, e_s_vec_Lasso)))
    
    y_res = c(HT = y_HT, Diff = y_diff, GREG = y_GREG, Lasso = y_Lasso, 
              SS = y_debiased, SSLasso = y_debiased_Lasso)
    
    # res = rbind(res, y_res)
    
    sigma_res = c(HT = sigma_HT, Diff = sigma_diff, GREG = sigma_GREG, Lasso = sigma_Lasso,
                  SS = sigma_Debiased, SSLasso = sigma_Debiased_Lasso)
    
    # res2 = rbind(res2, sigma_res)
    
    # res3 = rbind(res3, ifelse(abs(y_res - t_y) > 1.96 * sigma_res, 0, 1))
    
    list(res = y_res, 
         res2 = sigma_res,
         res3 = ifelse(abs(y_res - t_y) > 1.96 * sigma_res, 0, 1))
  }
  final_res1 = lapply(final_res, function(x) x[[1]])
  final_res2 = lapply(final_res, function(x) x[[2]])
  final_res3 = lapply(final_res, function(x) x[[3]])
  
  print(paste("# of failure:", sum(!sapply(final_res1, function(x) is.numeric(unlist(x)))))); final_res0 = final_res1
  final_res1 = final_res1[sapply(final_res1, function(x) is.numeric(unlist(x)))]
  res = do.call("rbind", final_res1)
  res2 = do.call("rbind", final_res2)
  res3 = do.call("rbind", final_res3)
  
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
  
  SE_res = cbind(SE_res, SE)
  RMSE_res = cbind(RMSE_res, RMSE)
  BIAS_res = cbind(BIAS_res, BIAS)
  RB_res = cbind(RB_res, BIAS2 / SE)
  CR_res = cbind(CR_res, colMeans(res3))
}
  stopCluster(cl)
  timenow2 = Sys.time()
  print("Running time")
  print(timenow2 - timenow1)

xtable(cbind(BIAS = BIAS_res[,1], SE = SE_res[,1], RMSE = RMSE_res[,1]) )
xtable(cbind(BIAS = BIAS_res[,length(seq_K)], SE = SE_res[,length(seq_K)], RMSE = RMSE_res[,length(seq_K)]) )

xtable(cbind(RB = RB_res[,1], CR = CR_res[,1]) )

matplot(t(RMSE_res[-c(1, 2),]), type = "l", col = hcl.colors(4, "Temps"), lty = 1, lwd = 2, ylim = c(min(SE_res[-c(1, 2),]), max(RMSE_res[-c(1, 2),])),
        xlab = "K", xaxt = "n", ylab = "", main = "solid line = RMSE, dashed line = SE")
axis(1, at = seq(seq_K), labels = seq_K, cex.axis = 0.7)
matlines(t(SE_res[-c(1, 2),]), type = "l", col = hcl.colors(4, "Temps"), lty = 2, lwd = 2)
legend("topright", rownames(RMSE_res[-c(1, 2),]), col = hcl.colors(4, "Temps"), lty = 1, cex = 0.7)


