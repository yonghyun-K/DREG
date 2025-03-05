N = 2000
n = 300
# p = ?
s = 5
SIMNUM = 50
# K = ?
# wls = T

# r = 0 # Noninformative
r = 0.75 # Informative

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
library(xtable)
library(mvtnorm)
library(glmnet)
library(ggplot2)
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

isInteractive = interactive()

if(!isInteractive){
  dir.create(timenow0)
  setwd(timenow0)
  
  sink(timenow, append=TRUE)
}

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(
    1:n - 1,
    nrow = n,
    ncol = n,
    byrow = TRUE
  ) -
    (1:n - 1))
  rho ^ exponent
}

BIAS_res = NULL
SE_res = NULL
RMSE_res = NULL
RB_res = NULL
CR_res = NULL
FDR_res = NULL
FNR_res = NULL
BIASapprox_res = NULL
# seq_p = c(10, 20, 30, 40, 50, 70, 80, 90, 100, 110, 120)
# seq_p = c(10, 20, 40, 80, 120)
seq_p = c(10, 120)

# seq_K = round(c(2, 200)) # K = N is jackknife
# seq_K = round(c(N)) # K = N is jackknife
seq_K = round(c(2)) # K = N is jackknife

# Set.seed for the multi-clusters ####
# cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
cl <- makeCluster(cores)
registerDoParallel(cl)

print(paste("N =", N))
print(paste("n =", n))
print(paste("s =", s))
print(paste("SIMNUM =", SIMNUM))
print(paste("r =", r))
print(paste("seq_K =", paste(seq_K)))

set.seed(2)
X0 = rmvnorm(n = N, rep(0, seq_p[length(seq_p)]), ar1_cor(seq_p[length(seq_p)], 0.2)) + 2
e = rnorm(N, 0, 1)

if (r == 0) {
  z = rnorm(n = N, 0, 1)
} else{
  z = rnorm(n = N, 0, sqrt((1 - r) / r)) + e
}
# z = rnorm(n = N, 0, sqrt(ratio)) + e
# z = e
# z = rnorm(n = N, 0, 1)
z_sorted_idx = order(z)
len = round(N / 4)

for(p in seq_p){
print(paste("p =", p))
X = X0[,1:p]
# X = pnorm(X)

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)

# X = X * c(rep(1, N), rep(-1, N))

# X = matrix(rnorm(N * p, 2, 1), nr = N, nc = p)
# X = scale(X, T, T)

beta = c(rep(1, s), rep(0, p - s))
mu = X %*% beta; if(p == seq_p[1]) print("linear model")
# mu = cbind(exp(1.25 * sin(X[,1:s])), X[,(s+1):ncol(X)]) %*% beta; if(p == seq_p[1]) print("nonlinear model")
# mu = cbind(X[,1:s]^(1:s %% 2 + 1) / 5, X[,(s+1):ncol(X)]) %*% beta; if(p == seq_p[1]) print("nonlinear model")
# mu = cbind(X[,1:s]^(1:s %% 3 + 1) / 20, X[,(s+1):ncol(X)]) %*% beta; if(p == seq_p[1]) print("nonlinear model")
# diag(var(X[,1:s]^(1:s %% 3 + 1) / 20))
y = mu + e
var(mu); var(e)
# var(exp(1.25 * sin(X[,2])))

t_y = sum(y)

# pi1 = e^2 / sum(e^2) * n

X = cbind(1, X)
beta = c(0, beta)

# cv_model <- cv.glmnet(X, y)
# best_lambda <- cv_model$lambda.min
# best_model <- glmnet(X, y, lambda = best_lambda)
# beta_hat_Lasso = as.vector(coef(best_model))[-1]

# pi1 = rep(n / N, N)
# pi1 = rep(n / N, N) + rnorm(N, 0, 0.0025)
# pi1 = rep(n / N, N) + e * n / N * 0.2

# pi1 = (y - min(y)) / (max(y) - min(y)) / 2 + 0.25
# pi1 = (e - min(e)) / (max(e) - min(e)) / 2 + 0.10
# pi1 = pi1 / sum(pi1) * n
# d1 = 1 / pi1

simnum = 0

wls = T
if(!isInteractive){
  # source("../sample.R")
  source("../sample_poi.R", local = TRUE)
  source("../run.R")
}else{
  # source("sample.R")
  source("sample_poi.R", local = TRUE)
  source("run.R")
}

registerDoRNG(seed = 11)
# for(simnum in 1:SIMNUM){
#   # print(simnum)
#   set.seed(simnum)
final_res <- foreach(
  simnum = 1:SIMNUM,
  .export = c("n", "N", "len", "z_sorted_idx", "isInteractive",
              "y", "X", "mu", "cv_model", "cv_model2", "seq_K", "e"),
  .packages = c("glmnet"),
  .errorhandling = "pass"
) %dopar% {
  
  wls = F
  if(!isInteractive){
    # source("../sample.R", local = TRUE)
    source("../sample_poi.R", local = TRUE)
    source("../run.R", local = TRUE)
  }else{
    # source("sample.R", local = TRUE)
    source("sample_poi.R", local = TRUE)
    source("run.R", local = TRUE)
  }
  
  # res2 = rbind(res2, sigma_res)
  
  # res3 = rbind(res3, ifelse(abs(y_res - t_y) > 1.96 * sigma_res, 0, 1))
  
  list_tmp = list(
    res = y_res,
    res2 = sigma_res,
    res3 = ifelse(abs(y_res - t_y) > 1.96 * sigma_res, 0, 1),
    res4 = sum(beta_hat_Lasso[-c(2 : (1 + s))] != 0) / (p + 1 - s), # FDR
    res5 = sum(beta_hat_Lasso[c(2 : (1 + s))] == 0) / (s), # FNR
    res6 = bias_res
  )
  
  wls = T
  if(!isInteractive){
    source("../run.R", local = TRUE)
  }else{
    source("run.R", local = TRUE)
  }

  list(
    res = c(list_tmp[[1]], y_res[-(1:2)]),
    res2 = c(list_tmp[[2]], sigma_res[-(1:2)]),
    res3 = c(list_tmp[[3]], ifelse(abs(y_res - t_y)[-(1:2)] > 1.96 * sigma_res[-(1:2)], 0, 1)),
    res4 = sum(beta_hat_Lasso[-c(2 : (1 + s))] != 0) / (p + 1 - s), # FDR
    res5 = sum(beta_hat_Lasso[c(2 : (1 + s))] == 0) / (s), # FNR
    res6 = c(list_tmp[[6]], bias_res[-(1:2)])
  )
}



final_res1 = lapply(final_res, function(x)
  x[[1]])
final_res2 = lapply(final_res, function(x)
  x[[2]])
final_res3 = lapply(final_res, function(x)
  x[[3]])
final_res4 = lapply(final_res, function(x)
  x[[4]])
final_res5 = lapply(final_res, function(x)
  x[[5]])
final_res6 = lapply(final_res, function(x)
  x[[6]])

print(paste("# of failure:", sum(
  !sapply(final_res1, function(x)
    is.numeric(unlist(x)))
)))
final_res0 = final_res1
final_res1 = final_res1[sapply(final_res1, function(x)
  is.numeric(unlist(x)))]
res = do.call("rbind", final_res1)
res2 = do.call("rbind", final_res2)
res3 = do.call("rbind", final_res3)
res4 = do.call("rbind", final_res4)
res5 = do.call("rbind", final_res5)
res6 = do.call("rbind", final_res6)

if(p == seq_p[1]) resk1 = res

BIAS = colMeans(res - t_y)
SE = apply(res, 2, function(x)
  sqrt(var(x) * (length(x) - 1) / length(x)))
RMSE = apply(res - t_y, 2, function(x)
  sqrt(mean(x ^ 2)))

# colMeans(res)
tmpdf21 = cbind(BIAS, SE, RMSE)
xtable(tmpdf21, digits = 3, caption = "Summary of point estimation")

BIAS2 = colMeans(res2 ^ 2) - SE ^ 2
SE2 = apply(res2, 2, function(x)
  sqrt(var(x) * (length(x) - 1) / length(x)))
RMSE2 = apply(res2 - SE, 2, function(x)
  sqrt(mean(x ^ 2)))

# cbind(BIAS2, REl_BIAS = BIAS2 / SE)
tmpdf22 = cbind(RB = BIAS2 / SE^2, CR = colMeans(res3))
print(cbind(tmpdf21, tmpdf22))
xtable(tmpdf22,
       digits = c(0, 4, 3),
       caption = "Summary of variance estimation")

print("FD summary")
print(summary(res4 * (p + 1 -s)))

print("FN summary")
print(summary(res5 * s))

# colnames(resk1)[6:8] <- paste(colnames(resk1)[6:8], "(K=", seq_K[1], ")", sep = "")
# colnames(res)[6:8] <- paste(colnames(res)[6:8], "(K=", seq_K[2], ")", sep = "")
# boxplot(res[,c(1,2,3,5,11,13,15,17,23,25)], col = rep(c(3,4,5), times = c(2,4,4)), main = paste("p =", p),
#         cex.axis = 0.5)
boxplot(res, main = paste("p =", p),
        cex.axis = 0.5)
abline(h = t_y, lty = 1, col = 2)
abline(v = c(2.5, 6.5), lty = 3)

SE_res = cbind(SE_res, SE)
RMSE_res = cbind(RMSE_res, RMSE)
BIAS_res = cbind(BIAS_res, BIAS)
RB_res = cbind(RB_res, BIAS2 / SE^2)
CR_res = cbind(CR_res, colMeans(res3))
FDR_res = cbind(FDR_res, res4)
FNR_res = cbind(FNR_res, res5)
BIASapprox_res = cbind(BIASapprox_res, colMeans(res6))
}

stopCluster(cl)
timenow2 = Sys.time()
print("Running time")
print(timenow2 - timenow1)

if(!interactive()) save.image(paste(timenow0, ".RData", sep = ""))

# xtable(cbind(BIAS = BIAS_res[,1], SE = SE_res[,1], RMSE = RMSE_res[,1]) )
# xtable(cbind(BIAS = BIAS_res[,length(seq_p)], SE = SE_res[,length(seq_p)], RMSE = RMSE_res[,length(seq_p)]) )

xtable(cbind(BIAS = BIAS_res[,1], BIAS_aprx = BIASapprox_res[,1], SE = SE_res[,1], RMSE = RMSE_res[,1],
             BIAS = BIAS_res[,length(seq_p)], BIAS_aprx = BIASapprox_res[,length(seq_p)], SE = SE_res[,length(seq_p)], RMSE = RMSE_res[,length(seq_p)]) )

# xtable(cbind(RB = RB_res[,1], CR = CR_res[,1]) )
# xtable(cbind(RB = RB_res[,length(seq_p)], CR = CR_res[,length(seq_p)]) )

xtable(cbind(RB = RB_res[,1], CR = CR_res[,1],
             RB = RB_res[,length(seq_p)], CR = CR_res[,length(seq_p)]))

colnames(FDR_res) <- seq_p

if(!isInteractive){
png("boxplot_FDR.png", width = 960, height = 560)
boxplot(FDR_res, xlab = "p", ylab = "FDR")
dev.off()
png("boxplot_FD.png", width = 960, height = 560)
boxplot(FDR_res * rep((seq_p + 1-s), each = nrow(FDR_res)), xlab = "p", ylab = "False Discoveries")
dev.off()
png("boxplot_FNR.png", width = 960, height = 560)
boxplot(FNR_res, xlab = "p", ylab = "FNR")
dev.off()
png("boxplot_FN.png", width = 960, height = 560)
boxplot(FNR_res * rep(s, each = nrow(FNR_res)), xlab = "p", ylab = "False Exclusions")
dev.off()

png("RMSE_linegraph.png", width = 960, height = 560)
includeidx = c(3,5,11,13)
# includeidx = c(3,4,5,9,10,11)
# includeidx = c(9,10,11)
matplot(t(RMSE_res[includeidx,]), type = "l", col = hcl.colors(length(includeidx), "Temps"), lty = 1, lwd = 2, ylim = c(min(SE_res[includeidx,]), max(RMSE_res[includeidx,])),
        xlab = "p", xaxt = "n", ylab = "", main = "solid line = RMSE, dashed line = SE")
axis(1, at = seq(seq_p), labels = seq_p, cex.axis = 0.7)
matlines(t(SE_res[includeidx,]), type = "l", col = hcl.colors(length(includeidx), "Temps"), lty = 2, lwd = 2)
legend("topleft", rownames(RMSE_res[includeidx,]), col = hcl.colors(length(includeidx), "Temps"), lty = 1, cex = 0.7)
dev.off()
}

# Assuming res, t_y, and res6 are already defined
df <- data.frame(
  x = c((res - t_y)[,3], (res - t_y)[,5], (res - t_y)[,7], (res - t_y)[,9]),
  y = c(res6[,3], res6[,5], res6[,7], res6[,9]),
  group = factor(rep(c(3, 5, 7, 9), each = nrow(res6)))
)
colors <- c("black", "red", "blue", "green")
names(colors) <- c("3", "5", "7", "9")
ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point() +
  scale_color_manual(values = colors, labels = rownames(RMSE_res[c(3,5,7,9),])) +
  theme_minimal() +
  geom_abline()+
  labs(color = "Legend", x = "Bias", y = "Bias_approx")

if(!interactive()) save.image(paste(timenow0, ".RData", sep = ""))

# Catour plot for the future
# funct <- function(x, y) (y + 0.001) / (x + 0.001)
# 
# grid <- expand.grid(x = seq(0,1,0.1), y = seq(0,1,0.1))
# 
# ggplot(grid, aes(x, y, z = funct(x,y))) +
#   geom_contour_filled() +  # Use filled contours
#   # scale_fill_viridis_c() +  # Use a perceptually uniform color scale
#   coord_cartesian(xlim = range(grid$x) / 10, ylim = range(grid$y)) +  # Keep xlim and ylim the same
#   theme_minimal()
