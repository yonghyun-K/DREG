N = 1000
N1 = 500
n = 20

set.seed(1)
y = rnorm(N)
z = y + rnorm(N)

res = NULL


for(simnum in 1:30000){
  Index_U1 = sample(1:N, N1)
  Index_U2 = (1:N)[-Index_U1]
  
  set.seed(simnum)
  
  Index_A = sample(1:N, n)
  
  Index1 = Index_A[Index_A %in% Index_U1]
  Index2 = Index_A[Index_A %in% Index_U2]
  
  # if(length(Index1) == 0 | length(Index2) == 0) stop()
  
  y1 = ifelse(length(Index1) == 0, 0, mean(y[Index1]))
  z2 = ifelse(length(Index2) == 0, 0, mean(z[Index2]))
  
  res = rbind(res, c(y1, z2))
}
cov(res)
-cov(y, z) / N
# -var(y) / N


