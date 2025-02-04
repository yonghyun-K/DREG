N = 100
n = 20
m = 10
set.seed(1)
y = rnorm(N)

res = NULL
for(simnum in 1:10000){
  set.seed(simnum)
  Index = sample(1:N, n)
  Index1 = sample(Index, m)
  Index2 = Index[!(Index %in% Index1)]
  
  y1 = mean(y[Index1])
  y2 = mean(y[Index2])
  res = rbind(res, c(y1, y2))
}
cov(res)
-var(y) / N


