# pi = pt(- e - 2, 3);
pi = pt(- e , 3);
pi = ifelse(pi >.7, .7, pi)

mean(pi)
delta = rbinom(N, 1, pi)
Index = which(delta == 1)
n = length(Index); #print(n)
d1 = 1 / pi

y_s = y[Index]
X_s = X[Index, , drop = F]
d1_s = d1[Index]
Omega = diag(d1_s^2 - d1_s)