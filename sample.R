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