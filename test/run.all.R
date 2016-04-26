rm(list = ls())

library(balanceHD)

# Data generating specs
n = 300
p = 600
tau = 0
nclust = 10
beta = 2 / (1:p) / sqrt(sum(1/(1:p)^2))
clust.ptreat = rep(c(0.1, 0.9), nclust/2)
	
# Generate data
cluster.center = 0.5 * matrix(rnorm(nclust * p), nclust, p)
cluster = sample.int(nclust, n, replace = TRUE)
X = cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, clust.ptreat[cluster])
Y = X %*% beta + rnorm(n, 0, 1) + tau * W

# Run residual balancing
tau.rb.all = residualBalance.ate(X, Y, W, estimate.se = TRUE)
tau.rb = tau.rb.all[1]
tau.rb.negative = residualBalance.ate(X, Y, W, allow.negative.weights = TRUE)
tau.rb.plain = residualBalance.ate(X, Y, W, fit.method = "none")

# Run inverse-propensity weighted methods
tau.ipw.l1 = ipw.ate(X, Y, W, fit.method = "none", alpha.prop = 1)
tau.ipw.elnet = ipw.ate(X, Y, W, fit.method = "none")
tau.ipw.rf = ipw.ate(X, Y, W, fit.method = "none", prop.method = "randomforest")
tau.ipw.resid = ipw.ate(X, Y, W, fit.method = "elnet")

# Run other baselines
tau.naive = naive.ate(X, Y, W)
tau.elnet = elnet.ate(X, Y, W)
tau.twostep = twostep.lasso.ate(X, Y, W)

results = c(Naive=tau.naive,
		Elnet=tau.elnet,
		Twostep = tau.twostep,
		ResidBalance=tau.rb,
		ResidBalanceNegative=tau.rb.negative,
		ResidBalancePlain=tau.rb.plain,
		IPW.Lasso = tau.ipw.l1,
		IPW.Elnet = tau.ipw.elnet,
		IPW.RF = tau.ipw.rf,
		IPW.Residual = tau.ipw.resid)

print("Error tau.hat - tau:")
print(results - tau)

print(paste0("95% CI for tau via approx. residual balancing: (", round(tau.rb.all[1] - 1.96 * tau.rb.all[2], 2), ", ", round(tau.rb.all[1] + 1.96 * tau.rb.all[2], 2), ")"))