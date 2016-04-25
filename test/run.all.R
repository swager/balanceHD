rm(list = ls())

library(balanceHD)

# Data generating specs
n = 300
p = 800
eps = 0.2
sigma = 1/10
beta = c(rep(10, 10), rep(1, 90), rep(0, p - 100))
delta.clust = rep(1/5, p)
tau = 0

# Generate data
CLUST = rbinom(n, 1, 0.5)
W = rbinom(n, 1, eps + (1 - 2*eps) * CLUST)	
X = matrix(rnorm(n * p), n, p) + outer(CLUST, delta.clust)
Y = X %*% beta + rnorm(n, 0, sigma) + tau * W

# Run residual balancing
tau.rb = residualBalance.ate(X, Y, W, fit.method = "elnet", alpha = 0.9)
tau.rb.negative = residualBalance.ate(X, Y, W, fit.method = "elnet", alpha = 0.9, allow.negative.weights = TRUE)
tau.rb.plain = residualBalance.ate(X, Y, W, fit.method = "none")

# Run inverse-propensity weighted methods
tau.ipw.l1 = ipw.ate(X, Y, W, fit.method = "none", prop.method = "elnet", alpha.prop = 1)
tau.ipw.elnet = ipw.ate(X, Y, W, fit.method = "none", prop.method = "elnet", alpha.prop = 0.5)
tau.ipw.rf = ipw.ate(X, Y, W, fit.method = "none", prop.method = "randomforest")
tau.ipw.resid = ipw.ate(X, Y, W, fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.5, alpha.prop = 0.5)

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

print(results - tau)
