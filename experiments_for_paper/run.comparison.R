run.comparison = function(XYW) {

	X = XYW$X
	Y = XYW$Y
	W = XYW$W

	# Run residual balancing
	tau.rb = residualBalance.ate(X, Y, W, target.pop = 1, fit.method = "elnet", alpha = 0.9, zeta = 0.5)
	tau.rb.plain = residualBalance.ate(X, Y, W, target.pop = 1, fit.method = "none", zeta = 0.5)

	# Run inverse-propensity weighted methods
	tau.ipw.elnet = ipw.ate(X, Y, W, target.pop = 1, fit.method = "none", prop.method = "elnet", alpha.prop = 0.5)
	tau.aipw = ipw.ate(X, Y, W, target.pop = 1, prop.weighted.fit = FALSE, targeting.method = "AIPW",
	                   fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)
	tau.ipw.weighted = ipw.ate(X, Y, W, target.pop = 1, prop.weighted.fit = TRUE, targeting.method = "AIPW",
	                           fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)
	tau.tmle = ipw.ate(X, Y, W, prop.weighted.fit = FALSE, targeting.method = "TMLE",
	                   target.pop = 1, fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)

	# Run other baselines
	tau.naive = naive.ate(X, Y, W)
	tau.elnet = elnet.ate(X, Y, W, target.pop = 1, alpha = 0.9)
	tau.twostep = twostep.lasso.ate(X, Y, W, target.pop = 1)

	results = c(Naive=tau.naive,
		Elnet=tau.elnet,
		Twostep = tau.twostep,
		ResidBalance=tau.rb,
		ResidBalancePlain=tau.rb.plain,
		IPW.Elnet = tau.ipw.elnet,
		IPW.Residual = tau.aipw,
		IPW.Weighted = tau.ipw.weighted,
		TMLE = tau.tmle)
	results

}

coverage.test = function(XYW) {

	X = XYW$X
	Y = XYW$Y
	W = XYW$W

	# Run residual balancing
	tau.rb = residualBalance.ate(X, Y, W, target.pop = 1, fit.method = "elnet", alpha = 0.9, zeta = 0.5, estimate.se = TRUE)
	tau.rb
}
