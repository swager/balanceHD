#' Estimate mean outcome at balance.target via residual balancing
#'
#' @param XW the input features for the sub-population of interest
#' @param YW the observed responses for the sub-population of interest
#' @param balance.target the desired center of the dataset
#' @param allow.negative.weights whether negative gammas are allowed for balancing
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fit.method the method used to fit mu(x) = E[YW | XW = x]
#' @param alpha tuning paramter for glmnet
#'
#' @return Estimate for E[YW | XW = balance.target]
#'
#' @export residualBalance.mean
residualBalance.mean = function(XW, YW,
		balance.target,
		allow.negative.weights = FALSE,
		zeta,
		fit.method = c("elnet", "none"),
		alpha) {
	
	fit.method = match.arg(fit.method)
	
	gamma = approx.balance(XW, balance.target, zeta = zeta, allow.negative.weights = allow.negative.weights)
	
	if (fit.method == "elnet") {
		
		lasso.fit = glmnet::cv.glmnet(XW, YW, alpha = alpha)
		mu.lasso = predict(lasso.fit, newx = matrix(balance.target, 1, length(balance.target)))
		
		residuals = YW - predict(lasso.fit, newx = XW)
		mu.residual = sum(gamma * residuals)
		
	} else if (fit.method == "none") {
	
		mu.lasso = 0
		mu.residual = sum(gamma * YW)
	
	} else {
		
		stop("Invalid choice of fitting method.")
	
	}
	
	mu.lasso + mu.residual
}

#' Estimate ATE via approximate residual balancing
#'
#' @param X the input features
#' @param Y the observed responses
#' @param W treatment/control assignment, coded as 0/1
#' @param target.pop which population should the treatment effect be estimated for?
#'            (0, 1): average treatment effect for everyone
#'            0: average treatment effect for controls
#'            1: average treatment effect for treated
#' @param allow.negative.weights whether negative gammas are allowed for balancing
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fit.method the method used to fit mu(x, w) = E[Y | X = x, W = w]
#' @param alpha tuning paramter for glmnet
#' @param scale.X whether non-binary features should be noramlized
#'
#' @return ATE estimate
#'
#' @export residualBalance.ate
residualBalance.ate = function(X, Y, W,
		target.pop=c(0, 1),
		allow.negative.weights = FALSE,
		zeta=0.5,
		fit.method = c("elnet", "none"),
		alpha=0.9,
		scale.X = TRUE) {
	
	fit.method = match.arg(fit.method)
	
	if (scale.X) {
		scl = apply(X, 2, sd, na.rm = TRUE)
		is.binary = apply(X, 2, function(xx) sum(xx == 0) + sum(xx == 1) == length(xx))
		scl[is.binary] = 1
		X.scl = scale(X, center = FALSE, scale = scl)
	} else {
		X.scl = X
	}
	
	# we want ATE for these indices
	target.idx = which(W %in% target.pop)
	balance.target = colMeans(X.scl[target.idx,])
	
	if (setequal(target.pop, c(0, 1))) {
		
		mu0 = residualBalance.mean(X.scl[W==0,], Y[W==0], balance.target, allow.negative.weights, zeta, fit.method, alpha)
		mu1 = residualBalance.mean(X.scl[W==1,], Y[W==1], balance.target, allow.negative.weights, zeta, fit.method, alpha)
		
	} else if (setequal(target.pop, c(1))) {
		
		mu0 = residualBalance.mean(X.scl[W==0,], Y[W==0], balance.target, allow.negative.weights, zeta, fit.method, alpha)
		mu1 = mean(Y[W==1])
		
	} else if (setequal(target.pop, c(0))) {
		
		mu0 = mean(Y[W==0])
		mu1 = residualBalance.mean(X.scl[W==1,], Y[W==1], balance.target, allow.negative.weights, zeta, fit.method, alpha)
		
	} else {
		stop("Invalid target.pop.")
	}

	tau.hat = mu1 - mu0
	tau.hat	
}
