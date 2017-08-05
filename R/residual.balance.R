#' Estimate mean outcome at balance.target via residual balancing
#'
#' @param XW the input features for the sub-population of interest
#' @param YW the observed responses for the sub-population of interest
#' @param balance.target the desired center of the dataset
#' @param allow.negative.weights whether negative gammas are allowed for balancing
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fit.method the method used to fit mu(x) = E[YW | XW = x]
#' @param alpha tuning paramter for glmnet
#' @param optimizer which optimizer to use for approximate balancing
#' @param use.dual whether balancing should be solved in dual form
#' @param verbose whether the optimizer should print progress information
#'
#' @return Estimate for E[YW | XW = balance.target], along with variance estimate
#'
#' @export residualBalance.mean
residualBalance.mean = function(XW, YW,
		balance.target,
		allow.negative.weights = FALSE,
		zeta,
		fit.method = c("elnet", "none"),
		alpha,
		optimizer = c("mosek", "pogs", "quadprog"),
		use.dual = NULL,
		verbose = FALSE) {
	
	fit.method = match.arg(fit.method)
	optimizer = match.arg(optimizer)
	if(is.null(use.dual)) {
	  use.dual = (optimizer %in% c("mosek"))
	}
	
	gamma = approx.balance(XW, balance.target, zeta = zeta, allow.negative.weights = allow.negative.weights, optimizer = optimizer, use.dual = use.dual, verbose=verbose)
	
	if (fit.method == "elnet") {
		
		lasso.fit = glmnet::cv.glmnet(XW, YW, alpha = alpha)
		mu.lasso = predict(lasso.fit, newx = matrix(balance.target, 1, length(balance.target)))
		
		residuals = YW - predict(lasso.fit, newx = XW)
		mu.residual = sum(gamma * residuals)
		
		var.hat = sum(gamma^2 * residuals^2) *
                # degrees of freedom correction
                  length(gamma) / max(1, length(gamma) - sum(coef(lasso.fit) != 0))
		
	} else if (fit.method == "none") {
	
		mu.lasso = 0
		mu.residual = sum(gamma * YW)
		
		var.hat = NA
	
	} else {
		
		stop("Invalid choice of fitting method.")
	
	}
	
	mu.hat = mu.lasso + mu.residual
	c(mu.hat, var.hat)
}


residualBalance.estimate.var = function(XW, YW, alpha, estimate.se) {

  # Don't waste time
  if (!estimate.se) return(NA)
  
  lasso.fit = glmnet::cv.glmnet(XW, YW, alpha = alpha)
  residuals = YW - predict(lasso.fit, newx = XW)
  mean(residuals^2) / max(1, length(YW) - sum(coef(lasso.fit) != 0))
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
#' @param estimate.se whether to return estimate of standard error
#' @param optimizer which optimizer to use for approximate balancing
#' @param use.dual whether balancing should be solved in dual form
#' @param verbose whether the optimizer should print progress information
#'
#' @return ATE estimate, along with (optional) standard error estimate
#'
#' @export residualBalance.ate
residualBalance.ate = function(X, Y, W,
		target.pop=c(0, 1),
		allow.negative.weights = FALSE,
		zeta=0.5,
		fit.method = c("elnet", "none"),
		alpha=0.9,
		scale.X = TRUE,
		estimate.se = FALSE,
		optimizer = c("mosek", "pogs", "quadprog"),
		use.dual = NULL,
		verbose = FALSE) {
	
	fit.method = match.arg(fit.method)
	optimizer = match.arg(optimizer)
	if(is.null(use.dual)) {
	  use.dual = (optimizer %in% c("mosek"))
	}
	
	if (estimate.se & fit.method == "none") {
		warning("Cannot estimate standard error with fit.method = none. Forcing estimate.se to FALSE.")
		estimate.se = FALSE
	}
	
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
		
		est0 = residualBalance.mean(X.scl[W==0,], Y[W==0], balance.target, allow.negative.weights, zeta, fit.method, alpha, optimizer=optimizer, use.dual = use.dual, verbose=verbose)
		est1 = residualBalance.mean(X.scl[W==1,], Y[W==1], balance.target, allow.negative.weights, zeta, fit.method, alpha, optimizer=optimizer, use.dual = use.dual, verbose=verbose)
		
	} else if (setequal(target.pop, c(1))) {
		
		est0 = residualBalance.mean(X.scl[W==0,], Y[W==0], balance.target, allow.negative.weights, zeta, fit.method, alpha, optimizer=optimizer, use.dual = use.dual, verbose=verbose)
		est1 = c(mean(Y[W==1]), residualBalance.estimate.var(X[W==1,], Y[W==1], alpha, estimate.se))
		
	} else if (setequal(target.pop, c(0))) {
		
		est0 = c(mean(Y[W==0]), residualBalance.estimate.var(X[W==0,], Y[W==0], alpha, estimate.se))
		est1 = residualBalance.mean(X.scl[W==1,], Y[W==1], balance.target, allow.negative.weights, zeta, fit.method, alpha, optimizer=optimizer, use.dual = use.dual, verbose=verbose)
		
	} else {
		
		stop("Invalid target.pop.")
		
	}

	tau.hat = est1[1] - est0[1]
	var.hat = est1[2] + est0[2]
	
	if (estimate.se) {
		return(c(tau.hat, sqrt(var.hat)))
	} else {
		return(tau.hat)
	}
}
