#' Estimate ATE via inverse propensity weighting
#'
#' @param X the input features
#' @param Y the observed responses
#' @param W treatment/control assignment, coded as 0/1
#' @param target.pop which population should the treatment effect be estimated for?
#'            (0, 1): average treatment effect for everyone
#'            0: average treatment effect for controls
#'            1: average treatment effect for treated
#' @param eps.threshold cap on the estimated propensities
#' @param fit.method the method used to fit mu(x, w) = E[Y | X = x, W = w]
#' @param alpha.fit tuning paramter for glmnet in the mu model
#' @param prop.method the method used to fit e(x) = P[W = 1 | X]
#' @param alpha.prop tuning paramter for glmnet in the propsenity model
#'
#' @return ATE estimate
#'
#' @export ipw.ate
ipw.ate = function(X, Y, W, target.pop=c(0, 1), eps.threshold = 1/20,
			fit.method = c("elnet", "none"), alpha.fit = 0.9,
			prop.method = c("elnet", "randomforest"), alpha.prop = 0.5,
			estimate.se=FALSE) {
	
	fit.method = match.arg(fit.method)
	prop.method = match.arg(prop.method)
	
	# First compute propensity of treatment
	if (prop.method == "elnet") {
		
		prop.fit = glmnet::cv.glmnet(X, W, family = "binomial", alpha = alpha.prop)
		prop.raw = predict(prop.fit, newx = X, type = "response")
		
	} else if (prop.method == "randomforest") {
		
		if (!requireNamespace("randomForest", quietly = TRUE)) {
    			stop("Package randomForest is needed for this function to work. Please install it.", call. = FALSE)
    		}
		
		prop.fit = randomForest::randomForest(X, factor(W))
		prop.raw = predict(prop.fit, type = "prob")[,2]
		
	} else {
		
		stop("Invalid choice for prop.method.")
		
	}
	
	# make sure no weights can be too large
	
	prop = prop.raw
	if (0 %in% target.pop) {
		prop = pmax(prop, eps.threshold)
	}
	
	if (1 %in% target.pop) {
		prop = pmin(prop, 1 - eps.threshold)
	}
	
	# Next turn this into propensity weights
	if (setequal(target.pop, c(0, 1))) {
		prop.weights = W / prop + (1 - W) / (1 - prop)
	} else if (setequal(target.pop, c(1))) {
		prop.weights = W + (1 - W) * prop / (1 - prop)
	} else if (setequal(target.pop, c(0))) {
		prop.weights = (1 - W) + W * (1 - prop) / prop
	} else {
		stop("Invalid target.pop.")
	}

	# normalize the weights
	prop.weights[W==0] = prop.weights[W==0] / sum(prop.weights[W==0])
	prop.weights[W==1] = prop.weights[W==1] / sum(prop.weights[W==1])
	
	# now fit a model to the outcomes, for W in refit.W
	
	if (setequal(target.pop, c(0, 1))) {
		refit.W = c(0, 1)
	} else if (setequal(target.pop, c(1))) {
		refit.W = c(0)
	} else if (setequal(target.pop, c(0))) {
		refit.W = c(1)
	} else {
		stop("Invalid target.pop.")
	}

	predictions = rep(0, length(Y))
	mu.main = c(0, 0)
	
	for (treatment.status in refit.W) {
		
		W.idx = which(W == treatment.status)
		if (fit.method == "elnet") {
			fit = glmnet::cv.glmnet(X[W.idx,], Y[W.idx], alpha = alpha.fit)
			predictions[W.idx] = predict(fit, newx = X[W.idx,])
			mu.main[treatment.status + 1] = mean(predict(fit, newx = X))
		} else if (fit.method == "none") {
			# do nothing
		} else {
			stop("Invalid choice for fit.method.")
		}
		
	}

	residuals = Y - predictions
	mu.residual = c(sum(prop.weights[W==0] * residuals[W==0]),
			sum(prop.weights[W==1] * residuals[W==1]))
				    
	mu.full = mu.main + mu.residual

	tau.hat = mu.full[2] - mu.full[1]
	
	if(estimate.se == FALSE) {
		tau.hat
	} else {
		var.hat = sum(prop.weights^2 * residuals^2)
		c(tau.hat, sqrt(var.hat))
	}
}
