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
#' @param prop.weighted.fit whether propensity weights should be used to as sample weights in outcome fit
#' @param targeting.method how to combine the outcome and propensity model fits.
#'
#' @return ATE estimate
#'
#' @export ipw.ate
ipw.ate = function(X, Y, W, target.pop=c(0, 1), eps.threshold = 1/20,
			fit.method = c("elnet", "none"), alpha.fit = 0.9,
			prop.method = c("elnet", "randomforest"), alpha.prop = 0.5,
			prop.weighted.fit = FALSE,
			targeting.method = c("AIPW", "TMLE"),
			estimate.se=FALSE) {
	
	fit.method = match.arg(fit.method)
	prop.method = match.arg(prop.method)
	targeting.method = match.arg(targeting.method)
	
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
	if (estimate.se) {
	  refit.W = c(0, 1)
	} else if (setequal(target.pop, c(0, 1))) {
		refit.W = c(0, 1)
	} else if (setequal(target.pop, c(1))) {
		refit.W = c(0)
	} else if (setequal(target.pop, c(0))) {
		refit.W = c(1)
	} else {
		stop("Invalid target.pop.")
	}

	mu.main = c(mean(Y[W == 0]),  mean(Y[W == 1]))
	predictions = mu.main[W + 1]
	target.idx = which(W %in% target.pop)
	
	if (fit.method != "none") {
	  for (treatment.status in refit.W) {
	    W.idx = which(W == treatment.status)
	    if (fit.method == "elnet") {
	      if (prop.weighted.fit) {
	        sample.weights = prop.weights[W.idx] / mean(prop.weights[W.idx])
	      } else {
	        sample.weights = rep(1, length(W.idx))
	      }
	      fit = glmnet::cv.glmnet(X[W.idx,], Y[W.idx], alpha = alpha.fit, weights = sample.weights)
	      predictions[W.idx] = predict(fit, newx = X[W.idx,])
	      mu.main[treatment.status + 1] = mean(predict(fit, newx = X[target.idx,]))
	    } else {
	      stop("Invalid choice for fit.method.")
	    }
	  }
	}

	residuals = Y - predictions
	
	if (targeting.method == "AIPW") {
	  
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
	  
	} else if (targeting.method == "TMLE") {
	  
	  # Next turn this into propensity weights
	  if (setequal(target.pop, c(0, 1))) {
	    mean.weight.0 = mean(1/(1 - prop)) / mean(1/(1 - prop[W==0]))
	    mean.weight.1 = mean(1/prop) / mean(1/(prop[W==1]))
	  } else if (setequal(target.pop, c(1))) {
	    mean.weight.0 = mean(prop/(1 - prop)) / mean(prop[W==0]/(1 - prop[W==0]))
	    mean.weight.1 = 1
	  } else if (setequal(target.pop, c(0))) {
	    mean.weight.0 = 1
	    mean.weight.1 = mean((1 - prop)/prop) / mean((1 - prop[W==1])/prop[W==1])
	  } else {
	    stop("Invalid target.pop.")
	  }
	  
	  # targeted MLE
	  eps.tmle.0 = lm(B ~ A + 0, data=data.frame(A=prop.weights[W==0] * sum(W==0), B=residuals[W==0]))
	  eps.tmle.1 = lm(B ~ A + 0, data=data.frame(A=prop.weights[W==1] * sum(W==1), B=residuals[W==1]))
	  delta.tmle.0 = predict(eps.tmle.0, newdata=data.frame(A=mean.weight.0))
	  delta.tmle.1 = predict(eps.tmle.1, newdata=data.frame(A=mean.weight.1))
	  tau.tmle = mu.main[2] - mu.main[1] + delta.tmle.1 - delta.tmle.0
	  
	  if(estimate.se == FALSE) {
	    tau.tmle
	  } else {
	    sigma2.tmle = sandwich::vcovHC(eps.tmle.0) * mean.weight.0^2 +
	      sandwich::vcovHC(eps.tmle.1) * mean.weight.1^2
	    c(tau.tmle, sqrt(sigma2.tmle))
	  }
	  
	} else {
	  stop("Invalid targeting method.")
	}
}
