#' No correction for confounding variables (difference in means)
#'
#' @param X the input features
#' @param Y the observed responses
#' @param W treatment/control assignment, coded as 0/1
#'
#' @return ATE estimate
#'
#' @export naive.ate
naive.ate = function(X, Y, W, estimate.se = FALSE) {
	tau.hat = mean(Y[W==1]) - mean(Y[W==0])
	
	if (!estimate.se) {
		return(tau.hat)
	} else {	
		se.hat = sqrt(var(Y[W==0]) / (sum(W==0) - 1) + var(Y[W==1]) / (sum(W==1) - 1))
		return(c(tau.hat, se.hat))
	}
}

#' Elastic net correction for confounding variables
#'
#' @param X the input features
#' @param Y the observed responses
#' @param W treatment/control assignment, coded as 0/1
#' @param target.pop which population should the treatment effect be estimated for?
#'            (0, 1): average treatment effect for everyone
#'            0: average treatment effect for controls
#'            1: average treatment effect for treated
#'
#' @return ATE estimate
#'
#' @export elnet.ate
elnet.ate = function(X, Y, W, target.pop=c(0, 1), alpha = 0.9, estimate.se = FALSE) {

	# we want ATE for these indices
	target.idx = which(W %in% target.pop)
	
	fit0 = glmnet::cv.glmnet(X[W == 0,], Y[W == 0], alpha = alpha)
	fit1 = glmnet::cv.glmnet(X[W == 1,], Y[W == 1], alpha = alpha)
	
	pred0 = coef(fit0)[1] + X %*% coef(fit0)[-1]
	pred1 = coef(fit1)[1] + X %*% coef(fit1)[-1]
	
	tau.hat = mean((pred1 - pred0)[target.idx])
	
	if (!estimate.se) {
		return(tau.hat)
	} else {	
		var0 = var((Y - pred0)[W==0]) / (sum(W==0) - 1)
		var1 = var((Y - pred1)[W==1]) / (sum(W==1) - 1)
		se.hat = sqrt(var0 + var1)
		return(c(tau.hat, se.hat))
	}
}

#' Refit lasso correction for confounding variables
#'
#' @param X the input features
#' @param Y the observed responses
#' @param W treatment/control assignment, coded as 0/1
#' @param target.pop which population should the treatment effect be estimated for?
#'            (0, 1): average treatment effect for everyone
#'            0: average treatment effect for controls
#'            1: average treatment effect for treated
#' @param fit.propensity should propensity model be used for variable selection?
#'
#' @return ATE estimate
#'
#' @export twostep.lasso.ate
twostep.lasso.ate = function(X, Y, W, target.pop=c(0, 1), fit.propensity = TRUE, estimate.se = FALSE) {

	# we want ATE for these indices
	target.idx = which(W %in% target.pop)
	
	lasso.fit0 = glmnet::cv.glmnet(X[W == 0,], Y[W == 0])
	lasso.fit1 = glmnet::cv.glmnet(X[W == 1,], Y[W == 1])
	strong.coef = which((coef(lasso.fit0)[-1] != 0) | (coef(lasso.fit1)[-1] != 0))
	
	if (fit.propensity) {
		lasso.prop = glmnet::cv.glmnet(X, W, family = "binomial")
		strong.coef = union(strong.coef, which(coef(lasso.prop)[-1] != 0))
	}
	
	if (length(strong.coef) == 0) {
		return(naive.ate(X, Y, W, estimate.se = TRUE))
	}
	
	twostep.raw = lapply(0:1, function(ww) {
		if (length(strong.coef) < 0.95 * sum(W == ww)) {
			coefs = strong.coef
		} else {
			coefs = which(strong.coef %in% coef(list(lasso.fit0, lasso.fit1)[[1 + ww]])[-1] != 0)
		}
	
		repeat {
			print(coefs)
			reg.df = data.frame(feat = X[,coefs], Y = Y, row.names = 1:nrow(X))
			center = apply(reg.df[target.idx,,drop=FALSE], 2, mean)
			center.df = data.frame(matrix(center, 1, ncol(reg.df)))
			names(center.df) = names(reg.df)
			fit.sel = lm(Y ~ ., data = reg.df[W == ww,,drop=FALSE])
			if(sum(is.na(coef(fit.sel))) == 0) {
				break
			} else {
				coefs = coefs[-sample.int(length(coefs), 1)]
			}
		}
		
		lm.pred = predict(fit.sel, newdata = center.df)
	
		if(estimate.se) {
			M = sandwich::sandwich(fit.sel)
			lm.var = c(1, center[-ncol(reg.df)]) %*% M %*%  c(1, center[-ncol(reg.df)]) 
		} else {
			lm.var = NA
		}
		
		c(lm.pred, lm.var)
	})
	
	tau.hat = twostep.raw[[2]][1] - twostep.raw[[1]][1]
	var.hat = twostep.raw[[2]][2] + twostep.raw[[1]][2]

	if(!estimate.se) {
		tau.hat
	} else {
		c(tau.hat, sqrt(var.hat))
	}
}
