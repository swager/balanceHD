#' No correction for confounding variables (difference in means)
#'
#' @param X the input features
#' @param Y the observed responses
#' @param W treatment/control assignment, coded as 0/1
#'
#' @return ATE estimate
#'
#' @export naive.ate
naive.ate = function(X, Y, W) {
	tau.hat = mean(Y[W==1]) - mean(Y[W==0])
	return(tau.hat)
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
elnet.ate = function(X, Y, W, target.pop=c(0, 1), alpha = 0.5) {

	# we want ATE for these indices
	target.idx = which(W %in% target.pop)
	
	fit0 = glmnet::cv.glmnet(X[W == 0,], Y[W == 0], alpha = alpha)
	fit1 = glmnet::cv.glmnet(X[W == 1,], Y[W == 1], alpha = alpha)
	
	pred0 = coef(fit0)[1] + X %*% coef(fit0)[-1]
	pred1 = coef(fit1)[1] + X %*% coef(fit1)[-1]
	
	mean((pred1 - pred0)[target.idx])
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
twostep.lasso.ate = function(X, Y, W, target.pop=c(0, 1), fit.propensity = TRUE) {

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
		
		return(naive.ate(X, Y, W, estimate.se = estimate.se))
		
	}
	
	reg.df = data.frame(feat = X[, strong.coef], Y = Y, row.names = 1:nrow(X))
	
	center = matrix(colMeans(reg.df[target.idx,]), 1, ncol(reg.df))
	colnames(center) = names(reg.df)
	center.df = data.frame(center)
	
	
	if (length(strong.coef) < 0.95 * sum(W == 0)) {
		fit0 = lm(Y ~ ., data = reg.df[W == 0,])
		lm.pred0 = predict(fit0, newdata = center.df, se.fit=TRUE)
		pred0 = c(fit=as.numeric(lm.pred0[[1]][1]), se.fit=lm.pred0[[2]][1])
	} else {
		lasso.pred0 = mean(predict(lasso.fit0, newx = X[target.idx,]))
		pred0 = c(fit=lasso.pred0, se.fit=NA)
	}
	
	if (length(strong.coef) < 0.95 * sum(W == 1)) {
		fit1 = lm(Y ~ ., data = reg.df[W == 1,])
		lm.pred1 = predict(fit1, newdata = center.df, se.fit=TRUE)
		pred1 = c(fit=as.numeric(lm.pred1[[1]][1]), se.fit=lm.pred1[[2]][1])
	} else {
		lasso.pred1 = mean(predict(lasso.fit1, newx = X[target.idx,]))
		pred1 = c(fit=lasso.pred1, se.fit=NA)
	}
	
	tau.hat = as.numeric(pred1["fit"] - pred0["fit"])
	tau.hat
}
