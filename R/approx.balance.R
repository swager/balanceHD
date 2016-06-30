#' Compute approximately balancing weights
#'
#' Returns the minimizer of:
#'   (1 - zeta) ||gamma||^2 + zeta ||M'gamma - balance.target||_infty^2 (*)
#'
#' @param M the feature matrix, see (*)
#' @param balance.target the target solution, see (*)
#' @param zeta tuning parameter, see (*)
#' @param allow.negative.weights are the gammas allowed to be negative?
#'
#' @return gamma, the minimizer of (*)
#'
#' @export approx.balance
approx.balance = function(M, balance.target, zeta = 0.5, allow.negative.weights = FALSE) {
	
	if (zeta <= 0 || zeta >= 1) {
		stop("approx.balance: zeta must be between 0 and 1")
	}
	
	# The system is effectively
	# minimize zeta * delta + (1 - zeta) * ||gamma||^2
	# subject to
	#   sum gamma = 1
	#   delta + (M'gamma)_j >= balance.target_j
	#   delta - (M'gamma)_j >= -balance.target_j
	#
	# The last two constraints mean that
	# delta = ||M'gamma - balance.target||_infty
	
	Dmat = diag(c(zeta, rep(1 - zeta, nrow(M))))
	dvec = rep(0, 1 + nrow(M))
	Amat = cbind(
		c(0, rep(1, nrow(M))),
		rbind(rep(1, ncol(M)), M ),
		rbind(rep(1, ncol(M)), -M))
	bvec = c(1, balance.target, -balance.target)
	
	if (!allow.negative.weights) {
		LB = 1/nrow(M)/10000
		Amat = cbind(Amat, rbind(rep(0, nrow(M)), diag(rep(1, nrow(M)))))
		bvec = c(bvec, rep(LB, nrow(M)))
	}
	
	balance.soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
	gamma = balance.soln$solution[-1]
	gamma
}
