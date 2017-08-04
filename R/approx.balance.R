#' Compute approximately balancing weights
#'
#' Returns the minimizer of:
#'   (1 - zeta) ||gamma||^2 + zeta ||M'gamma - balance.target||_infty^2 (*)
#'
#' @param M the feature matrix, see (*)
#' @param balance.target the target solution, see (*)
#' @param zeta tuning parameter, see (*)
#' @param allow.negative.weights are the gammas allowed to be negative?
#' @param optimizer Which optimizer to use? Mosek is a commercial solver, but free
#'                  academic licenses are available. Needs to be installed separately.
#'                  Pogs runs ADMM and may be useful for large problems, and
#'                  must be installed separately. Quadprog is the default
#'                  R solver.
#' @param verbose whether the optimizer should print progress information
#'
#' @return gamma, the minimizer of (*)
#'
#' @export approx.balance
approx.balance = function(M,
                          balance.target,
                          zeta = 0.5,
                          allow.negative.weights = FALSE,
                          optimizer = c("mosek", "pogs", "quadprog"),
                          verbose=FALSE) {
	
	if (zeta <= 0 || zeta >= 1) {
		stop("approx.balance: zeta must be between 0 and 1")
	}
  
  optimizer = match.arg(optimizer)
  
  if (optimizer == "mosek") {
    if (suppressWarnings(require("Rmosek", quietly = TRUE))) {
      gamma = approx.balance.mosek(M, balance.target, zeta, allow.negative.weights, verbose)
    } else {
      warning("The mosek optimizer is not installed. Using quadprog instead.")
      optimizer = "quadprog"
    }
  } else if (optimizer == "pogs") {
    if (suppressWarnings(require("pogs", quietly = TRUE))) {
      gamma = approx.balance.pogs(M, balance.target, zeta, allow.negative.weights)
    } else {
      warning("The POGS optimizer is not installed. Using quadprog instead.")
      optimizer = "quadprog"
    }
  }
  
  if (optimizer == "quadprog") {
    gamma = approx.balance.quadprog(M, balance.target, zeta, allow.negative.weights)
  }
  
  gamma
}

# Find approximately balancing weights using quadprog
approx.balance.quadprog = function(M,
                                   balance.target,
                                   zeta = 0.5,
                                   allow.negative.weights = FALSE) {
        # The system is effectively
        # minimize zeta * delta^2 + (1 - zeta) * ||gamma||^2
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

# Find approximately balancing weights using mosek
approx.balance.mosek = function(M,
                                balance.target,
                                zeta = 0.5,
                                allow.negative.weights = FALSE,
                                verbose = FALSE) {
	# The system is effectively
	# minimize zeta * delta^2 + (1 - zeta) * ||gamma||^2
	# subject to
	#   sum gamma = 1
	#    -Infty <= -delta + (M'gamma)_j <= balance.target_j
	#   balance.target_j <= delta + (M'gamma)_j <= Infty
	#   0 <= gamma <= gamma.max (lower bound is optional)
	#
	# The second and third constraints mean that
	# delta = ||M'gamma - balance.target||_infty
	Qmat = list(i = 1:(1 + nrow(M)),
	            j = 1:(1 + nrow(M)),
	            v = 2 * c(zeta, rep(1 - zeta, nrow(M))))
	cvec = rep(0, 1 + nrow(M))
	Amat = Matrix::Matrix(rbind(
	  	c(0, rep(1, nrow(M))),
		cbind(rep(-1, ncol(M)), t(M) ),
		cbind(rep(+1, ncol(M)), t(M))))
	Amat = as(Amat, "CsparseMatrix")
	buc = c(1, balance.target,  rep(Inf, ncol(M)))
	blc = c(1, rep(-Inf, ncol(M)), balance.target)
	
	gamma.max = 1/nrow(M)^(2/3)
	bux = c(Inf, rep(gamma.max, nrow(M)))
	
	if (allow.negative.weights) {
		blx = c(0, -rep(gamma.max, nrow(M)))
	} else {
	 	blx = rep(0, nrow(M) + 1)
	}
	
	mosek.problem <- list()
	mosek.problem$sense <- "min"
	mosek.problem$qobj <- Qmat
	mosek.problem$c <- cvec
	mosek.problem$A <-Amat
	mosek.problem$bc <- rbind(blc = blc, buc = buc)
	mosek.problem$bx <- rbind(blx = blx, bux = bux)
	
	if (verbose) {
		mosek.out = Rmosek::mosek(mosek.problem)
	} else {
		mosek.out = Rmosek::mosek(mosek.problem, opts=list(verbose=0))
	}
	
	delta = mosek.out$sol$itr$xx[1]
	gamma =  mosek.out$sol$itr$xx[1 + 1:nrow(M)]
	gamma
}


# Find approximately balancing weights using pogs
approx.balance.pogs = function(M,
                               balance.target,
                               zeta = 0.5,
                               allow.negative.weights = FALSE) {
  
  # Our original problem is the following:
  #
  # Minimize ||gamma||_2^2 + delta^2 / lambda^2, subject to
  # ||M'gamma - balance.target||_infty <= delta and sum gamma_i = 0.
  #
  # Equivalently, in notation recognized by POGS, we can write
  #
  # Minimize ||gamma||_2^2 + Z^2
  # subject to I[1:2p] <= 0, J[2p+1, 2p+2] = 1, where
  #
  #     (M'  -v -lambda)   (gamma)
  # I = (-M' v  -lambda) * (ONE  )
  #     (1'  0  0      )   (Z    )
  #     (0   1  0      )
  #
  # and v denotes balance.target.
	
	lambda = sqrt((1 - zeta) / zeta)
	g = list(h = kSquare())
	
	if (allow.negative.weights) {
	  f = list(h = c(kIndLe0(2 * ncol(M)), kIndEq0(2)),
	           b = c(rep(0, 2 * ncol(M)), 1, 1))
	  A = rbind(cbind(t(M), -balance.target, -lambda),
	            cbind(-t(M), balance.target, -lambda),
	            c(rep(1, nrow(M)), 0, 0),
	            c(rep(0, nrow(M)), 1, 0))
	} else {
	  f = list(h = c(kIndLe0(2 * ncol(M) + nrow(M)), kIndEq0(2)),
	           b = c(rep(0, 2 * ncol(M) + nrow(M)), 1, 1))
	  A = rbind(cbind(t(M), -balance.target, -lambda),
	            cbind(-t(M), balance.target, -lambda),
	            cbind(diag(-1, nrow(M)), 0, 0),
	            c(rep(1, nrow(M)), 0, 0),
	            c(rep(0, nrow(M)), 1, 0))
	}

  
  pogs.solution = pogs(A, f, g, params = list(rel_tol=1e-4, abs_tol=1e-5))
  
  gamma = pogs.solution$x[1:nrow(M)]
  gamma
}
