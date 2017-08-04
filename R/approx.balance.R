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
#' @param use.dual whether balancing should be solved in dual form
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
                          use.dual = NULL,
                          verbose = FALSE) {
	
	if (zeta <= 0 || zeta >= 1) {
		stop("approx.balance: zeta must be between 0 and 1")
	}
  
  optimizer = match.arg(optimizer)
  if(is.null(use.dual)) {
  	use.dual = (optimizer %in% c("pogs", "mosek"))
  }
  
  if (optimizer == "mosek") {
    if (suppressWarnings(require("Rmosek", quietly = TRUE))) {
    	if (use.dual) {
      		gamma = approx.balance.mosek.dual(M, balance.target, zeta, allow.negative.weights, verbose)
      	} else {
      		gamma = approx.balance.mosek(M, balance.target, zeta, allow.negative.weights, verbose)
      	}
    } else {
      warning("The mosek optimizer is not installed. Using quadprog instead.")
      optimizer = "quadprog"
    }
  } else if (optimizer == "pogs") {
    if (suppressWarnings(require("pogs", quietly = TRUE))) {
      if (use.dual) {
      	gamma = approx.balance.pogs.dual(M, balance.target, zeta, allow.negative.weights, verbose)
      } else {
      	gamma = approx.balance.pogs(M, balance.target, zeta, allow.negative.weights, verbose)
      }
    } else {
      warning("The POGS optimizer is not installed. Using quadprog instead.")
      optimizer = "quadprog"
    }
  }
  
  if (optimizer == "quadprog") {
  	if (use.dual) {warning("Dual solution not yet implemented for quadprog; using primal.")}
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

# Find approximately balancing weights using mosek, using dual of QP
approx.balance.mosek.dual = function(M,
                                balance.target,
                                zeta = 0.5,
                                allow.negative.weights = FALSE,
                                verbose = FALSE) {

	# The primal problem is:
	# Minimize 1/2 x' diag(qvec) x
	# subject to Ax <= b,
	# where the constraints indexted by "equality" are required to be equalities
	# 
	# Here, the vector x is interpreted as (max.imbalance, gamma)
	qvec.primal = 2 * c(zeta, rep(1 - zeta, nrow(M)))
	A.primal = rbind(
	  	c(0, rep(1, nrow(M))),
		cbind(rep(-1, ncol(M)), t(M) ),
		cbind(rep(-1, ncol(M)), -t(M)))
	b.primal = c(1, balance.target, -balance.target)
	equality.primal = c(TRUE, rep(FALSE, 2*ncol(M)))

	if (!allow.negative.weights) {
		A.primal = rbind(A.primal, cbind(0, diag(-1, nrow(M), nrow(M))))
		b.primal = c(b.primal, rep(0, nrow(M)))
		equality.primal = c(equality.primal, rep(FALSE, nrow(M)))
	}

	# The dual problem is then
	# Minimize 1/2 lambda A diag(1/qvec) A' lambda + b * lambda
	# Subject to lambda >= 0 for the lambdas corresponding to inequality constraints
	#
	# For computational purposes, we introduce additional variables mu = diag(1/sqrt(qvec)) A' lambda;
	# this adds extra linear constraints, but the objective is now 1/2 ||mu||_2^2.
	#
	# The solution to the optimization problem is (lambda, mu)
	nvar = length(qvec.primal)
	A.dual = Matrix::Matrix(cbind(diag(1/sqrt(qvec.primal)) %*% t(A.primal), diag(-1, nvar, nvar)))
	A.dual = as(A.dual, "CsparseMatrix")
	blc.dual = rep(0, nvar)
	buc.dual = rep(0, nvar)

	blx.dual = c(rep(0, length(equality.primal)), rep(-Inf, nvar))
	blx.dual[which(equality.primal)] = -Inf
	
	Q.dual = list(i=nrow(A.primal) + 1:nvar, j=nrow(A.primal) + 1:nvar, v=rep(1, nvar))
	c.dual = c(b.primal, rep(0, nvar))

	mosek.problem <- list()
	mosek.problem$sense <- "min"
	mosek.problem$qobj <- Q.dual
	mosek.problem$c <- c.dual
	mosek.problem$bx <- rbind(blx = blx.dual, bux = rep(Inf, length(blx.dual)))
	mosek.problem$A <- A.dual
	mosek.problem$bc <- rbind(blc = blc.dual, buc = buc.dual)
	
	if (verbose) {
		mosek.out = Rmosek::mosek(mosek.problem)
	} else {
		mosek.out = Rmosek::mosek(mosek.problem, opts=list(verbose=0))
	}
	
	primal = -diag(1/qvec.primal) %*% t(A.primal) %*% mosek.out$sol$itr$xx[1:nrow(A.primal)]
	delta = primal[1]
	gamma = primal[1 + 1:nrow(M)]
	gamma / sum(gamma)
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
                               allow.negative.weights = FALSE,
                               verbose = FALSE) {
  
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

  
  pogs.solution = pogs(A, f, g, params = list(rel_tol=1e-4, abs_tol=1e-5, verbose=2*as.numeric(verbose)))
  
  gamma = pogs.solution$x[1:nrow(M)]
  gamma
}

# Find approximately balancing weights using pogs
approx.balance.pogs.dual = function(M,
                               balance.target,
                               zeta = 0.5,
                               allow.negative.weights = FALSE,
                               verbose = FALSE) {
  
 	# The primal problem is:
	# Minimize 1/2 x' diag(qvec) x
	# subject to Ax <= b,
	# where the constraints indexted by "equality" are required to be equalities
	# 
	# Here, the vector x is interpreted as (max.imbalance, gamma)
	qvec.primal = 2 * c(zeta, rep(1 - zeta, nrow(M)))
	A.primal = rbind(
	  	c(0, rep(1, nrow(M))),
		cbind(rep(-1, ncol(M)), t(M) ),
		cbind(rep(-1, ncol(M)), -t(M)))
	b.primal = c(1, balance.target, -balance.target)
	equality.primal = c(TRUE, rep(FALSE, 2*ncol(M)))

	if (!allow.negative.weights) {
		A.primal = rbind(A.primal, cbind(0, diag(-1, nrow(M), nrow(M))))
		b.primal = c(b.primal, rep(0, nrow(M)))
		equality.primal = c(equality.primal, rep(FALSE, nrow(M)))
	}

	# The dual problem is then
	# Minimize 1/2 lambda A diag(1/qvec) A' lambda + b * lambda
	# Subject to lambda >= 0 for the lambdas corresponding to inequality constraint
	#
	# This is equivalent to
	# Minimize 1/2 ||c||_2^2 + d
	# Subject to c = diag(1/sqrt(qvec)) A' lambda
	#  d = b * lambda
	#  lambda >= 0

	A.pogs = rbind(diag(1/sqrt(qvec.primal)) %*% t(A.primal), b.primal)
	f.pogs = list(h = c(kSquare(ncol(A.primal)), kIdentity(1)))
	g.pogs = list(h = kIndGe0(length(equality.primal)))
	g.pogs$h[equality.primal] = kZero(sum(equality.primal))

  
  	pogs.solution = pogs(A.pogs, f.pogs, g.pogs, params = list(rel_tol=1e-6, abs_tol=1e-7, verbose=2*as.numeric(verbose)))
  
  	primal = -diag(1/sqrt(qvec.primal)) %*% t(A.primal) %*% pogs.solution$x
	delta = primal[1]
	gamma = primal[1 + 1:nrow(M)]

	gamma / sum(gamma)
}
