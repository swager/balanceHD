# balanceHD
Estimation of average treatment effects in high dimensions via approximate residual balancing, as proposed by Athey et al. (2016).

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("swager/balanceHD")
```

This package currently works with three optimizers: `mosek`, `pogs`, and `quadprog`.
Mosek is a commercial interior point solver,
pogs is a first-order optimizer, based on ADMM, while
quadprog is a standard `R` optimization library.
In general, we achieved best performance with mosek, and
recommend trying optimizers in the order listed above.
We found pogs to be somewhat slower than mosek on the problems
we tried. (Note that we offer two solution strategies based on pogs: `pogs` and
`pogs.dual`. We usually recommend the former, except when p is much larger than n.)
Finally, quadprog performors well on small problems, but can be much
slower for larger problems.

In terms of availability, the optimizer quadprog is easiest to access,
and is available directly from CRAN. Pogs needs to be installed separately,
but is still free. To install pogs, simply install one of the pre-compiled
binaries available from the project [repository](https://github.com/foges/pogs);
see [this page](https://github.com/foges/pogs/blob/master/src/interface_r/README.md)
for furhter instructions.
Finally, mosek is a commercial solver; however, academic
[licenses](https://www.mosek.com/resources/academic-license) are available for free.
One mosek has been installed, we call into it using the R package Rmosek.

Example usage:

```R
library(balanceHD)

n = 400
p = 1000
tau = 7
nclust = 10
beta = 2 / (1:p) / sqrt(sum(1/(1:p)^2))
clust.ptreat = rep(c(0.1, 0.9), nclust/2)

cluster.center = 0.5 * matrix(rnorm(nclust * p), nclust, p)
cluster = sample.int(nclust, n, replace = TRUE)
X = cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, clust.ptreat[cluster])
Y = X %*% beta + rnorm(n, 0, 1) + tau * W

tau.hat = residualBalance.ate(X, Y, W, estimate.se = TRUE)
print(paste("true tau:", tau))
print(paste("point estimate:", round(tau.hat[1], 2)))
print(paste0("95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2), ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"))
```

#### References
Susan Athey, Guido Imbens, and Stefan Wager.
<b>Approximate Residual Balancing: De-Biased Inference of Average Treatment Effects in High Dimensions.</b>
2016.
[<a href="http://arxiv.org/pdf/1604.07125.pdf">arxiv</a>]
