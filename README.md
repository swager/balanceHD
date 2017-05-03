# balanceHD
Estimation of average treatment effects in high dimensions via approximate residual balancing, as proposed by Athey et al. (2016).

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("swager/balanceHD")
```

The package currently works with two optimizers: `pogs` and `quadprog`. `quadprog` is a standard `R` optimization library, and is available directly from CRAN. `pogs` is a more performant optimizer, based on ADMM, which must be installed separately; see [this page](https://github.com/foges/pogs/blob/master/src/interface_r/README.md) for instructions. For large problems (say, n > 1000) we recommend using `pogs`.

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
<b>Efficient Inference of Average Treatment Effects in High Dimensions via Approximate Residual Balancing.</b>
2016.
[<a href="http://arxiv.org/pdf/1604.07125.pdf">arxiv</a>]
