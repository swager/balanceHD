n = 2000
p = 20
tau = 10
beta = 2 / (1:p) / sqrt(sum(1/(1:p)^2))

# Generate data
X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, 1/(1 + exp(-(X[,1] + X[,2])/4)))
Y = X %*% beta + rnorm(n, 0, 1) + tau * W

tau.aipw = ipw.ate(X, Y, W, fit.method = "elnet", estimate.se = TRUE)
tau.tmle = ipw.ate(X, Y, W, fit.method = "elnet", targeting.method = "TMLE", estimate.se = TRUE)
tau.weighted = ipw.ate(X, Y, W, fit.method = "elnet", prop.weighted.fit = TRUE, estimate.se = TRUE)

test_that("AIPW, TMLE and weighting match on easy problem", {
  expect_true(abs(tau.aipw[1] - tau.tmle[1]) <= 0.015)
  expect_true(abs(tau.aipw[2] - tau.tmle[2]) <= 0.005)
  expect_true(abs(tau.aipw[1] - tau.weighted[1]) <= 0.015)
  expect_true(abs(tau.aipw[2] - tau.weighted[2]) <= 0.005)
})

att.aipw = ipw.ate(X, Y, W, target.pop = 1, fit.method = "elnet", estimate.se = TRUE)
att.tmle = ipw.ate(X, Y, W, target.pop = 1, fit.method = "elnet", targeting.method = "TMLE", estimate.se = TRUE)
att.weighted = ipw.ate(X, Y, W, target.pop = 1, fit.method = "elnet", prop.weighted.fit = TRUE, estimate.se = TRUE)

test_that("AIPW, TMLE and weighting match on easy problem", {
  expect_true(abs(att.aipw[1] - att.tmle[1]) <= 0.015)
  expect_true(abs(att.aipw[2] - att.tmle[2]) <= 0.005)
  expect_true(abs(att.aipw[1] - att.weighted[1]) <= 0.015)
  expect_true(abs(att.aipw[2] - att.weighted[2]) <= 0.005)
})

n = 400
p = 1000
delta.clust = 40 / sqrt(n) * rep(c(1, rep(0, 9)), p/10)
beta.raw = 1/(9 + 1:p)
beta.main = 2 * beta.raw / sqrt(sum(beta.raw^2))

CLUST = rbinom(n, 1, 0.5)
W = rbinom(n, 1, 0.2 + 0.6 * CLUST)
X = matrix(rnorm(n * p), n, p) + outer(CLUST, delta.clust)
Y = X %*% beta.main + rnorm(n)

tau.elnet.cate = elnet.ate(X, Y, W, target.pop = c(0, 1), alpha = 0.9)
tau.aipw.cate = ipw.ate(X, Y, W, target.pop = c(0, 1), prop.weighted.fit = FALSE, targeting.method = "AIPW",
                        fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)

tau.elnet.catt = elnet.ate(X, Y, W, target.pop = 1, alpha = 0.9)
tau.aipw.catt = ipw.ate(X, Y, W, target.pop = 1, prop.weighted.fit = FALSE, targeting.method = "AIPW",
                   fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)

test_that("AIPW improves over IPW.", {
  expect_true(abs(tau.aipw.cate) < abs(tau.elnet.cate))
  expect_true(abs(tau.aipw.catt) < abs(tau.elnet.catt))
})

