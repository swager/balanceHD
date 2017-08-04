n = 500
p = 700
M = matrix(rnorm(n*p), n, p)
balance.target = c(1, rep(0, p - 1))

gamma.mosek.positive = approx.balance(M, balance.target = balance.target, optimizer = "mosek", allow.negative.weights = FALSE, zeta = 0.5)
imbalance.mosek.positive = t(M) %*% gamma.mosek.positive - balance.target

gamma.mosek.free = approx.balance(M, balance.target = balance.target, optimizer = "mosek", allow.negative.weights = TRUE, zeta = 0.5)
imbalance.mosek.free = t(M) %*% gamma.mosek.free - balance.target

gamma.pogs.positive = approx.balance(M, balance.target = balance.target, optimizer = "pogs", allow.negative.weights = FALSE, zeta = 0.5)
imbalance.pogs.positive = t(M) %*% gamma.pogs.positive - balance.target

gamma.pogs.free = approx.balance(M, balance.target = balance.target, optimizer = "pogs", allow.negative.weights = TRUE, zeta = 0.5)
imbalance.pogs.free = t(M) %*% gamma.pogs.free - balance.target

gamma.qp.positive = approx.balance(M, balance.target = balance.target, optimizer = "quadprog", allow.negative.weights = FALSE, zeta = 0.5)
imbalance.qp.positive = t(M) %*% gamma.qp.positive - balance.target

gamma.qp.free = approx.balance(M, balance.target = balance.target, optimizer = "quadprog", allow.negative.weights = TRUE, zeta = 0.5)
imbalance.qp.free = t(M) %*% gamma.qp.free - balance.target

test_that("positivity constraint works", {
  expect_true(all(gamma.mosek.positive >0))
  expect_true(all(gamma.pogs.positive > -10^(-4)))
  expect_true(all(gamma.qp.positive > 0))
})

test_that("removing positivity constraint improves objective", {
  expect_true((sum(gamma.mosek.positive^2) + max(abs(imbalance.mosek.positive))^2) >
                (sum(gamma.mosek.free^2) + max(abs(imbalance.mosek.free))^2))
  expect_true((sum(gamma.pogs.positive^2) + max(abs(imbalance.pogs.positive))^2) >
                (sum(gamma.pogs.free^2) + max(abs(imbalance.pogs.free))^2))
  expect_true((sum(gamma.qp.positive^2) + max(abs(imbalance.qp.positive))^2) >
                (sum(gamma.qp.free^2) + max(abs(imbalance.qp.free))^2))
})

test_that("gamma sums to 1", {
  expect_equal(sum(gamma.mosek.positive), 1)
  expect_equal(sum(gamma.mosek.free), 1)
  expect_equal(sum(gamma.pogs.positive), 1, tolerance = 5e-3)
  expect_equal(sum(gamma.pogs.free), 1, tolerance = 5e-5)
  expect_equal(sum(gamma.qp.positive), 1)
  expect_equal(sum(gamma.qp.free), 1)
})

test_that("optimizers match", {
  expect_equal(max(abs(gamma.qp.free - gamma.mosek.free)), 0, tolerance = 5e-4)
  expect_equal(max(abs(gamma.qp.positive - gamma.mosek.positive)), 0, tolerance = 5e-4)
  expect_equal(max(abs(gamma.qp.free - gamma.pogs.free)), 0, tolerance = 5e-3)
  expect_equal(max(abs(gamma.qp.positive - gamma.pogs.positive)), 0, tolerance = 5e-2)
})

gamma.mosek.free.zeta9 = approx.balance(M, balance.target = balance.target, optimizer = "mosek", allow.negative.weights = TRUE, zeta = 0.9)
imbalance.mosek.free.zeta9 = t(M) %*% gamma.mosek.free.zeta9 - balance.target

gamma.pogs.free.zeta9 = approx.balance(M, balance.target = balance.target, optimizer = "pogs", allow.negative.weights = TRUE, zeta = 0.9)
imbalance.pogs.free.zeta9 = t(M) %*% gamma.pogs.free.zeta9 - balance.target

gamma.qp.free.zeta9 = approx.balance(M, balance.target = balance.target, optimizer = "quadprog", allow.negative.weights = TRUE, zeta = 0.9)
imbalance.qp.free.zeta9 = t(M) %*% gamma.qp.free.zeta9 - balance.target

test_that("zeta tunes problem in right direction", {
  expect_true(max(abs(imbalance.mosek.free.zeta9)) < max(abs(imbalance.mosek.free)))
  expect_true(sum(gamma.mosek.free.zeta9^2) > sum(gamma.mosek.free^2))
  expect_true(max(abs(imbalance.pogs.free.zeta9)) < max(abs(imbalance.pogs.free)))
  expect_true(sum(gamma.pogs.free.zeta9^2) > sum(gamma.pogs.free^2))
  expect_true(max(abs(imbalance.qp.free.zeta9)) < max(abs(imbalance.qp.free)))
  expect_true(sum(gamma.qp.free.zeta9^2) > sum(gamma.qp.free^2))
})