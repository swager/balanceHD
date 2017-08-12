rm(list = ls())

library(balanceHD)
library(splines)

args=(commandArgs(TRUE))
num.samp = as.numeric(args[1])
NREP = as.numeric(args[2])

treat = as.numeric(read.csv("data/labor_treat.csv")$trainee)
site.raw = read.csv("data/labor_sample.csv")
site = factor(as.matrix(site.raw) %*% 1:4)
X = read.csv("data/labor_covariates.csv", header=FALSE)
Y.all = read.csv("data/labor_outcomes.csv", header=FALSE)
Y = rowMeans(Y.all[,1:36])

expand = function(vec) {
	mat = matrix(0, length(vec), 5)
	is.nonzero = which(vec > 0)
	mat[is.nonzero, 1] = 1
	mat[is.nonzero, 2:5] = ns(log(vec[is.nonzero]), 4)
	mat
}

# the original features # 11-30 were just boolean basis expansions
# for the first # 1-10
X.proc = Reduce(cbind, list(
	Reduce(cbind, lapply(1:10, function(iter) {
		expand(X[,iter])
	})),
	X[,31:34],
	Reduce(cbind, lapply(35:38, function(iter) {
		expand(X[,iter])
	})),
	X[,39:52],
	ns(X[,53], df = 5)
))

treat.agg = colSums(site.raw * treat)
treat.f = treat.agg / sum(treat.agg)

tau.by.site = sapply(1:4, function(ss) mean(Y[treat == 1 & site == ss]) - mean(Y[treat == 0 & site == ss]))
tau.star = sum(treat.f * tau.by.site)

results = replicate(NREP, {
	
idx = sample.int(length(treat), num.samp)

X.sub = as.matrix(X.proc[idx,])
Y.sub = Y[idx]
treat.sub = treat[idx]

tau.double = twostep.lasso.ate(X.sub, Y.sub, treat.sub, target.pop = 1, estimate.se=TRUE)
tau.arb = residualBalance.ate(X.sub, Y.sub, treat.sub, target.pop = 1, estimate.se = TRUE)
tau.ipw = ipw.ate(X.sub, Y.sub, treat.sub, target.pop = 1, estimate.se=TRUE)
tau.naive = naive.ate(X.sub, Y.sub, treat.sub, estimate.se = TRUE)

tau.site = sapply(1:4, function(ss) mean(Y.sub[treat.sub == 1 & site[idx] == ss]) - mean(Y.sub[treat.sub == 0 & site[idx] == ss]))
var.site = sapply(1:4, function(ss) var(Y.sub[treat.sub == 1 & site[idx] == ss]) / sum(treat.sub == 1 & site[idx] == ss) + var(Y.sub[treat.sub == 0 & site[idx] == ss])/ sum(treat.sub == 0 & site[idx] == ss))

treat.agg.sub = colSums(site.raw[idx,] * treat)
treat.f.sub = treat.agg.sub / sum(treat.agg.sub)

tau.oracle = c(sum(treat.f.sub * tau.site), sqrt(sum(treat.f.sub^2 * var.site)))

ret = c(DOUBLE=tau.double, ARB=tau.arb, IPW=tau.ipw, NAIVE=tau.naive, ORACLE=tau.oracle)

print(ret)
ret
})

err = results[c(1, 3, 5, 7, 9),] - tau.star
se = results[2 * (1:5),]

print(rowMeans(err^2))
print(rowMeans(abs(err/se) > 1.98))
print(rowMeans(abs(err/se) > 1.64))

save.image(paste0("results/output-", num.samp, ".RData"))
