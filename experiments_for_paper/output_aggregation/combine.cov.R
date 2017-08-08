library(xtable)
rm(list = ls())

filenames = list.files("results", pattern="res-.*-2.RData", full.names=TRUE)

ret = data.frame(t(sapply(filenames, function(fnm) {
		load(fnm)
		var.ratio = var(results[1,]) / mean(results[2,]^2)
		c(prop=prop.setup, beta=beta.setup, n=n, p=p, eps=eps, C=C, extra=extra.param, coverage=coverage, var.ratio= var.ratio)
})))

nvals = sort(as.numeric(as.character(unique(ret$n))))
pvals = sort(as.numeric(as.character(unique(ret$p))))

make.row = function(n, p) {
	ret.sub = ret[ret$n == n & ret$p == p,]
	ret.sub = ret.sub[order(as.numeric(as.character(ret.sub$C))),]
	ret.sub = ret.sub[order(as.numeric(as.character(ret.sub$eps)), decreasing = TRUE),]
	ret.sub = ret.sub[order(as.numeric(as.character(ret.sub$extra))),]
	ret.sub = ret.sub[order(as.numeric(as.character(ret.sub$beta))),]
	ret.sub = ret.sub[order(as.numeric(as.character(ret.sub$prop))),]
	c(n, p, sprintf("%.2f", round(ret.sub$coverage, 2)))
	#rbind(
	#	c(n, p, sprintf("%.2f", round(ret.sub$coverage, 2))),
	#	c("", "", sprintf("%.2f", round(ret.sub$var.ratio, 2))))
}

make.chunk = function(n) {
	Reduce(rbind, lapply(pvals, function(p) make.row(n, p)))
}

covs = Reduce(rbind, lapply(nvals, make.chunk))

#problem = c("", "", "many clust.", rep("", 3), "B/C/H", rep("", 3))
C.prop = c("", "", "$||\\beta||_0 = 10$", "", "$\\beta_j \\propto 1/j^2$", "", "$\\beta_j \\propto 1/j$", "")
C.main = c("$n$", "$p$", rep(c("$\\varepsilon = 0.25$", "$\\varepsilon = 0.1$"), 3))
results = data.frame(rbind(C.prop, C.main, covs))

xtab = xtable(results, align=c("r", "|", "r", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|"))
print(xtab, include.rownames = FALSE, include.colnames = FALSE, sanitize.text.function = identity, hline.after = c(0, 2, 5, 8, 11), file = "tables/coverage.tex")

#xtab = xtable(results, align=c("r", "|", "r", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 1), "|"))
#print(xtab, include.rownames = FALSE, include.colnames = FALSE, sanitize.text.function = identity)
