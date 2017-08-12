rm(list = ls())

filenames = list.files(path="results", pattern="output.*RData", full.names=TRUE)

load(filenames[1])

boot = 1000

treat.agg = colSums(site.raw * treat)
treat.f = treat.agg / sum(treat.agg)

tau.by.site = sapply(1:4, function(ss) mean(Y[treat == 1 & site == ss]) - mean(Y[treat == 0 & site == ss]))
tau.oracle.full = sum(treat.f * tau.by.site)

tau.oracle.boots = replicate(boot, {
	boot.idx = sample.int(length(treat), replace = TRUE)
	site.raw.boot = site.raw[boot.idx,]
	treat.boot = treat[boot.idx]
	Y.boot = Y[boot.idx]
	site.boot = site[boot.idx]
	
	treat.agg = colSums(site.raw.boot * treat.boot)
	treat.f = treat.agg / sum(treat.agg)
	tau.by.site = sapply(1:4, function(ss) mean(Y.boot[treat.boot == 1 & site.boot == ss]) - mean(Y.boot[treat.boot == 0 & site.boot == ss]))
	sum(treat.f * tau.by.site)
})


ret = lapply(filenames, function(fnm) {
	load(fnm)
	perf.boot = lapply(1:boot, function(bb) {
		err = results[c(1, 3, 5, 7, 9),] - tau.oracle.boots[bb]
		se = results[2 * (1:5),]
		perf = rbind(rowMeans(err^2), rowMeans(abs(err/se) > 1.98), rowMeans(abs(err/se) > 1.64))
	})
	perf.agg = Reduce(function(a,b) a + b, perf.boot) / length(perf.boot)
	
	list(num.samp, perf.agg)
})

save.image("summary.RData")
