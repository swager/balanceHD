clean = function(xx) {
  xxmin = min(xx)
  idx = which(xx <= 1.05 * xxmin)
  xtx = sapply(xx, function(xxx)sprintf("%.3f", round(xxx, 3)))
  xtx[idx] = paste(" \\bf", xtx[idx])
  xtx
}

process.fnm = function(filenames) {
	ret = sapply(filenames, function(fnm) {
		load(fnm)
		spec = c(beta=beta.setup, prop=prop.setup, n=n, p=p, eps=eps, C=C, extra=extra.param)
		
		tau = results[1,]
		estimates = results[-1,]
		err = t(t(estimates) - tau)
		rmse.est = sqrt(rowMeans(err^2))
		
		err.summary = clean(rmse.est)
		err.summary = err.summary[c(1, 2, 5, 4, 6, 7, 8, 9, 3)]
		c(spec, err.summary)
	})
	data.frame(t(ret))
}
