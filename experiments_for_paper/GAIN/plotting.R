rm(list = ls())
load("summary.RData")

samps.raw = sapply(ret, function(elem) elem[[1]])
out.raw = lapply(ret, function(elem) elem[[2]])

samps = sort(samps.raw)
out = out.raw[order(samps.raw)]

pdf("plots/labor_mse.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
mse = Reduce(rbind, lapply(out, function(oo) oo[1,]))
cols = c(2, 4, 5, 1, 3)
plot(NA, NA, xlim = range(samps), ylim = pmin(range(mse, na.rm = TRUE), 0.12), xlab = "n", ylab = "MSE", log="xy")
for(iter in 1:5) {
	lines(samps, mse[,iter], col=cols[iter], lwd = 3)
	abline(h=0.95, lwd = 1, lty = 3)
}
legend("topright", c("Oracle Adjustment", "Approx. Resid. Balance", "Double Select + OLS", "Augmented IPW", "No Adjustment (Naive)"), lwd = 3, col = c(3, 4, 2, 5, 1), cex=1.5)
par = pardef
dev.off()

pdf("plots/labor_coverage.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
covg = Reduce(rbind, lapply(out, function(oo) oo[2,]))
plot(NA, NA, xlim = range(samps), ylim = range(1 - covg, na.rm = TRUE), xlab = "n", ylab = "Coverage", log="x")
cols = c(2, 4, 5, 1, 3)
for(iter in 1:5) {
	lines(samps, 1 - covg[,iter], col=cols[iter], lwd = 3)
}
abline(h=0.95, lwd = 1, lty = 3)
legend("bottomleft", c("Oracle Adjustment", "Approx. Resid. Balance", "Double Select + OLS", "Augmented IPW", "No Adjustment (Naive)"), lwd = 3, col = c(3, 4, 2, 5, 1), cex=1.5)
par = pardef
dev.off()

pdf("plots/labor_coverage_90.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
covg2 = Reduce(rbind, lapply(out, function(oo) oo[3,]))
plot(NA, NA, xlim = range(samps), ylim = range(1 - covg2, na.rm = TRUE), xlab = "n", ylab = "Coverage", log="x")
cols = c(2, 4, 5, 1, 3)
for(iter in 1:5) {
	lines(samps, 1 - covg2[,iter], col=cols[iter], lwd = 3)
}
abline(h=0.9, lwd = 1, lty = 3)
legend("bottomleft", c("Oracle Adjustment", "Approx. Resid. Balance", "Double Select + OLS", "Augmented IPW", "No Adjustment (Naive)"), lwd = 3, col = c(3, 4, 2, 5, 1), cex=1.5)
par = pardef
dev.off()

