library(xtable)
rm(list = ls())
source("process.fnm.R")

filenames = list.files("results", pattern="res-.-3-.*-1.RData", full.names=TRUE)

out.raw = process.fnm(filenames)

betas = c("dense", "", "harmonic", "", "moderatly sparse", "", "very sparse", "")
props = rep(c("0.1", "0.25"), 4)
results = data.frame(Beta=betas, Propensity=props, out.raw[,8:16])

nms = c("Beta Model", "Overlap", "Naive", "Elastic Net", "Approximate Balance", "Approx. Residual Balance", "Inverse Propensity Weight", "Augmented IPW", "Weighted Elastic Net", "TMLE Elastic Net", "Select + OLS")
resT = data.frame(t(results))
colnames(resT) = ""
rownames(resT) = nms

xtab = xtable(resT, align=c("|", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|"))
print(xtab, include.rownames = TRUE, include.colnames = FALSE, sanitize.text.function = identity, hline.after = c(0, 2, 4, 6, 8, 10, 11), file = "tables/manyclust.tex")

#xtab = xtable(results, align=c("r", "|", "r", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 1), "|"))
#print(xtab, include.rownames = FALSE, include.colnames = FALSE, sanitize.text.function = identity)
