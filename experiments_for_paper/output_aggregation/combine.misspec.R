library(xtable)
rm(list = ls())
source("process.fnm.R")

filenames = list.files("results", pattern="res-.-4-.*-1.RData", full.names=TRUE)

out.raw = process.fnm(filenames)
out.raw$n = as.numeric(as.character(out.raw$n))
out.raw$p = as.numeric(as.character(out.raw$p))
out.raw = out.raw[order(out.raw$p),]
out.raw = out.raw[order(out.raw$n),]

nn = c("400", "", "", "", "", "1000", "", "", "", "")
pp = rep(c("100", "200", "400", "800", "1600"), 2)
results = data.frame(n=nn, p=pp, out.raw[,8:16])

nms = c("$n$", "$p$", "Naive", "Elastic Net", "Approximate Balance", "Approx. Residual Balance", "Inverse Propensity Weight", "Augmented IPW", "Weighted Elastic Net", "TMLE Elastic Net", "Select + OLS")
resT = data.frame(t(results))
colnames(resT) = ""
rownames(resT) = nms

xtab = xtable(resT, align=c("|", "r", "|", rep("c", 5), "|", rep("c", 5), "|"))
print(xtab, include.rownames = TRUE, include.colnames = FALSE, sanitize.text.function = identity, hline.after = c(0, 2, 4, 6, 8, 10, 11), file = "tables/misspec.tex")

#xtab = xtable(results, align=c("r", "|", "r", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 1), "|"))
#print(xtab, include.rownames = FALSE, include.colnames = FALSE, sanitize.text.function = identity)
