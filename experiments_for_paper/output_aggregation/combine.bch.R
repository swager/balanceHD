library(xtable)
rm(list = ls())
source("process.fnm.R")

filenames = list.files("results", pattern="res-1-6-100-200-.*-1.RData", full.names=TRUE)

out.raw = process.fnm(filenames)
out.raw = out.raw[order(as.numeric(as.character(out.raw$C))),]
out.raw = out.raw[order(as.numeric(as.character(out.raw$eps))),]
out.raw = out.raw[order(as.numeric(as.character(out.raw$extra))),]

out.raw = out.raw[which(out.raw$eps %in% c(1, 4)),]

props = c("sparse", rep("", 3), "dense", rep("", 3))
C.prop = c("1", "", "4", "", "1", "", "4", "")
C.main = rep(c("1", "4"), 4)
results = data.frame(Prop=props, C1=C.prop, C2=C.main, out.raw[,8:16])

nms = c("Propensity Model", "First Stage Sig. Strength", "Structure Sig. Strength", "Naive", "Elastic Net", "Approximate Balance", "Approx. Residual Balance", "Inverse Propensity Weight", "Augmented IPW", "Weighted Elastic Net", "TMLE Elastic Net", "Select + OLS")
resT = data.frame(t(results))
colnames(resT) = ""
rownames(resT) = nms

xtab = xtable(resT, align=c("|", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|"))
print(xtab, include.rownames = TRUE, include.colnames = FALSE, sanitize.text.function = identity, hline.after = c(0, 3, 5, 7, 9, 11, 12), file = "tables/bch.tex")

#xtab = xtable(results, align=c("r", "|", "r", "r", "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 1), "|"))
#print(xtab, include.rownames = FALSE, include.colnames = FALSE, sanitize.text.function = identity)
