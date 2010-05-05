library(mugnet)

data("lung")
data("lung.alive")
data("lung.dead")

priorOrder <- seq(1:nrow(lung.alive))

###################################################################
## BIC

eval.dead <- mgSearchOrder(lung.dead, NULL, 3, NULL, 1, 0, priorOrder, NULL, NULL, emIter=1000, stopDelta=0.01, selectMode="BIC", echo=TRUE)
eval.alive <- mgSearchOrder(lung.alive, NULL, 3, NULL, 1, 0, priorOrder, NULL, NULL, emIter=1000, stopDelta=0.01, selectMode="BIC", echo=TRUE)

eval <- eval.dead
for(i in 1:length(eval@nets)) {
  eval@nets[[i]]@likelihood <- sum(mgNodeLoglik(eval@nets[[i]], 1:eval@numnodes, lung.dead))
  eval@nets[[i]]@meta <- "PI3K pathway"
  eval@loglik[i] <- eval@nets[[i]]@likelihood
}
eval.dead <- eval
eval <- eval.alive
for(i in 1:length(eval@nets)) {
  eval@nets[[i]]@likelihood <- sum(mgNodeLoglik(eval@nets[[i]], 1:eval@numnodes, lung.alive))
  eval@nets[[i]]@meta <- "PI3K pathway"
  eval@loglik[i] <- eval@nets[[i]]@likelihood
}
eval.alive <- eval

bst0d <- cnFindBIC(eval.dead, ncol(lung.dead))
bst0d@meta <- "G0d, BIC"
bst0a <- cnFindBIC(eval.alive, ncol(lung.alive))
bst0a@meta <- "G0a, BIC"

###################################################################
## AIC

eval.dead <- mgSearchOrder(lung.dead, NULL, 3, NULL, 1, 0, priorOrder, NULL, NULL, emIter=1000, stopDelta=0.01, selectMode="AIC", echo=TRUE)
eval.alive <- mgSearchOrder(lung.alive, NULL, 3, NULL, 1, 0, priorOrder, NULL, NULL, emIter=1000, stopDelta=0.01, selectMode="AIC", echo=TRUE)
eval <- eval.dead
for(i in 1:length(eval@nets)) {
  eval@nets[[i]]@likelihood <- sum(mgNodeLoglik(eval@nets[[i]], 1:eval@numnodes, lung.dead))
  eval@nets[[i]]@meta <- "PI3K pathway"
  eval@loglik[i] <- eval@nets[[i]]@likelihood
}
eval.dead <- eval
eval <- eval.alive
for(i in 1:length(eval@nets)) {
  eval@nets[[i]]@likelihood <- sum(mgNodeLoglik(eval@nets[[i]], 1:eval@numnodes, lung.alive))
  eval@nets[[i]]@meta <- "PI3K pathway"
  eval@loglik[i] <- eval@nets[[i]]@likelihood
}
eval.alive <- eval

bst1d <- cnFindAIC(eval.dead)
bst1d@meta <- "G1d, AIC"
bst1a <- cnFindAIC(eval.alive)
bst1a@meta <- "G1a, AIC"

###################################################################
## 1-parent saturated

eval.dead <- mgSearchOrder(lung.dead, NULL, 3, NULL, 1, 0, priorOrder, NULL, NULL, emIter=1000, stopDelta=0.01, selectMode=0, echo=TRUE)
eval.alive <- mgSearchOrder(lung.alive, NULL, 3, NULL, 1, 0, priorOrder, NULL, NULL, emIter=1000, stopDelta=0.01, selectMode=0, echo=TRUE)
eval <- eval.dead
for(i in 1:length(eval@nets)) {
  eval@nets[[i]]@likelihood <- sum(mgNodeLoglik(eval@nets[[i]], 1:eval@numnodes, lung.dead))
  eval@nets[[i]]@meta <- "PI3K pathway"
  eval@loglik[i] <- eval@nets[[i]]@likelihood
}
eval.dead <- eval
eval <- eval.alive
for(i in 1:length(eval@nets)) {
  eval@nets[[i]]@likelihood <- sum(mgNodeLoglik(eval@nets[[i]], 1:eval@numnodes, lung.alive))
  eval@nets[[i]]@meta <- "PI3K pathway"
  eval@loglik[i] <- eval@nets[[i]]@likelihood
}
eval.alive <- eval

bst2d <- cnFind(eval.dead, 86)
bst2d@meta <- "G2d"
bst2a <- cnFind(eval.alive, 86)
bst2a@meta <- "G2a"

###################################################################
## Goodness-of-Fit
## Apply Kolmogorov-Smirnoff test

ps <- mgSamples(bst0a, 10000)
kspval <- rep(0,bst0a@numnodes)
for(i in 1:bst1a@numnodes) {
  res <- ks.test(lung[i,], ps[,i])
  kspval[i] <- res$p.value
}
names(kspval) <- bst2d@nodes
kspval

ps <- mgSamples(bst2a, 100000)
kspval <- rep(0,bst2a@numnodes)
for(i in 1:bst2a@numnodes) {
  res <- ks.test(lung[i,], ps[,i])
  kspval[i] <- res$p.value
}
names(kspval) <- bst2d@nodes
kspval

###################################################################
## probability plots

postscript("pmf_igf1.ps")
par(mfrow=c(1,2))
barplot(bst2a@probabilities[[2]][[1]], main="G2a(class alive)", xlab="categories", ylab="p(IGF1)")
barplot(bst2d@probabilities[[2]][[1]], main="G2d(class diseased)", xlab="categories")
dev.off()

## bst2d
postscript("pmf_bst2d.ps")
par(mfrow=c(5,3))
barplot(bst2d@probabilities[[2]][[1]], main="IGF1=1", ylab="p(PIK3CA|IGF1)")
barplot(bst2d@probabilities[[2]][[2]], main="IGF1=2")
barplot(bst2d@probabilities[[2]][[3]], main="IGF1=3")
barplot(bst2d@probabilities[[3]][[1]], main="PIK3CA=1", ylab="p(PTEN|PIK3CA)")
barplot(bst2d@probabilities[[3]][[2]], main="PIK3CA=2")
barplot(bst2d@probabilities[[3]][[3]], main="PIK3CA=3")
barplot(bst2d@probabilities[[4]][[1]], main="PIK3CA=1", ylab="p(PDK1|PIK3CA)")
barplot(bst2d@probabilities[[4]][[2]], main="PIK3CA=2")
barplot(bst2d@probabilities[[4]][[3]], main="PIK3CA=3")
barplot(bst2d@probabilities[[5]][[1]], main="PIK3CA=1", ylab="p(AKT1|PIK3CA)")
barplot(bst2d@probabilities[[5]][[2]], main="PIK3CA=2")
barplot(bst2d@probabilities[[5]][[3]], main="PIK3CA=3")
barplot(bst2d@probabilities[[8]][[1]], main="AKT1=1", ylab="p(GSK3A|AKT1)")
barplot(bst2d@probabilities[[8]][[2]], main="AKT1=2")
barplot(bst2d@probabilities[[8]][[3]], main="AKT1=3")
dev.off()

## bst2a
postscript("pmf_bst2a.ps")
par(mfrow=c(5,3))
barplot(bst2a@probabilities[[2]][[1]], main="IGF1=1", ylab="p(PIK3CA|IGF1)")
barplot(bst2a@probabilities[[2]][[2]], main="IGF1=2")
barplot(bst2a@probabilities[[2]][[3]], main="IGF1=3")
barplot(bst2a@probabilities[[3]][[1]], main="PIK3CA=1", ylab="p(PTEN|PIK3CA)")
barplot(bst2a@probabilities[[3]][[2]], main="PIK3CA=2")
barplot(bst2a@probabilities[[3]][[3]], main="PIK3CA=3")
barplot(bst2a@probabilities[[4]][[1]], main="PTEN=1", ylab="p(PDK1|PTEN)")
barplot(bst2a@probabilities[[4]][[2]], main="PTEN=2")
barplot(bst2a@probabilities[[4]][[3]], main="PTEN=3")
barplot(bst2a@probabilities[[5]][[1]], main="IGF1", ylab="p(AKT1|IGF1)")
barplot(bst2a@probabilities[[5]][[2]], main="IGF1=2")
barplot(bst2a@probabilities[[5]][[3]], main="IGF1=3")
barplot(bst2a@probabilities[[8]][[1]], main="AKT1=1", ylab="p(GSK3A|AKT1)")
barplot(bst2a@probabilities[[8]][[2]], main="AKT1=2")
barplot(bst2a@probabilities[[8]][[3]], main="AKT1=3")
dev.off()


par(mfrow=c(2,1))
ps <- mgSamples(bst2d, 10000)
hist(ps[,"IGF1"], 30, main="", xlab="G2d[IGF1]")
hist(lung.dead["IGF1",],30, main="", xlab="data[IGF1]")

cnDot(list(bst0a),"bst0a")
cnDot(list(bst0d),"bst0d")
cnDot(list(bst1a),"bst1a")
cnDot(list(bst1d),"bst1d")
cnDot(list(bst2a),"bst2a")
cnDot(list(bst2d),"bst2d")
