library(mugnet)

data(crabs)

eval <- mgSearchOrder(data=crabs, perturbations=NULL, maxCategories=4, nodeCategories=c(5, 3, 4, 4, 4), maxParentSet=2, maxComplexity=0, nodeOrder=c(1,2,3,5,4), emIterations = 10, selectMode="BIC", echo=TRUE)
eval
netlist <- eval@nets
compl<-rep(0, length(netlist))
log.lik<-rep(0, length(netlist))
for(i in 1:length(netlist)) {
	compl[i] <- netlist[[i]]@complexity
	log.lik[i] <- netlist[[i]]@likelihood
}

bnet <- cnFindBIC(eval, nrow(crabs))
bnet@meta <- "BIC selection"
anet <- cnFindAIC(eval)
anet@meta <- "AIC selection"

anet <- netlist[[12]]
cnPlot(anet)

cnDot(eval@nets, "crab")

vliks <- mgNodeLoglik(anet, 1:anet@numnodes, as.data.frame(crabs[,]))
vliks
anet

## plot likelihood vs complexity for the resulting list of networks
plot(compl, log.lik, main="A model selection curve.")
abline(v=anet@complexity,lty=2,col="red")
abline(v=bnet@complexity,lty=3,col="blue")

bst <- cnFind(eval, 121)
bst ## cmplexity = 98, good choice, the smallest for which satell has two parents - THIS IS THE CRITERION !!!

## WE SELECT REGRESSION OF: satell ~ weight + color
bst@meta <- "Final selection"
cnPlot(bst)

ps <- mgSamples(bst, 2000)
par(mfrow=c(2,1))
hist(ps$satell[ps$color<2.8], 50)
hist(ps$satell[ps$color>4], 50)
     
cnet <- cnReorderNodes(bnet,c(1,5,2,4,3))
ps <- mgPredict(cnet, ps)

ps <- crabs
ps$satell <- rep(NA, length(ps$satell))
ps <- mgPredict(bst, ps)
par(mfrow=c(1,2))
hist(ps$satell, 50)
hist(crabs$satell,50)

bnet <- cnFind(eval, 121)
mcov <- matrix(rep(1,bnet@numnodes*bnet@numnodes), nrow=bnet@numnodes)
for(i in 1:(bnet@numnodes-1))
  for(j in i:bnet@numnodes) {
    mcov[i,j] <- cov(ps[,i], ps[,j])
    mcov[j,i] <- mcov[i,j]
    }
mcov
mcov <- matrix(rep(1,bnet@numnodes*bnet@numnodes), nrow=bnet@numnodes)
for(i in 1:(bnet@numnodes-1))
  for(j in i:bnet@numnodes) {
    mcov[i,j] <- cov(crabs[,i], crabs[,j])
    mcov[j,i] <- mcov[i,j]
    }
mcov

# If you haven't already done so, convert weight to kilograms and
# the color scores to 1-4, instead of 2-5:
crabs$weight <- crabs$weight/1000
crabs$color <- crabs$color - 1
crabs$colfac <- factor(crabs$color) # Create a factor version of color
levels(crabs$colfac)  # Level "1" is "first"
# To make levels 1-3 get the dummy variables (like SAS), move
# level 4 to the front (make it the "reference level") using the
# relevel() function:
crabs$colfac <- relevel(crabs$colfac,ref="4")
levels(crabs$colfac)  # Now level "4" is "first"

## Fit the model for part (a)
crabs.pll.a <- glm(satell ~ colfac + weight, family=poisson(), data=crabs)
summary(crabs.pll.a)
# For part (b), prediction for medium light crabs (color=1):
predict(crabs.pll.a, newdata=data.frame(colfac="1",weight=2.44), type="response")
# Or by hand (showing you all the steps involved, but only the
# last one is really needed):
coef(crabs.pll.a)
coef(crabs.pll.a) * c(1,1,0,0,2.44)
sum(coef(crabs.pll.a) * c(1,1,0,0,2.44))
exp(sum(coef(crabs.pll.a) * c(1,1,0,0,2.44)))
# For part (b), dark crabs (color=4):
predict(crabs.pll.a, newdata=data.frame(colfac="4",weight=2.44), type="response")
exp(sum(coef(crabs.pll.a) * c(1,0,0,0,2.44)))
# Do the test for part (c):

## First fit the reduced model without color:
crabs.pll.c <- glm(satell ~ weight, family=poisson(), data=crabs)
anova(crabs.pll.c,crabs.pll.a,test="Chisq")

crabs.pll.c <- glm(satell ~ color, family=poisson(), data=crabs)
crabs.pll.cw <- glm(satell ~ color + weight, family=poisson(), data=crabs)
crabs.pll.cww <- glm(satell ~ color + weight + width, family=poisson(), data=crabs)
anova(crabs.pll.c, crabs.pll.cw, test="Chisq")
anova(crabs.pll.c, crabs.pll.cw, crabs.pll.cww, test="Chisq")
anova(crabs.pll.cw, crabs.pll.cww, test="Chisq")


