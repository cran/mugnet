library(mugnet)

data("letters")

colnames(lets) <- c(
                    "lettr", "x-box", "y-box", "width", "high", "onpix",
                    "x-bar", "y-bar", "x2bar", "y2bar", "xybar", "x2ybr",
                    "xy2br", "x-ege", "xegvy", "y-ege", "yegvx")
clas <- levels(lets$lettr)

nsamples <- nrow(lets)
numnodes <- ncol(lets)
train <- lets[1:10000,]
valid <- lets[1001:15000,]
test <- lets[15001:20000,]

evala <- vector("list", length(clas))
for(i in 1:length(clas)) {
  cat("\n", clas[i], "\n")
  letdata <- train[train$lettr==clas[i], 2:numnodes]
  evala[[i]] <- mgSearchOrder(letdata, NULL, maxCategories = 4, NULL, maxParentSet=2, maxComplexity=0, emIter=10, echo=TRUE)
}

net <- vector("list", length(clas))
targetcomplx <- 400
for(i in 1:length(clas)) {
  net[[i]] <- cnFind(evala[[i]], targetcomplx)
  net[[i]]@meta <- paste("letter", clas[i])
}

predict <- rep(NA, nrow(test))
for(j in 1:nrow(test)) {
  cat(j, "\n")
  liks <- rep(0, length(clas))
  for(i in 1:length(clas)) {
    ##cat(i," ")
    liks[i] <- mgLoglik(net[[i]], test[j,2:numnodes])
  }
  id <- which(liks == max(liks))
  if(length(id)==1)
    predict[j] <- clas[id]
}
sum(predict==test[,1])/nrow(test)
## complx=300, 82.26% good
## complx=400, 83.96% 
