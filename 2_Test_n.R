source("Functions_plotting.R")
source("Functions_JoineRML.R")
source("Functions_JMbayes2.R")
source("Functions_PCcox.R")
source("Functions_lcmm.R")
source("Functions_JM.R")
library(foreach)
library(doSNOW)

# Load simulated data
load("../Data/Simdata.RData")
simdata <- simdata.n

# Data in long format, including all time points.
simdata <- simdata[,c("feature","subject","time","intensity","event_time","group")]
colnames(simdata)[6] <- "status"
simdata$id <- as.numeric(as.factor(simdata[,"subject"]))
simdata$status <- ifelse(simdata[,"status"]=="case",1,0)

features <- grep("tp15",unique(simdata$feature), value=TRUE)

# Setup parallel calculations
cl <- makeCluster(16, type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(min=0, max=length(features), style=3)
progress <- function(g) setTxtProgressBar(pb, g)
opts <- list(progress=progress)

aucs <- foreach(i=1:length(features), .combine=rbind, .multicombine=TRUE, .packages=c("joineRML", "partlyconditional", "pROC"), .options.snow=opts) %dopar% {
  setTxtProgressBar(pb, i)
  
  # Select feature
  feature <- features[i]
  sel <- which(simdata$feature==feature)
  data.train <- simdata[sel,-1]
  
  # Data including only the first time point from each individual.
  ids <- unique(data.train[,"id"])
  data.train.id <- NULL
  for(i in 1:length(ids)) data.train.id<- rbind(data.train.id, data.train[data.train[,"id"]==ids[i],][1,])
  
  # Landmarks and horizons to test
  landmark <- c(5,10,15)
  horizon <- c(1:15)
  
  # JoineRML
  model.joineRML <- train.joineRML(data.train)
  predictions.joineRML <- pred.joineRML(data.train.id)
  auc.5 <- get.auc(eval.time=20, ladmark.time=5, predictions=predictions.joineRML, true.response=data.train.id[,"status"])
  auc.10 <- get.auc(eval.time=20, ladmark.time=10, predictions=predictions.joineRML, true.response=data.train.id[,"status"])
  auc.15 <- get.auc(eval.time=20, ladmark.time=15, predictions=predictions.joineRML, true.response=data.train.id[,"status"])
  auc.joineRML <- c(method="JoineRML",feature=feature,auc.5,auc.10,auc.15)
  
  # JMbayes2
  #model.JMbayes2 <- train.JMbayes2(data.train, data.train.id)
  #predictions.JMbayes2 <- pred.JMbayes2(data.train.id)
  #auc.5 <- get.auc(eval.time=20, ladmark.time=5, predictions=predictions.JMbayes2, true.response=data.train.id[,"status"])
  #auc.10 <- get.auc(eval.time=20, ladmark.time=10, predictions=predictions.JMbayes2, true.response=data.train.id[,"status"])
  #auc.15 <- get.auc(eval.time=20, ladmark.time=15, predictions=predictions.JMbayes2, true.response=data.train.id[,"status"])
  #auc.JMbayes2 <- c(method="JMbayes2",feature=feature,auc.5,auc.10,auc.15)
  
  # PCcox
  model.pcCox <- train.pcCox(data.train)
  predictions.pcCox <- pred.pcCox(data.train.id)
  auc.5 <- get.auc(eval.time=20, ladmark.time=5, predictions=predictions.pcCox, true.response=data.train.id[,"status"])
  auc.10 <- get.auc(eval.time=20, ladmark.time=10, predictions=predictions.pcCox, true.response=data.train.id[,"status"])
  auc.15 <- get.auc(eval.time=20, ladmark.time=15, predictions=predictions.pcCox, true.response=data.train.id[,"status"])
  auc.PCcox <- c(method="PCcox",feature=feature,auc.5,auc.10,auc.15)
  
  # lcmm
  #model.lcmm <- train.lcmm(data.train)
  #predictions.lcmm <- pred.lcmm(data.train.id)
  #auc.5 <- get.auc(eval.time=20, ladmark.time=5, predictions=predictions.lcmm, true.response=data.train.id[,"status"])
  #auc.10 <- get.auc(eval.time=20, ladmark.time=10, predictions=predictions.lcmm, true.response=data.train.id[,"status"])
  #auc.15 <- get.auc(eval.time=20, ladmark.time=15, predictions=predictions.lcmm, true.response=data.train.id[,"status"])
  #auc.lcmm <- c(method="lcmm",feature=feature,auc.5,auc.10,auc.15)
  
  # JM
  #model.JM <- train.JM(data.train, data.train.id)
  #predictions.JM <- pred.JM(data.train.id)
  #auc.5 <- get.auc(eval.time=20, ladmark.time=5, predictions=predictions.JM, true.response=data.train.id[,"status"])
  #auc.10 <- get.auc(eval.time=20, ladmark.time=10, predictions=predictions.JM, true.response=data.train.id[,"status"])
  #auc.15 <- get.auc(eval.time=20, ladmark.time=15, predictions=predictions.JM, true.response=data.train.id[,"status"])
  #auc.JM <- c(method="JM",feature=feature,auc.5,auc.10,auc.15)
  
  rbind(auc.joineRML, auc.PCcox)
}
close(pb)
stopCluster(cl)
aucs <- data.frame(aucs)
for(i in 3:5) aucs[,i] <- as.numeric(aucs[,i])

methods <- c("JoineRML","JMbayes2","PCcox","lcmm","JM")
colors <- colorspace::qualitative_hcl(length(methods), "Pastel 1")
names(colors) <- methods

pdf("../Results/Subjects.pdf", width=7, height=2)
par(mfrow=c(1,4))

plot(NULL, xlim=c(1,10), ylim=c(0.4,1), xlab="Difficulty", ylab="AUC", bty="n")
abline(h=0.5, lty=2, col="black")
lines(10:1, aucs[which(aucs$method=="JoineRML"),5], lwd=2, col=colors["JoineRML"])
lines(10:1, aucs[which(aucs$method=="PCcox"),5], lwd=2, col=colors["PCcox"])
mtext("landmark = 15")

plot(NULL, xlim=c(1,10), ylim=c(0.4,1), xlab="Difficulty", ylab="AUC", bty="n")
abline(h=0.5, lty=2, col="black")
lines(10:1, aucs[which(aucs$method=="JoineRML"),4], lwd=2, col=colors["JoineRML"])
lines(10:1, aucs[which(aucs$method=="PCcox"),4], lwd=2, col=colors["PCcox"])
mtext("landmark = 10")

plot(NULL, xlim=c(1,10), ylim=c(0.4,1), xlab="Difficulty", ylab="AUC", bty="n")
abline(h=0.5, lty=2, col="black")
lines(10:1, aucs[which(aucs$method=="JoineRML"),3], lwd=2, col=colors["JoineRML"])
lines(10:1, aucs[which(aucs$method=="PCcox"),3], lwd=2, col=colors["PCcox"])
mtext("landmark = 5")

frame()
legend("center", legend=methods, col=colors, lty=1, lwd=2, bty="n")

dev.off()

