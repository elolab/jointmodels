library(plyr)
library(foreach)
source("Functions_lcmm.R")
source("Functions_JoineRML.R")
source("Functions_PCcox.R")
source("Functions_plotting.R")

# Load data
load("../Data/Liu_et_al_data.RDAta")
data.meta <- data.frame(des_matrix, stringsAsFactors=FALSE); rm(des_matrix)
data <- data[,data.meta$Sample.Names]

dim(data)

# albumin: ALB
# IgG: IGHG1
# antitrypsin: SERPINA1
# IgA: IGHA1
# transferrin: TF
# haptoglobin HP
# fibrinogen FGA, FGB, FGG
# alpha2-macroglobulin: A2M
# alpha1-acid glycoprotein ORM
# IgM: IGHM
# apolipoprotein AI: APOA1
# apolipoprotein AII: APOA2
# complement C3: C3
# transthyretin: TTR or TBPA

depleted <- c("ALB","IGHG1","SERPINA1","IGHA1","TF","HP;HPR","FGA","FGB","FGG","A2M","ORM","IGHM","APOA1","APOA2","C3","TTR","TBPA")
sum(!(rownames(data) %in% depleted))

# Prepare metadata
data.meta <- data.meta[,2:4]
colnames(data.meta) <- c("group","time","id")
data.meta$time <- as.numeric(data.meta$time)
data.meta$id <- as.numeric(data.meta$id)
data.meta$status <- ifelse(data.meta$group=="TD",1,0)
data.meta$group <- ifelse(data.meta$group=="TD","case","control")
data.meta$subject <- paste(data.meta$group, data.meta$id, sep="_")
data.meta$event_time <- foreach(id=unique(data.meta$id), .combine=c) %do% {
  sel <- which(data.meta$id==id)
  rep(max(data.meta$time[sel]),length(sel))
}

# Metadata by id
data.meta.id <- ddply(data.meta, .(group, subject, id), numcolwise(mean))

# Set landmarks and horizons
landmark<- c(5,6,7,8)
horizon<- c(1:15)

# Function to supress cat output
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}

#--------------------------------------------------
# lcmm
#--------------------------------------------------

pb <- txtProgressBar(min=0, max=nrow(data), style=3)
lcmm.out <- foreach(i=1:nrow(data)) %do% {
  setTxtProgressBar(pb, i)
  tryCatch({
    data.train <- cbind(data.meta, intensity=data[i,])
    model.lcmm <- quiet(train.lcmm.randint(data.train))
    model.lcmm
  }, error=function(e) NA)
}
close(pb)
save(lcmm.out, file="../Data/Liu_lcmm.RData")

lcmm.ng <- c(t(unlist(sapply(lcmm.out, function(x) x["ng"]))))
table(lcmm.ng)

lcmm.min <- foreach(i=1:length(lcmm.out), .combine=c) %do% {
  tryCatch({
    min(postprob(lcmm.out[[i]])[[4]]["N",])
  }, error=function(e) NA)
}
table(lcmm.min)

sel <- which(lcmm.ng==2 & lcmm.min>5)
lcmm.genes <- rownames(data)[sel]
lcmm.aucs <- foreach(i=1:length(sel), .combine=rbind) %do% {
  tryCatch({
    data.test <- cbind(data.meta, intensity=data[sel[i],])
    predictions.lcmm <- pred.lcmm(lcmm.out[[sel[i]]], data.test=data.test, data.test.id=data.meta.id, landmark=landmark, horizon=horizon)
    sapply(landmark, function(x) get.auc(eval.time=13, landmark.time=x, predictions=predictions.lcmm, true.response=data.meta.id[,"status"]))
  }, error=function(e) rep(NA,4))
}

genes.lcmm <- rownames(data)[sel][which(apply(lcmm.aucs[,3:4],1,function(x) all(x>0.75)))]
genes.lcmm <- genes.lcmm[which(!(genes.lcmm %in% depleted))]

pdf("../Results/Liu_lcmm.pdf", width=4, height=3)
for(g in genes.lcmm) {
  data.train <- cbind(data.meta, intensity=data[g,])
  plotProfile(data.train, individuals=TRUE, scale=FALSE, main=g, landmark=landmark)
}
dev.off()

#--------------------------------------------------
# joineRML
#--------------------------------------------------

pb <- txtProgressBar(min=0, max=nrow(data), style=3)
joinerml.out <- foreach(i=1:nrow(data)) %do% {
  setTxtProgressBar(pb, i)
  tryCatch({
    data.train <- cbind(data.meta, intensity=data[i,])
    data.train <- data.train[complete.cases(data.train),]
    model.joinerml <- quiet(suppressMessages(train.joineRML.randint(data.train)))
    model.joinerml
  }, error=function(e) NA)
}
close(pb)
save(joinerml.out, file="../Data/Liu_joineRML.RData")

joinerml.p <- foreach(i=1:length(joinerml.out), .combine=c) %do% {
  tryCatch({
    summary(joinerml.out[[i]])$coefs.surv[4]
  }, error=function(e) NA)
}
sel <- which(joinerml.p<0.05)
joinerml.aucs <- foreach(i=1:length(sel), .combine=rbind) %do% {
  tryCatch({
    data.test <- cbind(data.meta, intensity=data[sel[i],])
    data.test <- data.test[complete.cases(data.test),]
    predictions.joinerml <- pred.joineRML(joinerml.out[[sel[i]]], data.test=data.test, data.test.id=data.meta.id, landmark=landmark, horizon=horizon)
    sapply(landmark, function(x) get.auc(eval.time=13, landmark.time=x, predictions=predictions.joinerml, true.response=data.meta.id[,"status"]))
  }, error=function(e) rep(NA,4))
}
joinerml.nmin <- foreach(i=1:length(sel), .combine=c) %do% {
  tryCatch({
    data.test <- cbind(data.meta, intensity=data[sel[i],])
    data.test.id <- ddply(data.test, .(group, subject, id), numcolwise(mean))
    min(table(data.test.id$status[which(!is.na(data.test.id$intensity))]))
  }, error=function(e) NA)
}

genes.joinerml <- rownames(data)[sel][which(apply(joinerml.aucs[,3:4],1,function(x) all(x>0.75)) & joinerml.nmin>5)]
genes.joinerml <- genes.joinerml[which(!(genes.joinerml %in% depleted))]

pdf("../Results/Liu_joineRML.pdf", width=4, height=3)
for(g in genes.joinerml) {
  data.train <- cbind(data.meta, intensity=data[g,])
  plotProfile(data.train, individuals=TRUE, scale=FALSE, main=g, landmark=landmark)
}
dev.off()

#--------------------------------------------------
# PCCox
#--------------------------------------------------

pb <- txtProgressBar(min=0, max=nrow(data), style=3)
pccox.out <- foreach(i=1:nrow(data)) %do% {
  setTxtProgressBar(pb, i)
  tryCatch({
    data.train <- cbind(data.meta, intensity=data[i,])
    model.pccox <- quiet(suppressMessages(train.pcCox(data.train)))
    model.pccox
  }, error=function(e) NA)
}
close(pb)
save(pccox.out, file="../Data/Liu_PCCox.RData")

pccox.p <- foreach(i=1:length(joinerml.out), .combine=c) %do% {
  tryCatch({
    summary(pccox.out[[i]]$model.fit)$coefficients["intensity","Pr(>|z|)"]
  }, error=function(e) NA)
}
sel <- which(pccox.p<0.01)
pccox.aucs <- foreach(i=1:length(sel), .combine=rbind) %do% {
  tryCatch({
    data.test <- cbind(data.meta, intensity=data[sel[i],])
    predictions.pccox <- pred.pcCox(pccox.out[[sel[i]]], data.test=data.test, data.test.id=data.meta.id, landmark=landmark, horizon=horizon, max.time=max(data.train$event_time[data.train$status==1]))
    sapply(landmark, function(x) get.auc(eval.time=13, landmark.time=x, predictions=predictions.pccox, true.response=data.meta.id[,"status"]))
  }, error=function(e) rep(NA,4))
}
pccox.nmin <- foreach(i=1:length(sel), .combine=c) %do% {
  tryCatch({
    data.test <- cbind(data.meta, intensity=data[sel[i],])
    data.test.id <- ddply(data.test, .(group, subject, id), numcolwise(mean))
    min(table(data.test.id$status[which(!is.na(data.test.id$intensity))]))
  }, error=function(e) NA)
}

genes.pccox <- rownames(data)[sel][which(apply(pccox.aucs[,3:4],1,function(x) all(x>0.75)) & pccox.nmin>5)]
genes.pccox <- genes.pccox[which(!(genes.pccox %in% depleted))]

pdf("../Results/Liu_PCCox.pdf", width=4, height=3)
for(g in genes.pccox) {
  data.train <- cbind(data.meta, intensity=data[g,])
  plotProfile(data.train, individuals=TRUE, scale=FALSE, main=g, landmark=landmark)
}
dev.off()

#--------------------------------------------------
# Overlap
#--------------------------------------------------

list.venn <- list("lcmm"=genes.lcmm, "joineRML"=genes.joinerml, "PCCox"=genes.pccox)

library(ggvenn)
pdf("../Results/Liu_venn.pdf", width=4, height=4)
ggvenn(list.venn, fill_color=c("#A8E1BF","#FFC5D0","#E2D4A8"), stroke_size=0.5,
       set_name_size=4,text_size = 1, show_elements=TRUE, show_percentage=FALSE)
dev.off()

genes.pccox
genes.joinerml
genes.lcmm

genes <- unique(c(genes.pccox, genes.joinerml, genes.lcmm))

aucs.lcmm <- foreach(i=1:length(genes), .combine=rbind) %do% {
  tryCatch({
    data.train <- cbind(data.meta, intensity=data[genes[i],])
    model.lcmm <- quiet(train.lcmm.randint(data.train))
    predictions.lcmm <- pred.lcmm(model.lcmm, data.test=data.train, data.test.id=data.meta.id, landmark=landmark, horizon=horizon)
    sapply(landmark, function(x) get.auc(eval.time=13, landmark.time=x, predictions=predictions.lcmm, true.response=data.meta.id[,"status"]))
  }, error=function(e) rep(NA,4))
}
rownames(aucs.lcmm) <- genes
colnames(aucs.lcmm) <- landmark

aucs.joinerml <- foreach(i=1:length(genes), .combine=rbind) %do% {
  tryCatch({
    data.train <- cbind(data.meta, intensity=data[genes[i],])
    data.train <- data.train[complete.cases(data.train),]
    model.joinerml <- quiet(train.joineRML.randint(data.train))
    predictions.joinerml <- pred.joineRML(model.joinerml, data.test=data.train, data.test.id=data.meta.id, landmark=landmark, horizon=horizon)
    sapply(landmark, function(x) get.auc(eval.time=13, landmark.time=x, predictions=predictions.joinerml, true.response=data.meta.id[,"status"]))
  }, error=function(e) rep(NA,4))
}
rownames(aucs.joinerml) <- genes
colnames(aucs.joinerml) <- landmark

aucs.pccox <- foreach(i=1:length(genes), .combine=rbind) %do% {
  tryCatch({
    data.train <- cbind(data.meta, intensity=data[genes[i],])
    model.pccox <- quiet(train.pcCox(data.train))
    predictions.pccox <- pred.pcCox(model.pccox, data.test=data.train, data.test.id=data.meta.id, landmark=landmark, horizon=horizon, max.time=max(data.train$event_time[data.train$status==1]))
    sapply(landmark, function(x) get.auc(eval.time=13, landmark.time=x, predictions=predictions.pccox, true.response=data.meta.id[,"status"]))
  }, error=function(e) rep(NA,4))
}
rownames(aucs.pccox) <- genes
colnames(aucs.pccox) <- landmark

aucs <- cbind(aucs.lcmm, aucs.pccox, aucs.joinerml)
aucs.avg <- (aucs.lcmm+aucs.pccox+aucs.joinerml)/3

library(pheatmap)
breaks <- seq(0.5, 1, by=0.05)
col <- colorRampPalette(c("white", "royalblue", "firebrick"))(length(breaks)-1)
#marks <- matrix("", nrow=nrow(aucs), ncol=ncol(aucs))
#marks[aucs>0.75] <- "*"

pdf("../Results/Liu_heatmap.pdf", width=4, height=3.3)
pheatmap(aucs, cluster_cols=FALSE, col=col, breaks=breaks, gaps_col=c(4,8))
#display_numbers=marks, number_color="white"
dev.off()

genes <- c("ALPL","ITGA5","CHST3","SAA1","PEF1","MATN2","HLA-B.2.2")
pdf("../Results/Liu_barplot.pdf", width=4, height=6)
par(mfrow=c(3,1))
barplot(t(aucs.lcmm[genes,]), beside=TRUE, ylim=c(0,1), ylab="AUC", col=colorRampPalette(c("white", "royalblue"))(4), las=2)
abline(h=0.75, lty=2)
barplot(t(aucs.pccox[genes,]), beside=TRUE, ylim=c(0,1), ylab="AUC", col=colorRampPalette(c("white", "royalblue"))(4), las=2)
abline(h=0.75, lty=2)
barplot(t(aucs.joinerml[genes,]), beside=TRUE, ylim=c(0,1), ylab="AUC", col=colorRampPalette(c("white", "royalblue"))(4), las=2)
abline(h=0.75, lty=2)
dev.off()

g <- "MGAT1"
data.train <- cbind(data.meta, intensity=data[g,])
plotProfile(data.train, individuals=TRUE, scale=FALSE, main=g)

