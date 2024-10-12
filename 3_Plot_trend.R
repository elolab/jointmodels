library(foreach)

# Gather results
files <- list.files("../Results/Simulation/", pattern=".RData")
results <- foreach(f=files, .combine=rbind) %do% {
  load(paste("../Results/Simulation/",f,sep=""))
  if(is.null(results.joineRML$AUC)) {results.joineRML$AUC <- c(NA,NA,NA)}
  if(is.null(results.lcmm$AUC)) {results.lcmm$AUC <- c(NA,NA,NA)}
  if(is.null(results.pcCox$AUC)) {results.pcCox$AUC <- c(NA,NA,NA)}
  if(is.null(results.JM$AUC)) {results.JM$AUC <- c(NA,NA,NA)}
  if(is.null(results.Cox$AUC)) {results.Cox$AUC <- c(NA,NA,NA)}
  
  if(is.null(results.joineRML$CI)) {results.joineRML$CI <- matrix(NA,3,3)}
  if(is.null(results.lcmm$CI)) {results.lcmm$CI <- matrix(NA,3,3)}
  if(is.null(results.pcCox$CI)) {results.pcCox$CI <- matrix(NA,3,3)}
  if(is.null(results.JM$CI)) {results.JM$CI <- matrix(NA,3,3)}
  if(is.null(results.Cox$CI)) {results.Cox$CI <- matrix(NA,3,3)}
  
  if(is.null(results.joineRML$p)) {results.joineRML$p <- NA}
  if(is.null(results.lcmm$p)) {results.lcmm$p <- NA}
  if(is.null(results.pcCox$p)) {results.pcCox$p <- NA}
  if(is.null(results.JM$p)) {results.JM$p <- NA}
  if(is.null(results.Cox$p)) {results.Cox$p <- NA}
  
  if(is.null(results.joineRML$time)) {results.joineRML$time <- NA}
  if(is.null(results.lcmm$time)) {results.lcmm$time <- NA}
  if(is.null(results.pcCox$time)) {results.pcCox$time <- NA}
  if(is.null(results.JM$time)) {results.JM$time <- NA}
  if(is.null(results.Cox$time)) {results.Cox$time <- NA}
  
  auc <- data.frame(rbind(results.joineRML$AUC, results.lcmm$AUC, results.pcCox$AUC, results.JM$AUC, results.Cox$AUC))
  ci.low <- data.frame(rbind(results.joineRML$CI[1,], results.lcmm$CI[1,], results.pcCox$CI[1,], results.JM$CI[1,], results.Cox$CI[1,]))
  ci.high <- data.frame(rbind(results.joineRML$CI[3,], results.lcmm$CI[3,], results.pcCox$CI[3,], results.JM$CI[3,], results.Cox$CI[3,]))
  res <- cbind(auc,ci.low,ci.high)
  res$p <- c(results.joineRML$p, results.lcmm$p, results.pcCox$p, results.JM$p, results.Cox$p)
  res$time <- c(results.joineRML$time, results.lcmm$time, results.pcCox$time, results.JM$time, results.Cox$time)
  res$method <- c("JoineRML","lcmm","pcCox","JM","Cox")
  res$file <- f
  res
}
results <- cbind(results, t(sapply(results$file, function(x) unlist(strsplit(x,"_|\\."))))[,2:4])
colnames(results)[14:16] <- c("subj","trend","var")
colnames(results)[1:3] <- c("auc5","auc10","auc15")
colnames(results)[4:6] <- c("auc5.ci.low","auc10.ci.low","auc15.ci.low")
colnames(results)[7:9] <- c("auc5.ci.high","auc10.ci.high","auc15.ci.high")

writexl::write_xlsx(results, path="../Results/Trend.xlsx")

methods <- c("JoineRML","pcCox","lcmm","JM","Cox")
colors <- colorspace::qualitative_hcl(length(methods), "Pastel 1")
colors <- adjustcolor(colors, offset=c(-0.3, -0.3, -0.3, 0))
names(colors) <- methods

landmarks <- c(3,2,1)
names(landmarks) <- c("15","10","5")

pdf("../Results/Trend.pdf", width=8, height=6)
par(mfrow=c(3,4))
for(type in 1:3) {
  for(var in unique(results$var)) {
    for (l in 1:length(landmarks)) {
      plot(NULL, xlim=c(1,5), ylim=c(0.4,1), xlab="Difficulty", ylab="AUC", bty="n")
      abline(h=0.5, lty=2, col="black")
      for(m in methods) {
        if(type==1 | type==3) {
          y <- results[which(results$method==m & results$var==var & results$subj=="subj10"),landmarks[l]]
          lines(5:1, y, lwd=1, lty=2, col=colors[m])
          points(5:1, y, pch=19, col=colors[m])
        }
        if(type==1 | type==2) {
          y <- results[which(results$method==m & results$var==var & results$subj=="subj50"),landmarks[l]]
          lines(5:1, y, lwd=1, col=colors[m])
          points(5:1, y, pch=19, col=colors[m])
        }
      }
      mtext(paste("landmark =", names(landmarks)[l]))
    }
    frame()
    legend("center", legend=methods, col=colors, lty=1, lwd=2, bty="n")
  }
}
dev.off()

pdf("../Results/Trend_CI.pdf", width=8, height=6)
par(mfrow=c(3,4))
for(type in 1:3) {
  for(var in unique(results$var)) {
    for (l in 1:length(landmarks)) {
      plot(NULL, xlim=c(1,5), ylim=c(0.4,1), xlab="Difficulty", ylab="AUC", bty="n")
      abline(h=0.5, lty=2, col="black")
      for(m in methods) {
        if(type==1 | type==3) {
          y <- results[which(results$method==m & results$var==var & results$subj=="subj10"),landmarks[l]]
          y.low <- results[which(results$method==m & results$var==var & results$subj=="subj10"),landmarks[l]+3]
          y.high <- results[which(results$method==m & results$var==var & results$subj=="subj10"),landmarks[l]+6]
          polygon(c(5:1, 1:5), c(y.low, rev(y.high)), col=adjustcolor(colors[m], alpha.f=0.2), border=NA)
          lines(5:1, y, lwd=1, lty=2, col=colors[m])
          points(5:1, y, pch=19, col=colors[m])
        }
        if(type==1 | type==2) {
          y <- results[which(results$method==m & results$var==var & results$subj=="subj50"),landmarks[l]]
          y.low <- results[which(results$method==m & results$var==var & results$subj=="subj50"),landmarks[l]+3]
          y.high <- results[which(results$method==m & results$var==var & results$subj=="subj50"),landmarks[l]+6]
          polygon(c(5:1, 1:5), c(y.low, rev(y.high)), col=adjustcolor(colors[m], alpha.f=0.2), border=NA)
          lines(5:1, y, lwd=1, col=colors[m])
          points(5:1, y, pch=19, col=colors[m])
        }
      }
      mtext(paste("landmark =", names(landmarks)[l]))
    }
    frame()
    legend("center", legend=methods, col=colors, lty=1, lwd=2, bty="n")
  }
}
dev.off()

pdf("../Results/Trend_pvalues.pdf", width=2, height=6)
par(mfrow=c(3,1))
for(type in 1:3) {
  for(var in unique(results$var)) {
    plot(NULL, xlim=c(1,5), ylim=c(0,30), xlab="Difficulty", ylab="-log10(p)", bty="n")
    abline(h=-log10(c(0.1, 0.01, 0.001)), lty=1, col="grey")
    for(m in methods) {
      if(type==1 | type==3) {
        lines(5:1, -log10(results[which(results$method==m & results$var==var & results$subj=="subj10"),"p"]), lwd=1, lty=2, col=colors[m])
        points(5:1, -log10(results[which(results$method==m & results$var==var & results$subj=="subj10"),"p"]), pch=19, col=colors[m])
      }
      if(type==1 | type==2) {
        lines(5:1, -log10(results[which(results$method==m & results$var==var & results$subj=="subj50"),"p"]), lwd=1, col=colors[m])
        points(5:1, -log10(results[which(results$method==m & results$var==var & results$subj=="subj50"),"p"]), pch=19, col=colors[m])
      }
    }
    #legend("topleft", legend=methods, col=colors, lty=1, lwd=2, bty="n")
  }
}
dev.off()


pdf("../Results/Trend_pvalues_comparison.pdf", width=2.5*4, height=2.5)
par(mfrow=c(1,4))

x <- -log10(results$p[which(results$method=="JoineRML" & results$subj=="subj50")])
y <- -log10(results$p[which(results$method=="pcCox" & results$subj=="subj50")])
plot(x, y, xlim=c(0,max(c(x,y))), ylim=c(0,max(c(x,y))), pch=21, bg="royalblue", bty="n", xlab="JoineRML -log10(p)", ylab="pcCox -log10(p)", main="50 subjects")
abline(0,1, lty=2)

x <- -log10(results$p[which(results$method=="JoineRML" & results$subj=="subj10")])
y <- -log10(results$p[which(results$method=="pcCox" & results$subj=="subj10")])
plot(x, y, xlim=c(0,max(c(x,y))), ylim=c(0,max(c(x,y))), pch=21, bg="royalblue", bty="n", xlab="JoineRML -log10(p)", ylab="pcCox -log10(p)", main="10 subjects")
abline(0,1, lty=2)

x <- -log10(results$p[which(results$method=="JoineRML" & results$subj=="subj10")])
y <- -log10(results$p[which(results$method=="JoineRML" & results$subj=="subj50")])
plot(x, y, xlim=c(0,max(c(x,y))), ylim=c(0,max(c(x,y))), pch=21, bg="royalblue", bty="n", xlab="10 subjects", ylab="50 subjects", main="JoineRML -log10(p)")
abline(0,1, lty=2)

x <- -log10(results$p[which(results$method=="pcCox" & results$subj=="subj10")])
y <- -log10(results$p[which(results$method=="pcCox" & results$subj=="subj50")])
plot(x, y, xlim=c(0,max(c(x,y))), ylim=c(0,max(c(x,y))), pch=21, bg="royalblue", bty="n", xlab="10 subjects", ylab="50 subjects", main="pcCox -log10(p)")
abline(0,1, lty=2)

dev.off()

# Plot running times
methods <- unique(results$method)[c(2,1,4,3,5)]
results.time <- results[,c(5,6,8,9,10)]
results.time <- tidyr::spread(results.time, key=c("method"), value="time")
results.time <- results.time[complete.cases(results.time),]

pdf("../Results/Trend_running_time.pdf", width=5, height=5)
par(bty="n", mfrow=c(1,2))
for(s in c("subj10","subj50")) {
  boxplot(results.time[which(results.time$subj==s),methods], outline=FALSE, xaxt="n", ylim=c(0,30), ylab="Seconds")
  mtext(s)
  axis(1, at=1:5, labels=methods, las=2, tick=FALSE)
}
dev.off()

# Covert to long format by landmarks
results.auc <- reshape2::melt(results, id.vars=colnames(results)[-c(1:3)])

# Test differences between Cox and others
methods <- unique(results.auc$method)
wilcox <- foreach(m=1:4, .combine=cbind) %do% {
  a <- which(results.auc$method=="Cox")
  b <- which(results.auc$method==methods[m])
  suppressWarnings(wilcox.test(results.auc[a,"value"], results.auc[b,"value"], paired=TRUE, alternative="less")$p.value)
}
names(wilcox) <- methods[1:4]
max(wilcox)

# Test differences between 10 and 50
a <- which(results.auc$subj=="subj10")
b <- which(results.auc$subj=="subj50")
suppressWarnings(wilcox.test(results.auc[a,"value"], results.auc[b,"value"], paired=TRUE, alternative="less")$p.value)

# Test differences to joineRML and lcmm
methods <- unique(results.auc$method)
wilcox <- foreach(m=2:5, .combine=cbind) %do% {
  a <- which(results.auc$method=="JoineRML")
  b <- which(results.auc$method==methods[m])
  suppressWarnings(wilcox.test(results.auc[a,"value"], results.auc[b,"value"], paired=TRUE, alternative="greater")$p.value)
}
names(wilcox) <- methods[2:5]
max(wilcox)

# Test differences between Cox and others (highest variance)
methods <- unique(results$method)
wilcox <- foreach(m=1:4, .combine=cbind) %do% {
  foreach(auc=c("auc5","auc10","auc15"), .combine=c) %do% {
    a <- which(results$method=="Cox" & results$var=="var3")
    b <- which(results$method==methods[m] & results$var=="var3")
    suppressWarnings(wilcox.test(results[a,auc], results[b,auc], paired=TRUE, alternative="less")$p.value)
  }
}
colnames(wilcox) <- methods[1:4]
rownames(wilcox) <- c("auc5","auc10","auc15")
max(wilcox)

