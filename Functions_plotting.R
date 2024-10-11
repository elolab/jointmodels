library(foreach)

#--------------------------------------------------
# Plotting
#--------------------------------------------------

# Function for loess smoothing
smooth.loess <- function(x, y, n=100){
  lo <- suppressWarnings(loess(y ~ x))
  xl <- seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE), (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))/n)
  pred <- predict(lo,xl, se=TRUE)
  yl <- pred$fit
  return(list(x=xl, y=yl, y.low=pred$fit - qt(0.975,pred$df)*pred$se, y.high=pred$fit + qt(0.975,pred$df)*pred$se))
}

# Function to plot profile
plotProfile <- function(data, individuals=FALSE, scale=FALSE, landmark=NULL, ...) {
  if(!is.null(landmark)) {
    sel <- which(data$time<=max(landmark))
    data <- data[sel,]
  }
  if(scale) {
    means <- tapply(data$intensity, data$id, function(x) mean(x, na.rm=TRUE))
    for(status in unique(data$status)) {
      ids <- unique(data$id[which(data$status==status)])
      means[ids] <- means[ids] - mean(means[ids], na.rm=TRUE)
    }
    ids <- unique(data$id)
    for(id in ids) {
      sel <- which(data$id==id)
      data$intensity[sel] <- data$intensity[sel] - means[id]
    }
  }
  colors <- c("firebrick", "dodgerblue4")
  names(colors) <- c("case","control")
  # Prepare plot
  plot(NULL, xlab="Time", ylab="Value", xlim=range(c(data$time,data$event_time)), ylim=range(data$intensity, na.rm=TRUE), bty="n", ...)
  # Add landmarks
  if(!is.null(landmark)) {
    abline(v=landmark, col=adjustcolor("black", alpha.f=0.5))
  }
  # Plot all values
  points(data$time, data$intensity, pch=19, cex=0.5, col=adjustcolor(colors[data$group],alpha.f=0.25))
  # Plot group averages
  for(group in c("case","control")) {
    sel <- which(data$group==group)
    lo <- smooth.loess(data$time[sel], data$intensity[sel])
    use <- which(!is.na(lo$y))
    lo <- lapply(lo, function(x) x[use])
    polygon(c(lo$x, rev(lo$x)), c(lo$y.low, rev(lo$y.high)), col=adjustcolor(colors[group],alpha.f=0.5), border=NA)
    lines(lo$x, lo$y, col=colors[group], lty=2, lwd=3)
  }
  # Plot events
  data.event <- foreach(s=unique(data$subject), .combine=rbind) %do% {
    sel <- which(data$subject==s)
    sel <- sel[which.min(abs(data$time[sel]-data$event_time[sel][1]))]
    data[sel,]
  }
  points(data.event$event_time, data.event$intensity, pch=ifelse(data.event$group=="case",4,1), col=adjustcolor(colors[data.event$group],alpha.f=0.75))
  # Plot individuals
  if(individuals==TRUE) {
    for(id in unique(data$subject)) {
      sel <- which(data$subject==id)
      if (length(sel)>1)  {
        ord <- order(data$time[sel])
        lines(data$time[sel][ord], data$intensity[sel][ord], lwd=2, col=adjustcolor(colors[data$group[sel][1]],alpha.f=0.25))
      }
      sel <- which(data.event$subject==id)
      lines(x=c(data.event$time[sel], data.event$event_time[sel]),
            y=c(data.event$intensity[sel],data.event$intensity[sel]),
            col=adjustcolor(colors[data.event$group[sel][1]],alpha.f=0.25), lwd=2, lty=3)
    }
  }
  #legend("top", names(colors), col=colors, pch=15, bty="n", horiz=TRUE)
}

#--------------------------------------------------
# Model evaluation
#--------------------------------------------------

plot.boxplots <- function(eval.time, predictions, cases, controls, mtext=NULL) {
  eval <- paste("eval", eval.time, sep=".")
  landmarks <- as.numeric(sapply(strsplit(row.names(predictions[[1]]), "[.]"), function(x) x[2]))
  par(mfrow=c(1,length(landmarks)))
  for(l in 1:length(landmarks)) {
    temp<- unlist(lapply(predictions, function(x) ifelse(!is.null(nrow(x)),x[l,eval],NA)))
    boxplot(list(control=temp[grep(controls,names(predictions))], case=temp[grep(cases,names(predictions))]), col=c("cyan","red"), ylab=paste("Disease risk at",eval.time,"years"), las=2, main=paste("Landmark",landmarks[l],"years"), ylim=c(0,1))
    mtext(mtext)
  }
}

plot.risks <- function(predictions, cases, controls, mtext=NULL) {
  landmarks <- as.numeric(sapply(strsplit(row.names(predictions[[1]]), "[.]"), function(x) x[2]))
  horizons <- as.numeric(sapply(strsplit(colnames(predictions[[1]]), "[.]"), function(x) x[2]))
  par(mfrow=c(2,2))
  for(l in 1:length(landmarks)) {
    risks<- matrix(nrow=length(predictions), ncol=ncol(predictions[[1]]))
    for(i in 1:length(predictions)) {
      if(!is.null(nrow(predictions[[i]]))) risks[i,]<- predictions[[i]][l,]
    }
    row.names(risks) <- names(predictions)
    colnames(risks) <- colnames(predictions[[1]])
    colors<- rep("gray", length(predictions))
    colors[grep(cases,names(predictions))]<- "black"
    plot(horizons, risks[1,], type="b", pch=20, col=colors[1], ylim=c(0,1), xlab="Time", ylab="Disease risk", las=1, main=paste("Landmark",landmarks[l],"years"))
    for(i in 2:length(predictions)) tryCatch(lines(horizons, risks[i,], type="b", pch=20, col=colors[i]), error=function(e) NA)
    mtext(mtext)
  }
}

plot.roc <- function(eval.time, predictions, true.response, mtext=NULL) {
  library(pROC)
  par(mfrow=c(2,2))
  eval <- paste("eval", eval.time, sep=".")
  landmarks <- as.numeric(sapply(strsplit(row.names(predictions[[1]]), "[.]"), function(x) x[2]))
  for(l in 1:length(landmarks)) {
    temp <- unlist(lapply(predictions, function(x) ifelse(!is.null(nrow(x)),x[l,eval],NA)))
    roc.pred <- pROC::roc(response=true.response, predictor=temp, direction="<")
    pROC::plot.roc(roc.pred, main=paste("Landmark",landmarks[l],"y, disease risk",eval.time,"y"))
    legend("bottomright", legend=paste("AUC",round(roc.pred$auc,3)), bty="n")
    mtext(mtext)
  }
}

get.auc <- function(eval.time, landmark.time, predictions, true.response) {
  library(pROC)
  eval.time <- paste("eval", eval.time, sep=".")
  landmark.time <- paste("landmark", landmark.time, sep=".")
  pred <- unlist(lapply(predictions, function(x) ifelse(!is.null(nrow(x)),x[landmark.time,eval.time],NA)))
  c(t(auc(pROC::roc(response=true.response, predictor=pred, direction="<"))))
}

get.ci <- function(eval.time, landmark.time, predictions, true.response) {
  library(pROC)
  eval.time <- paste("eval", eval.time, sep=".")
  landmark.time <- paste("landmark", landmark.time, sep=".")
  pred <- unlist(lapply(predictions, function(x) ifelse(!is.null(nrow(x)),x[landmark.time,eval.time],NA)))
  c(t(ci(pROC::roc(response=true.response, predictor=pred, direction="<"))))
}

