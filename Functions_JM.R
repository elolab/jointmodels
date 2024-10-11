library(JM)

#--------------------------------------------------
# JM
#--------------------------------------------------

train.JM <- function(data.train, data.train.id) {
  set.seed(12345)
  fitLME <- lme(intensity ~ time, random = ~ time | id, data = data.train)
  fitSURV <- coxph(Surv(event_time, status) ~ 1, data = data.train.id, x = TRUE, model = TRUE)
  fit.JM <- jointModel(fitLME, fitSURV, timeVar="time", method="weibull-PH-aGH")
  return(fit.JM)
}

test.JM <- function(data.test, model.JM, landmark, horizon) {
  set.seed(12345)
  result <- matrix(nrow=length(landmark), ncol=length(horizon))
  row.names(result) <- paste("landmark", landmark, sep=".")
  colnames(result) <- paste("eval", landmark[1] + horizon, sep=".")
  for(l in 1:length(landmark)) {
    newd <- data.test[data.test[,"time"]<= landmark[l],,drop=FALSE]
    temp <- survfitJM(model.JM, newdata=newd, idVar="id", survTimes=landmark[l]+horizon)
    temp <- temp$summaries[[1]]
    temp <- temp[temp[,"times"] <= landmark[1] + max(horizon),,drop=FALSE]
    result[l,] <- 1 - c(rep(NA,length(horizon)-nrow(temp)), temp[,"Median"])
  }
  return(result)
}

pred.JM <-  function(model.JM, data.test, data.test.id, landmark, horizon) {
  result <- list()
  for(i in 1:nrow(data.test.id)) {
    test <- data.test[data.test[,"subject"]==data.test.id[i,"subject"],c("id","time","intensity")]
    result[[i]] <- tryCatch(test.JM(test, model.JM, landmark, horizon), error=function(e) NA)
  }
  names(result) <- data.test.id[,"subject"]
  return(result)
}
