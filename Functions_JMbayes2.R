library(JMbayes2)

#--------------------------------------------------
# JMbayes2
#--------------------------------------------------

train.JMbayes2 <- function(data.train, data.train.id) {
  set.seed(12345)
  fitLME <- lme(intensity ~ time, random = ~ time | id, data = data.train)
  fitSURV <- coxph(Surv(event_time, status) ~ 1, data = data.train.id, x = TRUE, model = TRUE)
  fit.JMbayes2 <- jm(fitSURV, fitLME, time_var = "time", n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
  return(fit.JMbayes2)
}

test.JMbayes2 <- function(data.test, model.JMbayes2, landmark, horizon) {
  set.seed(12345)
  result<- matrix(nrow=length(landmark), ncol=length(horizon))
  row.names(result) <- paste("landmark", landmark, sep=".")
  colnames(result) <- paste("eval", landmark[1] + horizon, sep=".")
  for(l in 1:length(landmark)) {
    newd <- data.test[data.test[,"time"]<= landmark[l],,drop=FALSE]
    temp <- predict(model.JMbayes2, newdata = newd, process = "event", times = landmark[l]+horizon, return_newdata = TRUE)
    temp <- temp[temp[,"time"] <= landmark[1] + max(horizon),,drop=FALSE]
    temp <- temp[temp[,"time"] >= landmark[l] + min(horizon),,drop=FALSE]
    result[l,] <- c(rep(NA,landmark[l] - landmark[1]), temp[,"pred_CIF"], rep(NA,length(horizon) - (landmark[l] - landmark[1]) - nrow(temp)))
  }
  return(result)
}

pred.JMbayes2 <- function(model.JMbayes2, data.test, data.test.id, landmark, horizon) {
  result <- list()
  for(i in 1:nrow(data.test.id)) {
    test <- data.test[data.test[,"subject"]==data.test.id[i,"subject"],c("id","time","intensity")]
    result[[i]] <- tryCatch(test.JMbayes2(test, model.JMbayes2, landmark, horizon), error=function(e) NA)
  }
  names(result) <- data.test.id[,"subject"]
  return(result)
}
