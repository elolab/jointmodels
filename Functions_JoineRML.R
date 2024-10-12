library(joineRML)

#--------------------------------------------------
# JoineRML
#--------------------------------------------------

# Train
train.joineRML<- function(data.train) {
  set.seed(12345)
  fit.joineRML <- mjoint(formLongFixed = list(intensity ~ time),
                         formLongRandom = list(~ time | id), 
                         formSurv = Surv(event_time, status) ~ 1,
                         data = data.train,
                         timeVar = "time")     
  return(fit.joineRML)
}

# Train with random intercept
train.joineRML.randint <- function(data.train) {
  set.seed(12345)
  fit.joineRML <- mjoint(formLongFixed = list(intensity ~ time),
                         formLongRandom = list(~ 1 | id), 
                         formSurv = Surv(event_time, status) ~ 1,
                         data = data.train,
                         timeVar = "time")     
  return(fit.joineRML)
}

# Test
test.joineRML <- function(data.test, model.joineRML, landmark, horizon) {
  set.seed(12345)
  result <- matrix(nrow=length(landmark), ncol=length(horizon))
  row.names(result) <- paste("landmark", landmark, sep=".")
  colnames(result) <- paste("eval", landmark[1] + horizon, sep=".")
  max.time <- max(model.joineRML$survData[model.joineRML$survData[,"status"]==1,"event_time"])
  for(l in 1:length(landmark)) {
    newd <- data.test[data.test[,"time"]<= landmark[l],,drop=FALSE]
    sel <- which(horizon + landmark[l] < max.time)
    nas <- length(horizon) - length(sel)
    temp <- dynSurv(model.joineRML, newdata = newd, horizon = horizon[sel] + (landmark[l]-newd[nrow(newd),"time"]))
    temp <- temp$pred
    temp <- temp[temp[,"u"] <= landmark[1] + max(horizon),,drop=FALSE]
    result[l,] <- 1 - c(rep(NA,landmark[l] - landmark[1]), temp[,"surv"], rep(NA,length(horizon) - (landmark[l] - landmark[1]) - nrow(temp)))
  }
  return(result)
}

# Predict
pred.joineRML <- function(model.joineRML, data.test, data.test.id, landmark, horizon) {
  result <- list()
  for(i in 1:nrow(data.test.id)) {
    test <- data.test[data.test[,"subject"]==data.test.id[i,"subject"],c("id","time","intensity")]
    result[[i]] <- tryCatch(test.joineRML(test, model.joineRML, landmark, horizon), error=function(e) NA)
  }
  names(result) <- data.test.id[,"subject"]
  return(result)
}
