library(lcmm)

#--------------------------------------------------
# lcmm
#--------------------------------------------------

train.lcmm<- function(data.train) {
  set.seed(12345)
  lcmm.m1 <- Jointlcmm(intensity ~ time,
                       random = ~ time,
                       survival = Surv(event_time, status) ~ 1,
                       hazard = "Weibull", subject = "id", data = data.train, ng = 1)
  # gridsearch with 15 iterations from 30 random departures
  lcmm.m2 <- gridsearch(rep = 30, maxiter = 15, minit = lcmm.m1,
                        Jointlcmm(intensity ~ time,
                                  mixture = ~ time,
                                  random = ~ time,
                                  survival = Surv(event_time, status) ~ 1, 
                                  hazard = "Weibull", subject = "id", data = data.train, ng = 2, verbose = FALSE))
  lcmm.m <- lcmm.m1
  if(lcmm.m2$BIC < lcmm.m1$BIC) lcmm.m <- lcmm.m2
  return(lcmm.m)
}

train.lcmm.randint <- function(data.train) {
  set.seed(12345)
  lcmm.m1 <- Jointlcmm(intensity ~ time,
                       random = ~ 1,
                       survival = Surv(event_time, status) ~ 1,
                       hazard = "Weibull", subject = "id", data = data.train, ng = 1)
  # gridsearch with 15 iterations from 30 random departures
  lcmm.m2 <- gridsearch(rep = 30, maxiter = 15, minit = lcmm.m1,
                        Jointlcmm(intensity ~ time,
                                  mixture = ~ time,
                                  random = ~ 1,
                                  survival = Surv(event_time, status) ~ 1, 
                                  hazard = "Weibull", subject = "id", data = data.train, ng = 2, verbose = FALSE))
  lcmm.m<- lcmm.m1
  if(lcmm.m2$BIC < lcmm.m1$BIC) lcmm.m<- lcmm.m2
  return(lcmm.m)
}

test.lcmm <- function(data.test, model.lcmm, landmark, horizon) {
  set.seed(12345)
  dynp.case <- dynpred(model.lcmm, data.test, landmark = landmark, var.time = "time", horizon = horizon)
  result <- matrix(nrow=length(landmark), ncol=length(horizon))
  row.names(result) <- paste("landmark", landmark, sep=".")
  colnames(result) <- paste("eval", landmark[1] + horizon, sep=".")
  for(l in 1:length(landmark)) {
    temp<- dynp.case$pred[dynp.case$pred[,"landmark"]==landmark[l],,drop=FALSE]
    temp<- temp[temp[,"landmark"] + temp[,"horizon"] <= landmark[1] + max(horizon),,drop=FALSE]
    result[l,] <- c(rep(NA,length(horizon)-nrow(temp)), temp[,"pred"])
  }
  return(result)
}

pred.lcmm <-  function(model.lcmm, data.test, data.test.id, landmark, horizon) {
  result <- list()
  for(i in 1:nrow(data.test.id)) {
    test <- data.test[data.test[,"subject"]==data.test.id[i,"subject"],c("id","time","intensity")]
    result[[i]] <- tryCatch(test.lcmm(test, model.lcmm, landmark, horizon), error=function(e) NA)
  }
  names(result) <- data.test.id[,"subject"]
  return(result)
}
