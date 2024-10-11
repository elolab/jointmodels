#devtools::install_github("mdbrown/partlyconditional")
library(partlyconditional)

#--------------------------------------------------
# PCcox
#--------------------------------------------------

train.pcCox <- function(data.train) {
  set.seed(12345)
  pc.cox <-  PC.Cox(id = "id", stime = "event_time", status = "status", measurement.time = "time", predictors =c("time", "intensity"), data = data.train)
  return(pc.cox)
}

test.pcCox<- function(data.test, model.pcCox, landmark, horizon, max.time) {
  set.seed(12345)
  
  var.missing <- setdiff(model.pcCox$variable.names, colnames(data.test))
  data.test.append <- data.test
  for(i in 1:length(var.missing)) data.test.append <- cbind(data.test.append, rep(max(landmark)+max(horizon),nrow(data.test)))
  colnames(data.test.append) <- c(colnames(data.test), var.missing)
  
  result <- matrix(nrow=length(landmark), ncol=length(horizon))
  row.names(result) <- paste("landmark", landmark, sep=".")
  colnames(result) <- paste("eval", landmark[1] + horizon, sep=".")
  for(l in 1:length(landmark)) {
    newd <- data.test.append[data.test.append[,"time"]<= landmark[l],,drop=FALSE]
    sel <- which(horizon + landmark[l] < max.time)
    temp <- predict(model.pcCox, newdata = newd, prediction.time = horizon[sel])
    temp <- as.numeric(as.character(temp[nrow(temp), (ncol(newd)+2):(ncol(newd)+2+length(horizon[sel])-1)]))
    result[l,] <- c(rep(NA,landmark[l] - landmark[1]), temp, rep(NA,length(horizon) - (landmark[l] - landmark[1]) - length(temp)))
  }
  return(result)
  return(result)
}

pred.pcCox <-  function(model.pcCox, data.test, data.test.id, landmark, horizon, max.time) {
  result <- list()
  for(i in 1:nrow(data.test.id)) {
    test <- data.test[data.test[,"subject"]==data.test.id[i,"subject"],c("id","time","intensity")]
    result[[i]]<- tryCatch(test.pcCox(test, model.pcCox, landmark, horizon, max.time), error=function(e) NA)
  }
  names(result)<- data.test.id[,"subject"]
  return(result)
}
