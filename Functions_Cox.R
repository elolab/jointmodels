#devtools::install_github("mdbrown/partlyconditional")
library(survival)

#--------------------------------------------------
# Cox
#--------------------------------------------------

# Train
train.Cox <- function(data.train) {
  set.seed(12345)
  model.cox <-  coxph(Surv(time, status) ~ intensity, data=data.train)
  return(model.cox)
}

# Test
test.Cox <- function(data.test, model.Cox) {
  set.seed(12345)
  newd <- cbind(data.test, status=rep(0,nrow(data.test)))
  newd <- newd[1,,drop=FALSE]
  result <- 1 - predict(model.Cox, newdata = newd, type="survival")
  return(result)
}

# Predict
pred.Cox <-  function(model.Cox, data.test, data.test.id) {
  result <- list()
  for(i in 1:nrow(data.test.id)) {
    test <- data.test[data.test[,"subject"]==data.test.id[i,"subject"],c("id","time","intensity")]
    result[[i]] <- test.Cox(test, model.Cox)
  }
  names(result) <- data.test.id[,"subject"]
  return(result)
}

# Subset data, utility function
subsetData <- function(data, limit) {
  subset <- foreach(s=unique(data$subject), .combine=rbind) %do% {
    dat <- data[which(data$subject==s),]
    dat <- dat[which(dat$time<limit),]
    dat[which.max(dat$time),]
  }
  return(subset)
}

