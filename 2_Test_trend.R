library(plyr)
source("Functions_plotting.R")
source("Functions_JoineRML.R")
source("Functions_PCcox.R")
source("Functions_lcmm.R")
source("Functions_JM.R")
source("Functions_Cox.R")

# Landmarks and horizons to test
landmark <- c(5,10,15)
horizon <- c(1:15)

# Run for all files
files <- list.files("../Data/Simulation")
for(f in files) {
  
  # Load simulated data
  load(paste("../Data/Simulation/",f,sep=""))
  
  # Data including only the first time point from each individual.
  data.train.id <- ddply(data.train, .(group, subject, id), numcolwise(mean))
  data.test.id <- ddply(data.test, .(group, subject, id), numcolwise(mean))
  
  # JoineRML
  results.joineRML <- list()
  tryCatch({
    results.joineRML$time <- c(t(system.time(model.joineRML <- train.joineRML(data.train))[3]))
    predictions.joineRML <- pred.joineRML(model.joineRML, data.test, data.test.id, landmark, horizon)
    results.joineRML$AIC <- summary(model.joineRML)$AIC
    results.joineRML$BIC <- summary(model.joineRML)$BIC
    results.joineRML$p <- summary(model.joineRML)$coefs.surv[4]
    results.joineRML$AUC <- sapply(landmark, function(x) get.auc(eval.time=20, landmark.time=x, predictions=predictions.joineRML, true.response=data.test.id[,"status"]))
    names(results.joineRML$AUC) <- landmark
    results.joineRML$CI <- sapply(landmark, function(x) get.ci(eval.time=20, landmark.time=x, predictions=predictions.joineRML, true.response=data.test.id[,"status"]))
    colnames(results.joineRML$CI) <- landmark
  }, error=function(e) {})
  
  # lcmm
  results.lcmm <- list()
  tryCatch({
    results.lcmm$time <- c(t(system.time(model.lcmm <- train.lcmm(data.train))[3]))
    predictions.lcmm <- pred.lcmm(model.lcmm, data.test, data.test.id, landmark, horizon)
    results.lcmm$BIC <- summarytable(model.lcmm)[,4]
    results.lcmm$summary <- summary(model.lcmm)
    results.lcmm$AUC <- sapply(landmark, function(x) get.auc(eval.time=20, landmark.time=x, predictions=predictions.lcmm, true.response=data.test.id[,"status"]))
    names(results.lcmm$AUC) <- landmark
    results.lcmm$CI <- sapply(landmark, function(x) get.ci(eval.time=20, landmark.time=x, predictions=predictions.lcmm, true.response=data.test.id[,"status"]))
    colnames(results.lcmm$CI) <- landmark
  }, error=function(e) {})
  
  # JM
  results.JM <- list()
  tryCatch({
    results.JM$time <- c(t(system.time(model.JM <- train.JM(data.train, data.train.id))[3]))
    predictions.JM <- pred.JM(model.JM, data.test, data.test.id, landmark, horizon)
    results.JM$AIC <- summary(model.JM)$AIC
    results.JM$BIC <- summary(model.JM)$BIC
    results.JM$p <- summary(model.JM)$"CoefTable-Event"["Assoct","p-value"]
    results.JM$summary <- summary(model.JM)
    results.JM$AUC <- sapply(landmark, function(x) get.auc(eval.time=20, landmark.time=x, predictions=predictions.JM, true.response=data.test.id[,"status"]))
    names(results.JM$AUC) <- landmark
    results.JM$CI <- sapply(landmark, function(x) get.ci(eval.time=20, landmark.time=x, predictions=predictions.JM, true.response=data.test.id[,"status"]))
    colnames(results.JM$CI) <- landmark
  }, error=function(e) {})
  
  # pcCox
  results.pcCox <- list()
  tryCatch({
    results.pcCox$time <- c(t(system.time(model.pcCox <- train.pcCox(data.train))[3]))
    predictions.pcCox <- pred.pcCox(model.pcCox, data.test, data.test.id, landmark, horizon, max.time=max(data.train$event_time[data.train$status==1]))
    results.pcCox$p <- summary(model.pcCox$model.fit)$coefficients["intensity","Pr(>|z|)"]
    results.pcCox$AUC <- sapply(landmark, function(x) get.auc(eval.time=20, landmark.time=x, predictions=predictions.pcCox, true.response=data.test.id[,"status"]))
    names(results.pcCox$AUC) <- landmark
    results.pcCox$CI <- sapply(landmark, function(x) get.ci(eval.time=20, landmark.time=x, predictions=predictions.pcCox, true.response=data.test.id[,"status"]))
    colnames(results.pcCox$CI) <- landmark
  }, error=function(e) {})
  
  # Baseline Cox
  results.Cox <- list()
  tryCatch({
    results.Cox$time <- c(t(system.time(model.Cox <- train.Cox(data.train.id))[3]))
    results.Cox$p <- summary(model.Cox)$coefficients[5]
    results.Cox$AUC <- foreach(limit=landmark, .combine=c) %do% {
      model.Cox <- train.Cox(subsetData(data.train, limit=limit))
      predictions.Cox <- pred.Cox(model.Cox, subsetData(data.test,limit=limit), data.test.id)
      c(t(auc(pROC::roc(response=data.test.id[,"status"], predictor=unlist(predictions.Cox), direction="<"))))
    }
    results.Cox$CI <- foreach(limit=landmark, .combine=cbind) %do% {
      model.Cox <- train.Cox(subsetData(data.train, limit=limit))
      predictions.Cox <- pred.Cox(model.Cox, subsetData(data.test,limit=limit), data.test.id)
      c(t(ci(pROC::roc(response=data.test.id[,"status"], predictor=unlist(predictions.Cox), direction="<"))))
    }
    colnames(results.Cox$CI) <- landmark
  }, error=function(e) {})
  
  # Save
  save(results.joineRML, results.lcmm, results.JM, results.pcCox, results.Cox, file=paste("../Results/Simulation/",f,sep=""))
  
}

