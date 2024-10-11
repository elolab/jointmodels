source("Functions_plotting.R")
library(foreach)
library(reshape2)

simulate.tp <- function(n.tp, tp.min, tp.max) {
  sapply(0:(n.tp-1), function(i) {
    min <- tp.min + i*((tp.max-tp.min)/n.tp)
    max <- tp.min + (i+1)*((tp.max-tp.min)/n.tp)
    runif(1, min=min, max=max)
  })
}

simulate.surv <- function(groupnames=c("case","control"), n.subjects=5, n.tp=5, tp.min=0, tp.max=10, event1.mean=15, event1.sd=1, event2.mean=15, event2.sd=1, event1.fix=20) {
  sim <- foreach(i=1:n.subjects, .combine=rbind) %do% {
    a <- data.frame(group=groupnames[1],
                    subject=paste(groupnames[1],i,sep="_"),
                    time=simulate.tp(n.tp, tp.min, tp.max),
                    event_time=rnorm(n=1, mean=event1.mean, sd=event1.sd))
    b <- data.frame(group=groupnames[2],
                    subject=paste(groupnames[2],i,sep="_"),
                    time=simulate.tp(n.tp, tp.min, tp.max),
                    event_time=rnorm(n=1, mean=event2.mean, sd=event2.sd))
    rbind(a,b)
  }
  
  sim$id <- as.numeric(as.factor(sim[,"subject"]))
  sim$status <- ifelse(sim$group=="case",1,0)
  
  sel <- which(sim$status==1)
  if(max(sim$event_time[sel])<event1.fix) {
    sel <- sel[which(sim$event_time[sel] == max(sim$event_time[sel]))]
    sim$event_time[sel] <- event1.fix
  }
  
  return(sim)
}

simulate.intensity <- function(survdata, groupnames=c("case","control"), mean=0, sd=0, sd.subject=0, diff=0, trend1=0, trend2=0) {
  survdata$intensity <- rnorm(n=nrow(survdata), mean=mean, sd=sd)
  survdata$intensity <-  survdata$intensity + ifelse(survdata$group==groupnames[1],diff,-diff)
  survdata$intensity <-  survdata$intensity + survdata$time * ifelse(survdata$group==groupnames[1],trend1,-trend1)
  survdata$intensity <-  survdata$intensity + survdata$time^2 * ifelse(survdata$group==groupnames[1],trend2,-trend2)
  subjects <- unique(survdata$subject)
  sd.subject <- rnorm(n=length(subjects), mean=0, sd=sd.subject)
  for(s in 1:length(subjects)) {
    sel <- which(survdata$subject==subjects[s])
    survdata$intensity[sel] <- survdata$intensity[sel] + sd.subject[s]
  }
  neg <- which(survdata$intensity<0)
  if (length(neg)>0) {
    survdata$intensity[neg] <- 0
    if(length(neg)>(0.01*length(survdata$intensity))) {
      warning("More than 1% of values were replaced with 0")
    }
  }
  return(survdata)
}

# Run simulations
trend1 <- seq(0, 0.1, length.out=5)
var <- seq(0.5, 1.5, length.out=3)
subj <- c(10,50)

for(s in subj) {
  for(i in 1:length(trend1)) {
    for(j in 1:length(var)) {
      set.seed(12345)
      # Training data
      survdata.train <- simulate.surv(n.subjects=s, n.tp=10, tp.min=0, tp.max=15, event1.mean=18, event1.sd=1, event2.mean=18, event2.sd=1, event1.fix=20.1)
      data.train <- simulate.intensity(survdata.train, mean=6, sd=var[j], sd.subject=0.1, diff=0, trend1=trend1[i], trend2=0)
      # Test data
      survdata.test <- simulate.surv(n.subjects=s, n.tp=10, tp.min=0, tp.max=15, event1.mean=18, event1.sd=1, event2.mean=18, event2.sd=1, event1.fix=20.1)
      data.test <- simulate.intensity(survdata.test, mean=6, sd=var[j], sd.subject=0.1, diff=0, trend1=trend1[i], trend2=0)
      # Save
      filename <- paste("subj",s,"_trend",i,"_var",j,sep="")
      save(data.train, data.test, file=paste("../Data/Simulation/simulation_",filename,".RData",sep=""))
    }
  }
}

# Plot the generated data
pdf("../Results/Simulation/Simulation_trend.pdf", width=10, height=15)
par(mfrow=c(5,3))
for(pattern in c("_subj50","_subj10_")) {
  files <- list.files("../Data/Simulation", pattern=pattern)
  for(f in files) {
    load(paste("../Data/Simulation/",f,sep=""))
    plotProfile(data.test, individuals=TRUE, ylim=c(2,10))
    mtext(f)
  }
}
dev.off()

# Plot single examples
files <- list.files("../Data/Simulation")
load("../Data/Simulation/simulation_subj50_trend3_var3.RData")
pdf("../Results/Simulation/Simulation_subj50_trend3_var3.pdf", width=3, height=3)
plotProfile(data.test, individuals=TRUE, ylim=c(2,10))
dev.off()

load("../Data/Simulation/simulation_subj50_trend3_var2.RData")
pdf("../Results/Simulation/Simulation_subj50_trend3_var2.pdf", width=3, height=3)
plotProfile(data.test, individuals=TRUE, ylim=c(2,10))
dev.off()

load("../Data/Simulation/simulation_subj50_trend3_var1.RData")
pdf("../Results/Simulation/Simulation_subj50_trend3_var1.pdf", width=3, height=3)
plotProfile(data.test, individuals=TRUE, ylim=c(2,10))
dev.off()

load("../Data/Simulation/simulation_subj50_trend5_var1.RData")
pdf("../Results/Simulation/Simulation_subj50_trend5_var1.pdf", width=3, height=3)
plotProfile(data.test, individuals=FALSE, ylim=c(2,10))
dev.off()


