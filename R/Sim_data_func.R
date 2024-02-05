
Sim_data_func <- function(Y, d, t, X, sim.size=2000){
  
  n= length(t)
  mydata <- data.frame(Y, d, X) # make dataframe for model fitting
  trt.ind <- as.numeric(t) # create treatment indicator variable
  
  gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model
  var <- paste(cox.variables(X), collapse = "+") ## cox variables for time to event models
  
  ## fit models
  t.fit = mgcv::gam(as.formula(paste("trt.ind ~", gam.var)), data=X, family=binomial) ## treatment model
  
  data_t1 <- mydata[trt.ind==1,] # get data for treated group
  
  S_t1 = survreg(as.formula(paste("Surv(time = Y, event = d) ~", var)), dist="weibull", data = data_t1) 
  
  G_t1 = survreg(as.formula(paste("Surv(time = Y, event = 1-d) ~", var)), dist="weibull", data = data_t1) 
  
  data_t0 <- mydata[trt.ind!=1,] # get data for control group
  
  S_t0 = survreg(as.formula(paste("Surv(time = Y, event = d) ~", var)), dist="weibull", data = data_t0) 
  
  G_t0 = survreg(as.formula(paste("Surv(time = Y, event = 1-d) ~", var)), dist="weibull", data = data_t0) 
  
  
  ## get random sample from data
  indices <- sample.int(n, size=sim.size, replace=T)
  X.new <- X[indices,]
  
  ## predict treatment assignment given the random sample
  prob.treat.assign <- predict (t.fit, newdata=X.new, type="response") ## prob of receiving treatment
  t.new <- rbinom(n=sim.size,size=1,prob=prob.treat.assign)
  
  
  X.new.t1 <- X.new[t.new==1,]
  
  time.fun <- function(p, survreg.scale, survreg.lp){
    shape <- 1/survreg.scale 
    scale <- exp(survreg.lp)
    ans <- qweibull(p, shape = shape, scale = scale, lower.tail = F)
  }
  
  lp.t1 <- predict(S_t1, newdata = X.new.t1, type = "lp")
  event.time.t1 <- sapply(1:nrow(X.new.t1), function(x){
    ecdf <- runif(1,0,1)
    time.fun(ecdf, survreg.scale = S_t1$scale, survreg.lp = lp.t1[x])
  })
  
  lp.t1.censor <- predict(G_t1, newdata = X.new.t1, type = "lp")
  censor.time.t1 <- sapply(1:nrow(X.new.t1), function(x){
    ecdf <- runif(1,0,1)
    time.fun(ecdf, survreg.scale = G_t1$scale, survreg.lp = lp.t1.censor[x])
  })
  
  X.new.t0 <- X.new[t.new==0,]
  
  lp.t0 <- predict(S_t0, newdata = X.new.t0, type = "lp")
  event.time.t0 <- sapply(1:nrow(X.new.t0), function(x){
    ecdf <- runif(1,0,1)
    time.fun(ecdf, survreg.scale = S_t0$scale, survreg.lp = lp.t0[x])
  })
  
  lp.t0.censor <- predict(G_t0, newdata = X.new.t0, type = "lp")
  censor.time.t0 <- sapply(1:nrow(X.new.t0), function(x){
    ecdf <- runif(1,0,1)
    time.fun(ecdf, survreg.scale = G_t0$scale, survreg.lp = lp.t0.censor[x])
  })
  
  
  Obs_data <- data.frame(rbind(cbind(event.time=event.time.t1,censoring.time=censor.time.t1,t=1,X.new.t1,indices=indices[t.new==1]),
                                cbind(event.time=event.time.t0,censoring.time=censor.time.t0,t=0,X.new.t0,indices=indices[t.new==0])))

  
  Obs_data$Y <- pmin(Obs_data$event.time,Obs_data$censoring.time)
  Obs_data$delta <- as.numeric(Obs_data$event.time<=Obs_data$censoring.time)
  
  Obs_data <- Obs_data[sample(1:nrow(Obs_data)), ] ## shuffle the data
  
  return(Obs_data)
  
  
}
