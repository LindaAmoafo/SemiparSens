F_s_hat <- function(Y, d, t, X, X.new, trt=1,s, gamma=0){
  
  d <- ifelse(Y>=150,1,d)
  Y <- pmin(Y,150)
  
  n= length(t)
  mydata <- data.frame(Y, d, X) # make dataframe for model fitting
  data_t <- mydata[t==trt,] # get data for treated group
  data_t0 <- mydata[t!=trt,] # get data for control group
  trt.ind <- as.numeric(t==trt) # create treatment indicator variable
  dt <- d[t==trt]
  dt0 <- d[t!=trt]
  
  if(is.null(s)){
    s <- unique(sort(Y[trt.ind==1 & d==1], decreasing =FALSE))
  } else{s=s} ## Use entire treatment death times in treatment group if s is not given.
  
  gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model
  var <- paste(cox.variables(X), collapse = "+") ## cox variables for time to event models
  
  
  ## fit models
  t.fit = mgcv::gam(as.formula(paste("trt.ind ~", gam.var)), data=X, family=binomial) ## treatment model
  
  S_t = coxph( as.formula(paste("Surv(time = Y, event = d) ~", var)),  ties="breslow",
               data = data_t) ## fit cox.ph model for treatment group with death status
  
  G_t = coxph( as.formula(paste("Surv(time = Y, event = 1-d) ~", var)),
               ties="breslow",data = data_t) ## fit cox.ph model for treatment group with censored status
  
  G_t0 = coxph( as.formula(paste("Surv(time = Y, event = 1-d) ~", var)),
                ties="breslow",data = data_t0) ## fit cox.ph model for control group with censored status
  
  
  ## get predictions for pi and survival probabilities
  pi <- predict (t.fit, newdata=X.new[trt.ind==1,], type="response") ## prob of receiving treatment for treated group only
  
  S_t_pred.obj <- survfit(S_t, stype=1, newdata=X.new, se.fit = FALSE) ## get survfit object for prediction with death for all covariates in the dataset
  pred.Yt <- Y[trt.ind==1]
  S_t_pred.sum <- summary(S_t_pred.obj, times=Y[trt.ind==1], extend = TRUE)  ## prediction object to get survival probabilities at all times in treatment time points
  S_t_s_mat<- summary(S_t_pred.obj, times=s, extend = TRUE)$surv ## get survival probabilities for each s to be used later in the loop. (sxn matrix)
  
  ## Working of probabilities at each person's lagged time with prediction with all covariates
  sort.Yt <- sort(Y[trt.ind==1], decreasing =FALSE) ## sort all the treated time points
  S_t_pred.sum.Yt <- S_t_pred.sum$surv[,trt.ind==1] ## get survival probabilities at treated covariates (nt x nt matrix)
  my.S_t_surv <- rbind(1, S_t_pred.sum.Yt[-nrow(S_t_pred.sum.Yt),]) ## lag the survival probabilities
  S_t_y_mat=my.S_t_surv[match(pred.Yt,sort.Yt),] ## match the lagged survival probability for each person
  F_t_y_mat=1-S_t_y_mat
  S_t_y_mat_dt = S_t_y_mat[dt==0,] ## filter S_t_y_mat to only censored observations
  F_t_y_mat_dt=1-S_t_y_mat_dt
  S_t_y=diag(S_t_y_mat) ## get diagonals to represent each persons lagged survival probability at their covariates
  F_t_y=1-S_t_y ## get failure probabilities
  
  
  G_t_pred.obj <- survfit(G_t, stype=1,  newdata=X.new[trt.ind==1,],se.fit = FALSE) ## get survfit object for prediction with censored model
  G_t_pred.sum <- summary(G_t_pred.obj, times=pred.Yt, extend = TRUE) $ surv ## get censored probabilities
  my.G_t_surv <- rbind(1, G_t_pred.sum[- nrow(G_t_pred.sum),]) ## lag the censored probabilities
  G_t_surv_mat=my.G_t_surv[match(pred.Yt,sort.Yt),] ## match the lagged censored probability for each person
  G_t_weight_mat <- ifelse((1/G_t_surv_mat)>100, 100, 1/G_t_surv_mat) ## truncate at weight =100
  G_t_weight_mat_dt <- G_t_weight_mat[dt==0,]
  G_t_weight=diag(G_t_weight_mat) ## get diagonals to represent each persons censored probability at their covariates
  
  beta_t <- G_t$coefficients
  ind.X <- model.matrix(as.formula(paste("~", var)), data = X.new)[trt.ind==1,-1]
  
  if(anyNA(beta_t)){
    no_coeff <- which(is.na(beta_t))
    beta_t <- beta_t[-no_coeff]
    ind.X <- ind.X[,-no_coeff]
  }
  
  ## Do same thing for treatment=0
  pred.Yt0 <- Y[trt.ind==0]
  sort.Yt0 <- sort(pred.Yt0, decreasing =FALSE)
  G_t0_pred.obj <- survfit(G_t0, stype=1,  newdata=X.new[trt.ind==0,],se.fit = FALSE) ## get survfit object for prediction with censored model
  G_t0_pred.sum <- summary(G_t0_pred.obj, times=pred.Yt0, extend = TRUE) $ surv ## get censored probabilities
  my.G_t0_surv <- rbind(1, G_t0_pred.sum[- nrow(G_t0_pred.sum),]) ## lag the censored probabilities
  G_t0_surv_mat=my.G_t0_surv[match(pred.Yt0,sort.Yt0),] ## match the lagged censored probability for each person
  G_t0_weight_mat <- ifelse((1/G_t0_surv_mat)>100, 100, 1/G_t0_surv_mat) ## truncate at weight =100
  G_t0_weight_mat_dt <- G_t0_weight_mat[dt0==0,]
  G_t0_weight=diag(G_t0_weight_mat) ## get diagonals to represent each persons censored probability at their covariates
  
  beta_t0 <- G_t0$coefficients
  ind.X_t0 <- model.matrix(as.formula(paste("~", var)), data = X.new)[trt.ind==0,-1]
  
  if(anyNA(beta_t0)){
    no_coeff_t0 <- which(is.na(beta_t0))
    beta_t0 <- beta_t0[-no_coeff_t0]
    ind.X_t0 <- ind.X_t0[,-no_coeff_t0]
  }
  
  Y.i.j_ind <- sapply(pred.Yt, function(x) as.numeric(pred.Yt[dt==0]<=x))
  Yt0.i.j_ind <- sapply(pred.Yt0, function(x) as.numeric(pred.Yt0[dt0==0]<=x))
  
  denn <- sapply(pred.Yt[dt==0],function(x) sum(c(exp(ind.X  %*% beta_t)) * as.numeric(pred.Yt>=x))) ## denominator with trt.ind==1
  
  denn_t0 <- sapply(pred.Yt0[dt0==0],function(x) sum(c(exp(ind.X_t0  %*% beta_t0)) * as.numeric(pred.Yt0>=x)))## denominator with trt.ind==0
  
  big.den <- (1/n) * sum(c(G_t_weight,G_t0_weight) - 
                           c(rowSums(t(G_t_weight_mat_dt * Y.i.j_ind / denn) * c(exp(ind.X  %*% beta_t))),
                             rowSums(t(G_t0_weight_mat_dt * Yt0.i.j_ind / denn_t0) * c(exp(ind.X_t0  %*%beta_t0)))))
  
  ## objects to store outputs
  
  Est <- vector(length=length(s))
  var_error<- vector(length=length(s))
  var_error.second.inc<- vector(length=length(s))
  
  
 
   for(i in 1:length(s)){
    S_t_s <- S_t_s_mat[i,trt.ind==1]
    F_t_s <- 1-S_t_s
    S_t0_s <- S_t_s_mat[i,trt.ind==0]
    F_t0_s <- 1-S_t0_s
    s.ind <- as.numeric(pred.Yt<=s[i])
    
    ## working phi stuff
    temp <- (1-pi)/pi * (exp(gamma)/ (S_t_s + exp(gamma)*F_t_s)^2)
    phi_t1 <- s.ind * (1+temp)     
    phi_t2 <- -1 * (F_t_s*temp) ## second term
    
    phi_t3 <- ( (F_t0_s*exp(gamma)) / (S_t0_s + exp(gamma)*F_t0_s) ) ## third part as in formula
    
    ## Finishing it all up with the weight
    T1.phi_t1 <- G_t_weight*dt*phi_t1
    T1.phi_t2 <- G_t_weight*dt*phi_t2
    T1.phi_t3 <- G_t0_weight*dt0*phi_t3    
    
    ## Working on number 2 in the fomula
    l_t1 <- s.ind * ((F_t_s - F_t_y) / S_t_y) * (1+temp)
    
    T2.l_t1 <- G_t_weight*(1-dt)* l_t1
    T2.l_t2 <- G_t_weight*(1-dt)* phi_t2
    T2.l_t3 <- G_t0_weight*(1-dt0)*phi_t3
    
    
    ## Working on number 3 in the fomula
    T3.l_t1 <- colSums( ( G_t_weight_mat_dt * as.numeric(pred.Yt[dt==0] <= s[i]) * 
                            t(((F_t_s - t(F_t_y_mat_dt))/t(S_t_y_mat_dt)) *
                                c((1+temp) *c(exp(ind.X  %*% beta_t))) ) * Y.i.j_ind ) / denn, na.rm=T)     
    
    T3.l_t2 <- phi_t2 * rowSums(t(G_t_weight_mat_dt * Y.i.j_ind / denn) * c(exp(ind.X  %*% beta_t)))
    
    
    T3.l_t3 <- phi_t3 * rowSums(t(G_t0_weight_mat_dt * Yt0.i.j_ind / denn_t0) * c(exp(ind.X_t0  %*% beta_t0)))
    
    
    Est.temp <- ((1/n)*sum(c(T1.phi_t1 + T1.phi_t2 + T2.l_t1 + T2.l_t2 - T3.l_t1 - T3.l_t2, T1.phi_t3 + T2.l_t3 - T3.l_t3)))/big.den
    
    
    var_error.first <- ( (1/n) *sum(c(G_t_weight,G_t0_weight) - 
                                        c(rowSums(t(G_t_weight_mat_dt * Y.i.j_ind / denn) * c(exp(ind.X%*%beta_t))),
                                          rowSums(t(G_t0_weight_mat_dt * Yt0.i.j_ind / denn_t0) * 
                                                    c(exp(ind.X_t0 %*% beta_t0))))) )^(-2)
    
    var_error.second <- (1/n^2)*sum((c(T1.phi_t1 + T1.phi_t2 + T2.l_t1 + T2.l_t2 - T3.l_t1 - T3.l_t2, 
                                     T1.phi_t3 + T2.l_t3 - T3.l_t3) - 
                                     (c(G_t_weight,G_t0_weight) * Est.temp - 
                                        c(rowSums(t(G_t_weight_mat_dt * Y.i.j_ind / denn) * c(exp(ind.X%*%beta_t))),
                                          rowSums(t(G_t0_weight_mat_dt * Yt0.i.j_ind / denn_t0) * 
                                                    c(exp(ind.X_t0 %*% beta_t0))))* Est.temp))^2)
    
    
    
    var_error.temp <- var_error.first * var_error.second
    Est[i] <- Est.temp
    var_error[i] <- var_error.temp
    var_error.second.inc[i] <- var_error.second
    
  }
  
  No.neg.pava.estimate <- ifelse(Est<0, 0, Est)
  All.G.pava.estimate <- pava(No.neg.pava.estimate ,decreasing=FALSE,long.out = FALSE)
  Lower.Limit <- All.G.pava.estimate - qnorm(0.975)*sqrt(var_error)
  Upper.Limit <- All.G.pava.estimate + qnorm(0.975)*sqrt(var_error)
  
  # Logit transformation
  Trans_Est <- log(All.G.pava.estimate/(1-All.G.pava.estimate))
  MOE <- qnorm(0.975) * sqrt((1/(All.G.pava.estimate*(1-All.G.pava.estimate)))^2 * var_error)
  Lower.Limit.logit <- exp(Trans_Est - MOE)/(1+exp(Trans_Est - MOE))
  Upper.Limit.logit <- exp(Trans_Est + MOE)/(1+exp(Trans_Est + MOE))
  
  return(list(big.den=big.den,
              Est=Est,
              pava.Est=All.G.pava.estimate,
              var_error=var_error,
              var_error.second.inc=var_error.second.inc,
              Trans_Est =Trans_Est,
              MOE=MOE,
              Lower.Limit.logit=Lower.Limit.logit,
              Upper.Limit.logit=Upper.Limit.logit))
  
  
}

## Induced counterfactual Curves

counterfactual <- function(Y,d,t,X,trt=0,s=NULL, gamma){
  d <- ifelse(Y>=150,1,d)
  Y <- pmin(Y,150)
  
  n= length(t)
  mydata <- data.frame(Y, d, X) # make dataframe for model fitting
  data_t <- mydata[t==trt,] # get data for treated group
  trt.ind <- as.numeric(t==trt) # create treatment indicator variable
  
  gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model
  var <- paste(cox.variables(X), collapse = "+") ## cox variables for time to event models
  
  if(is.null(s)){
    s <- unique(sort(Y[trt.ind==1 & d==1], decreasing =FALSE))
  } else{s=s} ## Use entire treatment death times in treatment group if s is not given.
  
  ## fit models
  t.fit = mgcv::gam(as.formula(paste("trt.ind ~", gam.var)), data=X, family=binomial) ## treatment model
  S_t = coxph( as.formula(paste("Surv(time = Y, event = d) ~", var)),  ties="breslow",
               data = data_t) ## fit cox.ph model for treatment group with death status
  
  ## get predictions for pi and survival probabilities
  pi <- predict (t.fit, newdata=X, type="response") ## prob of receiving treatment
  
  S_t_pred.obj <- survfit(S_t, stype=1, newdata=X, se.fit = FALSE) ## get survfit object for prediction with death 
  S_t_s <- summary(S_t_pred.obj, times=s, extend = TRUE)$surv  ## prediction object to get survival probabilities
  F_t_s <- 1-S_t_s
  
  true_est <- (1/n) * rowSums(t(t(F_t_s)*c(pi)) + t(t((F_t_s*exp(gamma))/(S_t_s + (F_t_s*exp(gamma)))) * c(1-pi)))
  
  prob_fit <- rowMeans(S_t_s)
  prop_t <- sum(trt.ind)/n
  
  l <- cbind(counter=1-((true_est - ((1-prob_fit) * prop_t))/(1-prop_t)), prob_fit=prob_fit)
  return(l)
  
}




