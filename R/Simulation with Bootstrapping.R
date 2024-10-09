# Set my usual working directory. 
# Change according to who is testing!!
#setwd("/uufs/chpc.utah.edu/common/home/u1271266/CHPC/Multiple/divide_and_conquer/Stratify_Analysis/New/EBRT/BootstrapCI/Gamma_0")

library(survival)
library(RcppZiggurat)
library(dplyr)
library(rms)
library(data.table)
library(Rfast)
library(Rcpp)
library(collapse)
library(splines)
#library(table1)
library(parallel)
library(Iso)
library(rlist)
library(doSNOW)

## Code for one simulation setup to be parallelized on the high performance computer. 


epsilom <- 0.00000001
data_9_11_22 <- read.csv("data_9_11_22.csv") 
# clinical T-stage is TNM_CLIN_T
# PSA is CS_SITESPECIFIC_FACTOR_1
# GS is CS_SITESPECIFIC_FACTOR_8
# classify as high risk on the basis of either clinical T stage T3 or greater, or biopsy Gleason score of 8 to 10, or pretreatment PSA. 20 ng/dL.

High_Risk_data <- data_9_11_22 %>% filter (TNM_CLIN_T %in% c("c3", "c3A", "c3B", "c3C", "c4") | 
                                             CS_SITESPECIFIC_FACTOR_8  %in% c(8,9,10) |
                                             (CS_SITESPECIFIC_FACTOR_1 >= 200 &  CS_SITESPECIFIC_FACTOR_1<=980))

High_Risk_dat <- High_Risk_data %>%
  mutate(days2trx= ifelse(days2trx==0, .5, days2trx)) 

# Recategorise some of the data
dat_High <- High_Risk_dat %>% mutate(
  CS_SITESPECIFIC_FACTOR_8= ifelse(CS_SITESPECIFIC_FACTOR_8<=6, 6, CS_SITESPECIFIC_FACTOR_8), ## GS less than 6 set to 6
  RACE = ifelse(!(RACE %in% c(1,2)), 3, RACE), #combine non-white & non-black race group
  SPANISH_HISPANIC_ORIGIN = ifelse(SPANISH_HISPANIC_ORIGIN %in% c(1:8), 1, SPANISH_HISPANIC_ORIGIN),## combine hispanic groups
  #tstage =ifelse(TNM_CLIN_T == 88, NA, ifelse(TNM_CLIN_T %in% c("c3", "c3A", "c3B", "c3C", "c4"),2,1)),
  tstage=ifelse(TNM_CLIN_T == 88, NA,TNM_CLIN_T),
  CDCC_TOTAL_BEST = ifelse(CDCC_TOTAL_BEST ==0,"0",">=1")) %>% 
  
  select(c("AGE", "RACE" , "SPANISH_HISPANIC_ORIGIN",  "INSURANCE_STATUS" ,"MED_INC_QUAR_00",
           "NO_HSD_QUAR_00","CDCC_TOTAL_BEST","YEAR_OF_DIAGNOSIS" , "CS_SITESPECIFIC_FACTOR_1",
           "CS_SITESPECIFIC_FACTOR_8", "RX_SUMM_SURG_PRIM_SITE","RX_SUMM_SURGRAD_SEQ", "RX_SUMM_HORMONE", "fupmonth", "id","trx",
           "days2trx", "death" ,"tstage")) %>% na.omit()

names(dat_High)<-c("age","race","spanish","insurance","income","education",
                   "deyo","dyear","psa","gs", "surgsite","surgseq","hormone",
                   "fupmonth","id","trx","days2trx","death","tstage")


dat_High<- dat_High %>% mutate(dyearcat=ifelse(dyear %in% c(2004:2010), 1,2)) %>% 
  mutate_at(.vars=c("race","spanish","insurance","income","education",
                    "deyo","gs","trx","tstage","dyearcat"), as.factor)


#  psa rescale by 1/10 
dat_High$psa <- dat_High$psa/10

# use average days per month (365/12 = 30.42) to obtain time-to-treatment 
dat_High<- subset(dat_High, gs %in% c(6,7,8,9,10), drop=TRUE) #delete gs = 988, 998, 999
dat_High$gs <- factor(dat_High$gs, levels = c(6,7,8,9,10))
dat_High$month2trx <- dat_High$days2trx/30.42
dat_High$time <- dat_High$fupmonth - dat_High$month2trx
dat_High[dat_High$death==0,]$time<- dat_High[dat_High$death==0,]$time+epsilom ## add epsilom to censored times


# use the largest category as reference
dat_High$race <- relevel(factor(dat_High$race),ref='1')
dat_High$spanish <- relevel(factor(dat_High$spanish),ref='0')
dat_High$insurance <- relevel(factor(dat_High$insurance),ref='3')
dat_High$income <- relevel(factor(dat_High$income),ref='4')
dat_High$education <- relevel(factor(dat_High$education),ref='4')
dat_High$deyo <-relevel(factor(dat_High$deyo),ref='0')
#dat_High$gs <- relevel(factor(dat_High$gs),ref='8')
#dat_High$tstage <- relevel(factor(dat_High$tstage),ref='1')
dat_High$dyearcat <- relevel(factor(dat_High$dyearcat),ref='2')


## Now restrict to earlier diagnosis and further exclude Low risk (GS <=6 & T-stage = cT1) 
dat_early_diag<- dat_High[dat_High$dyearcat==1,] %>% filter(gs %in% c(7,8,9,10) | tstage %in% c("c2",  "c2A",  "c2B",  "c2C",   "c3",  "c3A",  "c3B",  "c3C",   "c4")) %>% 
  mutate(tstage=ifelse(tstage %in% c("c3", "c3A", "c3B", "c3C", "c4"),2,1)) %>% 
  mutate(tstage=relevel(factor(tstage),ref='1')) %>% 
  mutate(race= ifelse(race==1, "White","Non-White"),
         insurance = ifelse(insurance == 3, "Medicare",
                            ifelse(insurance == 1, "Private Insurance/Managed Care", "Other")),
         gs = ifelse(gs %in% c(6,7), "<=7", ifelse(gs== 8, 8, ">=9")))


source("Helper Functions.R")
source("Sim_data_func.R")

source("F_s_hat.R")





# First set seeds. 
if(length(args <- commandArgs(T))>0){  
  stopifnot(length(args)==1)
  case.id <- as.integer(args[[1]])
  message("running for parameter set ", case.id)
} 

Y <- dat_early_diag$time
d <- dat_early_diag$death
X=data.frame(dplyr::select(dat_early_diag, c(age,psa,tstage, deyo,gs,insurance, income, education, race))) 
X.new <- X
t= as.numeric(dat_early_diag$trx)-1
#t= 1-t ## for RP
trt=1 ## For EBRT+AD
gamma=0

s <- seq(0,130,length.out=30)

set.seed(case.id+100)
sim.size <- 1000

sim_dat <- Sim_data_func(Y=Y, d=d, t=t, X=X, sim.size=sim.size)

Y.sim <- sim_dat$Y
d.sim <- sim_dat$delta
X.sim <- data.frame(dplyr::select(sim_dat, c(age,psa,tstage, deyo,gs,insurance, income, education, race)))
t.sim <- sim_dat$t

## Original Estimate from simulated data
get_est <- F_s_hat(Y=Y.sim, d=d.sim, t=t.sim, X=X.sim, X.new=X.sim, trt=trt,s=s,gamma=gamma) 


## Setup bootstrap
nboot <-300
Boot.Est <- matrix(nrow=length(s), ncol=nboot)
Boot.SE <- matrix(nrow=length(s), ncol=nboot)
for (i in 1:nboot){
  boot.indices <- sample(1:sim.size, replace = T)
  Y.boot <- Y.sim[boot.indices]
  d.boot <- d.sim[boot.indices]
  X.boot <- X.sim[boot.indices,]
  t.boot <- t.sim[boot.indices]
  
  get_boot_est <- F_s_hat(Y=Y.boot, d=d.boot, t=t.boot, X=X.boot, X.new=X.boot, trt=trt,s=s,gamma=gamma)
  
  Boot.Est[,i] <- get_boot_est$pava.Est
  Boot.SE[,i] <-  sqrt(get_boot_est$var_error)
  
}
## Bootstrap confidence interval
Boot.CI <- apply(Boot.Est, 1, \(x)quantile(x, probs=c(0.025, 0.975)))

## Bootstrap-t confidence interval
abs_value <- sapply(1:ncol(Boot.Est), \(x) abs((Boot.Est[,x] - get_est$pava.Est)/ Boot.SE[,x])  )
abs_value <- ifelse(is.infinite(abs_value), NaN, abs_value)
C_j <- apply(abs_value, 1, \(x)quantile(x, probs=0.95, na.rm=T))

## Include logit transform CI
Trans_Est <- log(get_est$pava.Est/(1-get_est$pava.Est))
MOE <- C_j * sqrt((1/(get_est$pava.Est*(1-get_est$pava.Est)))^2 * get_est$var_error)
Symmetric.Boot.Lower.logit <- exp(Trans_Est - MOE)/(1+exp(Trans_Est - MOE))
Symmetric.Boot.Upper.logit <- exp(Trans_Est + MOE)/(1+exp(Trans_Est + MOE))

temp.save <- cbind(big.den=get_est$big.den,
                   pava.Est=get_est$pava.Est,
                   Est=get_est$Est,
                   Std_error=sqrt(get_est$var_error),
                   Lower.Limit=get_est$Lower.Limit,
                   Upper.Limit=get_est$Upper.Limit,
                   Trans_Est=get_est$Trans_Est,
                   MOE=get_est$MOE,
                   Lower.Limit.logit=get_est$Lower.Limit.logit,
                   Upper.Limit.logit=get_est$Upper.Limit.logit,
                   Boot.Lower.Limit= Boot.CI["2.5%",],
                   Boot.Upper.Limit= Boot.CI["97.5%",],
                   Symmetric.Boot.Lower.Limit = get_est$pava.Est - C_j*sqrt(get_est$var_error),
                   Symmetric.Boot.Upper.Limit = get_est$pava.Est + C_j*sqrt(get_est$var_error),
                   Symmetric.Boot.Lower.logit = Symmetric.Boot.Lower.logit,
                   Symmetric.Boot.Upper.logit = Symmetric.Boot.Upper.logit,
                   C_j = C_j)

write.csv(temp.save,paste("file",case.id, " .csv", sep=""))
