library(survival)
library(dplyr)
library(rms)
library(splines)

## Read in Data
epsilom <- 0.00000001
data_9_11_22 <- read.csv("data_9_11_22.csv") 


## Data manupulation
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
  tstage=ifelse(TNM_CLIN_T == 88, NA,TNM_CLIN_T),
  CDCC_TOTAL_BEST = ifelse(CDCC_TOTAL_BEST ==0,"0",">=1")
  ) %>% 
  
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
dat_High$dyearcat <- relevel(factor(dat_High$dyearcat),ref='2')


## Now restrict to earlier diagnosis and further exclude Low risk (GS <=6 & T-stage = cT1) 
dat_early_diag<- dat_High[dat_High$dyearcat==1,] %>% filter(gs %in% c(7,8,9,10) | tstage %in% c("c2",  "c2A",  "c2B",  "c2C",   "c3",  "c3A",  "c3B",  "c3C",   "c4")) %>% 
  mutate(tstage=ifelse(tstage %in% c("c3", "c3A", "c3B", "c3C", "c4"),2,1)) %>% 
  mutate(tstage=relevel(factor(tstage),ref='1')) %>% 
  mutate(race= ifelse(race==1, "White","Non-White"),
         insurance = ifelse(insurance == 3, "Medicare",
                            ifelse(insurance == 1, "Private Insurance/Managed Care", "Other")),
         gs = ifelse(gs %in% c(6,7), "<=7", ifelse(gs== 8, 8, ">=9")))


Y <- dat_early_diag$time
d <- dat_early_diag$death
X=data.frame(dplyr::select(dat_early_diag, c(age,psa,tstage, deyo,gs,insurance, income, education, race))) 
X.new <- X
t= as.numeric(dat_early_diag$trx)-1
trt=1

RP <- F_s_hat(Y=Y, d=d, t=t, X=X, X.new=X, trt=0,s=NULL,gamma=gamma)
EBRT <- F_s_hat(Y=Y, d=d, t=t, X=X, X.new=X, trt=1,s=NULL,gamma=gamma)
