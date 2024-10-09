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



## Using entire data to create Reference survival curves.
ncdb_puf <- read_sas("NCDB_PUF_DATA_Aug-15-2022/NCDB_PUF_DATA_Aug-15-2022/ncdb_puf.sas7bdat")

dat0_a <- ncdb_puf %>%
  dplyr::mutate(id = row_number()) |> 
  dplyr::filter(substr(TNM_CLIN_N, 1, 2) %in% c('c0') &
                  substr(TNM_CLIN_M, 1, 2) %in% c('c0'))

dat0_c <- dat0_a %>%
  filter(HISTOLOGY %in% c(8140, 8010) &
           BEHAVIOR != 2 &
           GRADE != 4 &
           DIAGNOSTIC_CONFIRMATION <= 6 &
           substr(SEQUENCE_NUMBER, 1, 2) == '00')

# For treatment only include this
dat0_d <- dat0_c %>%
  mutate(trx = ifelse(RX_SUMM_SURG_PRIM_SITE %in% c(50, 70), 1,
                      ifelse(PHASE_I_RT_MODALITY %in% c(1, 2, 3, 4, 5, 6) & RX_SUMM_HORMONE == 1 & PHASE_I_RT_VOLUME %in% c(64, 65), 2,"Other")))


dat0_e <- dat0_d %>%
  filter(!is.na(DX_LASTCONTACT_DEATH_MONTHS))

#For treatment only include this
dat0_f <- dat0_e %>%
  mutate(trxtime = case_when(
    trx == 1 ~ DX_DEFSURG_STARTED_DAYS,
    trx == 2 ~ DX_HORMONE_STARTED_DAYS))

dat0_g <- dat0_e %>%
  filter(!(#is.na(trxtime) | 
    TNM_CLIN_T == "")) %>%
  mutate(death = 1 - PUF_VITAL_STATUS,
         fupmonth = DX_LASTCONTACT_DEATH_MONTHS,
         #days2trx = trxtime,
         #dose_tot = RAD_REGIONAL_DOSE_CGY + RAD_BOOST_DOSE_CGY
  ) %>%
  filter(#DX_LASTCONTACT_DEATH_MONTHS*30.25>trxtime &
    RACE != 99 &
      SPANISH_HISPANIC_ORIGIN != 9 &
      !is.na(CROWFLY) &
      !is.na(FACILITY_TYPE_CD) &
      !is.na(FACILITY_LOCATION_CD))

dat0_b <- dat0_g %>%
  filter(YEAR_OF_DIAGNOSIS %in% c(2004:2010)) |> 
  filter(CS_SITESPECIFIC_FACTOR_1 <= 980) |> 
  filter(!(TNM_CLIN_T %in% c("CX","X",'cX', 'c0',88, "pA","pIS"))) |> 
  filter(!(CS_SITESPECIFIC_FACTOR_8 %in% c(988, 998, 999))) |> 
  mutate(CS_SITESPECIFIC_FACTOR_1 = ifelse (CS_SITESPECIFIC_FACTOR_1 >980, NA, CS_SITESPECIFIC_FACTOR_1)) |> 
  mutate(Risk_group = ifelse( ((CS_SITESPECIFIC_FACTOR_1 >= 200 & CS_SITESPECIFIC_FACTOR_1 <= 980) &
                                 ( CS_SITESPECIFIC_FACTOR_8 %in% c(8, 9, 10))  ) 
                              | (substr(TNM_CLIN_T, 1, 2) %in% c('c3', 'c4')), "High Risk", "Low Risk"),
         #Risk_group2 = ifelse(CS_SITESPECIFIC_FACTOR_8 %in% c(8, 9, 10), "High Risk", "Low Risk"), 
         # Risk_group2 = ifelse(substr(TNM_CLIN_T, 1, 2) %in% c('c3', 'c4'), "High Risk", "Low Risk"),
         #Risk_group2 = ifelse((CS_SITESPECIFIC_FACTOR_1 >= 200 & CS_SITESPECIFIC_FACTOR_1 <= 980), "High Risk", "Low Risk"),
         Risk_group2 = ifelse(TNM_CLIN_T %in% c('1', '1A', "1C","2","2C","c1","c1A", "c1A2","c1B", "c1B1", "c1B2", "c1C","C1C","cA1", "c2", "C2", "c2A", "C2A") &
                                (CS_SITESPECIFIC_FACTOR_1 < 100) & CS_SITESPECIFIC_FACTOR_8 <= 6, "Low Risk", ifelse(substr(TNM_CLIN_T, 1, 2) %in% c('c3', 'c4') |
                                                                                                                       (CS_SITESPECIFIC_FACTOR_1 >= 200 & CS_SITESPECIFIC_FACTOR_1 <= 980) |
                                                                                                                       CS_SITESPECIFIC_FACTOR_8 %in% c(8, 9, 10), "High Risk", "Intermediate Low Risk")),
         Risk_group3 = ifelse(TNM_CLIN_T %in% c('1', '1A', "1C","2","2C","c1","c1A", "c1A2","c1B", "c1B1", "c1B2", "c1C","C1C","cA1", "c2", "C2", "c2A", "C2A") &
                                (CS_SITESPECIFIC_FACTOR_1 < 100) & CS_SITESPECIFIC_FACTOR_8 <= 6, "Low Risk", ifelse(substr(TNM_CLIN_T, 1, 2) %in% c('c3', 'c4') , "High Risk", "Intermediate Low Risk")),
         #days2trx= ifelse(days2trx==0, .5, days2trx),
         #month2trx = days2trx/30.42,
         time = fupmonth #- month2trx
  ) |> 
  mutate(time = ifelse (death==0,time+0.0001,time)) 

# Recategorise some of the data
dat0_bcat <- dat0_b |> 
  mutate(RACE= factor(ifelse(RACE==1, "White", "Non-White"), levels=c("White","Non-White")),
         INSURANCE_STATUS = factor(ifelse(INSURANCE_STATUS == 3, "Medicare", ifelse(INSURANCE_STATUS == 1, "Private Insurance/Managed Care", "Other")), levels=c("Private Insurance/Managed Care", "Medicare", "Other")),
         MED_INC_QUAR_00= factor(MED_INC_QUAR_00, levels=c("1", "2", "3","4"), labels= c("<30,000", "30,000-34,999", "35,000-49,999", ">45,000")),
         NO_HSD_QUAR_00 = factor(NO_HSD_QUAR_00, levels=c("1", "2", "3","4"), labels= c("<14", "14-19.99", "20-28.99", ">29")),
         #CS_SITESPECIFIC_FACTOR_8= ifelse(CS_SITESPECIFIC_FACTOR_8>=8, 1, 0), ## indicator for gs for high risk >=8
         CS_SITESPECIFIC_FACTOR_8 = factor(ifelse(CS_SITESPECIFIC_FACTOR_8 <=7, "<=7", ifelse(CS_SITESPECIFIC_FACTOR_8== 8, 8, ">=9")), levels=c("<=7", "8", ">=9")),
         #CS_SITESPECIFIC_FACTOR_1 = ifelse(CS_SITESPECIFIC_FACTOR_1 >= 200, 1, 0), ## indicator for psa for high risk >=200,
         CS_SITESPECIFIC_FACTOR_1 = CS_SITESPECIFIC_FACTOR_1/10,
         tstage = ifelse(TNM_CLIN_T %in% c("c3", "c3A", "c3B", "c3C", "c4"),"1","0")) |> ## indicator for tstage for high risk >=cT3
  mutate(tstage = factor(tstage, levels=c("0", "1"), label= c("<=cT2", ">=cT3")),
         CDCC_TOTAL_BEST = ifelse(CDCC_TOTAL_BEST ==0,"0",">=1"),
         TNM_CLIN_N= ifelse(substr(TNM_CLIN_N, 1, 2) %in% c('c0'),"0","1"),
         TNM_CLIN_M=ifelse(substr(TNM_CLIN_M, 1, 2) %in% c('c0'), "0","1")
  ) %>% 
  
  select(c("AGE", "RACE" , "SPANISH_HISPANIC_ORIGIN",  "INSURANCE_STATUS" ,"MED_INC_QUAR_00","NO_HSD_QUAR_00","CDCC_TOTAL_BEST", "CS_SITESPECIFIC_FACTOR_1",
           "CS_SITESPECIFIC_FACTOR_8", "RX_SUMM_SURG_PRIM_SITE","RX_SUMM_SURGRAD_SEQ", "RX_SUMM_HORMONE", "time","fupmonth", "id","death","trx" ,"tstage", Risk_group, TNM_CLIN_N, TNM_CLIN_M, Risk_group2, Risk_group3)) %>% na.omit()

names(dat0_bcat)<-c("age","race","spanish","insurance","income","education",
                    "deyo","psa","gs", "surgsite","surgseq","hormone", "time",
                    "fupmonth","id","death", "trx", "tstage", "Risk_group", 
                    "TNM_CLIN_N", "TNM_CLIN_M", "Risk_group2", "Risk_group3")



## Survival curves for reference
pred_time <- seq(0, 130, length.out = 30)

newdata_low <- dat0_bcat |> filter(Risk_group3 == "Intermediate Low Risk") |> filter(TNM_CLIN_N == "0") |> filter(TNM_CLIN_M == "0") 
S_t_low = coxph( as.formula(paste("Surv(time = time, event = death) ~", var)),  ties="breslow",
                 data = newdata_low)
fit_low <- survfit(S_t_low, stype=1, newdata=newdata_low, se.fit = FALSE)
low <- rowMeans(summary(fit_low, times = pred_time, extend = TRUE)$surv)

newdata_high <- dat0_bcat |> filter(Risk_group == "High Risk")  |> filter(TNM_CLIN_N == "0") |> filter(TNM_CLIN_M == "0") 
S_t_high = coxph( as.formula(paste("Surv(time = time, event = death) ~", var)),  ties="breslow",
                  data = newdata_high)
fit_high <- survfit(S_t_high, stype=1, newdata=newdata_high, se.fit = FALSE)
high <- rowMeans(summary(fit_high, times = pred_time, extend = TRUE)$surv)

ref_dat <- tibble(`Survival Probability` = c(low, high), 
                  Time = rep(pred_time, 2),
                  Estimates = rep(c("Lower Risk", "Higher Risk"), each = 30),
                  Treatment = rep(c("EBRT + AD", "RP"), each = 30 ))
