library(tidyverse)
library(lubridate)
library(broom)
library(boot)



load("enclave.eligible.Rdata")


# cases
enclave.eligible %>%
  group_by(MASS.cat=="MASS 3 or less") %>%
  tally() %>%
  mutate(pct = n/sum(n))

# severe outcomes COVID
enclave.eligible  %>%
  tally(covid_admit_death)
  
    #compare MASS 4 or more vs 3 or less
    enclave.eligible  %>%
      group_by(MASS.cat=="MASS 3 or less") %>%
      tally(covid_admit_death) %>%
      mutate(pct = n/sum(n))


# hospitalizations within 14d
enclave.eligible  %>%
  tally(covid_admit)

# deaths with 28d
enclave.eligible  %>%
  tally(death.4wks)

      #compare MASS 4 or more vs 3 or less
      enclave.eligible  %>%
        group_by(MASS.cat=="MASS 3 or less") %>%
        tally(death.4wks) %>%
        mutate(pct = n/sum(n))

# severe outcomes by timing of diagnosis
enclave.eligible  %>%
  group_by(outpatient.coviddx) %>%
  tally(covid_admit_death)

# outpatient COVID treatment
enclave.eligible  %>%
  tally(covid_treat)

# outpatient COVID treatment
enclave.eligible  %>%
  summarize(paxlovid =sum(paxlovid),
            molnupiravir = sum(molnupiravir),
            outpatient.remdesivir = sum(outpatient.remdesivir),
            bebtelovimab = sum(bebtelovimab))

# follow-up time
(enclave.eligible  %>%
 tally(days.fu))/365.25


##to use Chapman nearly unbiased estimator (NUE) adding + 1 to a, b, c as described in Hook and Regal, Table 1
a<- enclave.eligible  %>% #a
  group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  count() %>%
  ungroup() %>%
  select(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  mutate(outpatient.coviddx = 1, covid_admit_death = 1)

b<- enclave.eligible  %>% #b
  group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  count() %>%
  ungroup() %>%
  select(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  mutate(outpatient.coviddx = 1, covid_admit_death = 0)

c <- enclave.eligible  %>% #c
  group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  count() %>%
  ungroup() %>%
  select(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  mutate(outpatient.coviddx = 0, covid_admit_death = 1) 


# joining additional rows (+ 1 of the NUE) for each strata, limit to individuals with severe outcomes (second capture)
inpatients<- bind_rows(enclave.eligible , a, b, c) %>% filter(covid_admit_death == 1 )

save(inpatients, file="inpatients.Rdata")

rm(a,b,c)

# logistic model of probability of preceding outpatient dx
denom.fit <- glm(
  outpatient.coviddx ~ as.factor(MASS.cat) + as.factor(age.cat) + as.factor(mod.vax.status)  
  + as.factor(race.eth.collapsed)  + as.factor(highADI) + as.factor(covid_treat) ,
  family = binomial(),
  data = inpatients)

tidy(denom.fit, conf.int = TRUE, exponentiate = TRUE)

# assess significance of apriori selected predictors in univariate
gtsummary::tbl_regression(
  glm(outpatient.coviddx ~ as.factor(MASS.cat),
      family = binomial(),
      data = inpatients), exponentiate = TRUE) %>%
  gtsummary::add_global_p() # <0.001

gtsummary::tbl_regression(
  glm(outpatient.coviddx ~ as.factor(age.cat),
      family = binomial(),
      data = inpatients), exponentiate = TRUE) %>%
  gtsummary::add_global_p() # <0.001

gtsummary::tbl_regression(
  glm(outpatient.coviddx ~ as.factor(mod.vax.status) ,
      family = binomial(),
      data = inpatients), exponentiate = TRUE) %>%
  gtsummary::add_global_p() # <0.001

gtsummary::tbl_regression(
  glm(outpatient.coviddx ~ as.factor(race.eth.collapsed) ,
      family = binomial(),
      data = inpatients), exponentiate = TRUE) %>%
  gtsummary::add_global_p() # 0.02

gtsummary::tbl_regression(
  glm(outpatient.coviddx ~ as.factor(highADI) ,
      family = binomial(),
      data = inpatients), exponentiate = TRUE) %>%
  gtsummary::add_global_p() # <0.001

gtsummary::tbl_regression(
  glm(outpatient.coviddx ~ as.factor(covid_treat) ,
      family = binomial(),
      data = inpatients), exponentiate = TRUE) %>%
  gtsummary::add_global_p() # <0.001



#probability of outpatient diagnosis, then taking inverse
pd.outpt <- predict(denom.fit, type = "response")
inpatients$outpt.wts <- 1/pd.outpt

# developing strata-specific weights and the final aspect of the Chapman NUE
weights<- inpatients %>%
  group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  summarise(admit_death = n(),
            outpt.wts = mean(outpt.wts) - 1/n() ) %>% ungroup() # mean() used to pull a single value, not the mean, includes final -1 observation for Chapman formula


## calculate projected population and incidence in pseudopopulation 
enclave.eligible.wt <- left_join(
  enclave.eligible, weights, by = c("MASS.cat", "age.cat", "mod.vax.status", "race.eth.collapsed", "highADI", "covid_treat") ) %>% 
  mutate(outpt.wts = if_else(outpatient.coviddx == 1 & covid_admit_death == 0, outpt.wts, 1)) # use weights for outpatient dx (Chapman NUE) to estimate denominator, no weights for numerator as all observed

save(enclave.eligible.wt, file="enclave.eligible.wt.Rdata")

# projected cases
est.cases <-enclave.eligible.wt  %>%
  summarise(n = sum(if_else(outpatient.coviddx == 1, outpt.wts, 0), na.rm = TRUE))
est.cases

# projected individuals by comorbidity strata
enclave.eligible.wt  %>%
  group_by(MASS.cat) %>%
  summarise(n = sum(outpt.wts)) %>%
  ungroup() %>%
  mutate(pct = n / sum(n))



## outcomes
# descriptive statistics, overall
enclave.eligible.wt %>%
  summarize(covid_admit_death = sum(covid_admit_death*outpt.wts, na.rm = TRUE),
            death.4wks = sum(death.4wks*outpt.wts, na.rm = TRUE),
            covid_admit = sum(covid_admit*outpt.wts, na.rm = TRUE),
            total = sum(if_else(outpatient.coviddx == 1, outpt.wts, 0), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(death.pct = death.4wks/total,
         hosp.pct = covid_admit/total,
         severe.covid.pct = covid_admit_death/total)


# descriptive statistics in untreated, immunocompromised and/or MASS 6+ vs others
enclave.eligible.wt %>%
  filter(covid_treat == 0 ) %>%
  mutate(outpt.wts = if_else(outpatient.coviddx == 1 & covid_admit_death == 0, outpt.wts, 1)) %>%
  mutate(mass6_ic = if_else(MASS.cat ==  "MASS 4 and 5" | MASS.cat ==  "MASS 6 or greater" | 
                              MASS.cat == "Severe immunocompromise", 1, 0)) %>%
  group_by(mass6_ic) %>%
  summarize(n = sum(if_else(outpatient.coviddx == 1, outpt.wts, 0), na.rm=TRUE),
            risk = sum(covid_admit_death*outpt.wts, na.rm=TRUE)
            /sum(if_else(outpatient.coviddx == 1, outpt.wts, 0), na.rm=TRUE),
            mortality = sum(death.4wks*outpt.wts, na.rm=TRUE)
            /sum(if_else(outpatient.coviddx == 1, outpt.wts, 0), na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(pct.cases = n/sum(n),
         pct.deaths = (n*mortality)/sum(n*mortality))


## bootstrap 95% CI for total number of increased risk COVID cases

# create function
boot_totalcases <- function(d, i) {
  d2 <- d[i,]
  return(floor(sum(d2$outpt.wts, na.rm=TRUE)))
}

bootcases <- boot(enclave.eligible.wt, 
                  boot_totalcases, R=2099)  
tidy(bootcases, conf.int = TRUE,  conf.level = .95) 



## bootstrap 95%CI for strata-specific estimates, used in table 2

### No outpatient therapy for COVID
# create function to get percent covid_admit_death
boot_freq.severe <- function(d, i) {
  d2 <- d[i,]
  return(sum(d2$covid_admit_death , na.rm=TRUE)/sum(if_else(d2$outpatient.coviddx == 1, d2$outpt.wts, 0), na.rm=TRUE))
}
# limited data set for untreated, outpatient covid diagnosed patients (weights to full population)
enclave.eligible.wt.boot <- enclave.eligible.wt %>% 
  filter(covid_treat == 0)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat3.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat3.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
cat3.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat5.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat5.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
cat5.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat6.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat6.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
cat6.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
ihc.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
ihc.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
ihc.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)


enclave.eligible.wt %>%
  filter(covid_treat == 0) %>%
  group_by(MASS.cat, mod.vax.status) %>%
  summarise(severecovid = sum(covid_admit_death),
            total = floor(sum(outpt.wts))) %>%
  mutate(risk = severecovid/total)


severecovid.risk.est.untreat.all  <- rbind(cat3.unvax, cat3.vax.gt8mo, cat3.vax.lt8mo, 
                                           cat5.unvax, cat5.vax.gt8mo, cat5.vax.lt8mo,  
                                           cat6.unvax, cat6.vax.gt8mo, cat6.vax.lt8mo, 
                                           ihc.unvax, ihc.vax.gt8mo, ihc.vax.lt8mo) %>%
  mutate(statistic = signif(100*statistic, 2),
         conf.low = signif(100*conf.low , 2),
         conf.high = signif(100*conf.high , 2),
         risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
  select(MASS.cat, mod.vax.status, risk) %>%
  rename("MASS score" = MASS.cat, 
         "Vaccination status" = mod.vax.status,
         "No outpatient COVID treatment" = risk) 
save(severecovid.risk.est.untreat.all, file="severecovid.risk.est.untreat.all.Rdata")

rm(cat3.unvax, cat3.vax.gt8mo, cat3.vax.lt8mo, 
   cat5.unvax, cat5.vax.gt8mo, cat5.vax.lt8mo,  
   cat6.unvax, cat6.vax.gt8mo, cat6.vax.lt8mo, 
   ihc.unvax, ihc.vax.gt8mo, ihc.vax.lt8mo)

# TREATED
# create function to get percent covid_admit_death
boot_freq.severe <- function(d, i) {
  d2 <- d[i,]
  return(sum(d2$covid_admit_death, na.rm=TRUE)/sum(if_else(d2$outpatient.coviddx == 1, d2$outpt.wts, 0), na.rm=TRUE))
}
# limited data set for untreated, outpatient covid diagnosed patients (weights to full population)
enclave.eligible.wt.boot <- enclave.eligible.wt %>% 
  filter(covid_treat == 1)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat3.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat3.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
cat3.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat5.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat5.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
cat5.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat6.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
cat6.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
cat6.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
                 boot_freq.severe, R=2099)  
ihc.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
                 boot_freq.severe, R=2099)  
ihc.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)

bootfreq <- boot(enclave.eligible.wt.boot %>%
                   filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
                 boot_freq.severe, R=2099)  
ihc.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
  mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
  select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)



severecovid.risk.est.treat.all  <- rbind(cat3.unvax.tx, cat3.vax.gt8mo.tx, cat3.vax.lt8mo.tx, 
                                         cat5.unvax.tx, cat5.vax.gt8mo.tx, cat5.vax.lt8mo.tx,  
                                         cat6.unvax.tx, cat6.vax.gt8mo.tx, cat6.vax.lt8mo.tx, 
                                         ihc.unvax.tx, ihc.vax.gt8mo.tx, ihc.vax.lt8mo.tx) %>%
  mutate(statistic = signif(100*statistic, 2),
         conf.low = signif(100*conf.low , 2),
         conf.high = signif(100*conf.high , 2),
         risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
  select(MASS.cat, mod.vax.status, risk) %>%
  rename("MASS score" = MASS.cat, 
         "Vaccination status" = mod.vax.status,
         "Prescribed outpatient COVID treatment" = risk) 
save(severecovid.risk.est.treat.all, file="severecovid.risk.est.treat.all.Rdata")

rm(cat3.unvax.tx, cat3.vax.gt8mo.tx, cat3.vax.lt8mo.tx, 
   cat5.unvax.tx, cat5.vax.gt8mo.tx, cat5.vax.lt8mo.tx,  
   cat6.unvax.tx, cat6.vax.gt8mo.tx, cat6.vax.lt8mo.tx, 
   ihc.unvax.tx, ihc.vax.gt8mo.tx, ihc.vax.lt8mo.tx)

### sensitivity analyses

## non-parametric estimation of weights
nonpara.wts<- enclave.eligible %>%
  group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI,
           covid_treat) %>%
  summarize(a = sum(if_else(outpatient.coviddx == 1 & covid_admit_death == 1, 1, 0)),
            b = sum(if_else(outpatient.coviddx == 1 & covid_admit_death == 0, 1, 0)),
            c = sum(if_else(outpatient.coviddx == 0 & covid_admit_death == 1, 1, 0))) %>%
  ungroup() %>%
  mutate(NUE1 = ((b + 1)*(c + 1)/(a + 1))-1,
         NUE2 = a + b + c + ((b*c)/(a+1)),
         MLE1 = if_else(a != 0, a + b + c + (b*c)/a, 0),
         MLE2 = if_else(a != 0, (a + b)*(a + c)/a, 0),
         NUEd = (b*c)/(a+1),
         MLEd = if_else(a != 0,(b*c)/(a), 0),
         wt.nonpara = 1/((a + 1) / (a + c + 2))) %>% #summarise(NUE1 = sum(NUE1), NUE2 = sum(NUE2), MLE1 = sum(MLE1), MLE2 = sum(MLE2), NUEd = sum(NUEd), MLEd = sum(MLEd)) %>%
  select(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat, wt.nonpara)

enclave.nonpara.wt <- left_join(
  enclave.eligible, nonpara.wts, by = c("MASS.cat", "age.cat", "mod.vax.status", "race.eth.collapsed", "highADI", "covid_treat") ) %>%
  mutate(outpt.wts = if_else(outpatient.coviddx == 1 & covid_admit_death == 0, wt.nonpara, 1)) %>% # use weights for outpatient dx (Chapman NUE) to estimate denominator, no weights for numerator as all observed
  mutate(mod.vax.status = fct_relevel(mod.vax.status, "Unvaccinated or fewer than 3 doses", "Vaccinated, last dose or COVID-19 ≥ 8 months prior"))

# projected cases
est.cases.nonpara <- enclave.nonpara.wt  %>%
  summarise(n = sum(outpt.wts))

pct.nonpara.vs.para<- est.cases.nonpara/est.cases

rm(est.cases.nonpara, bootcases, denom.fit, bootfreq, inpatients)

# est.para <- cbind(enclave.eligible.wt %>%
#                        filter(covid_treat == 0) %>%
#                        group_by(MASS.cat, mod.vax.status) %>%
#                        summarize(covid_admit_death = round(sum(covid_admit_death),0),
#                                  cases = round(sum(outpt.wts),0),
#                                  untreat.risk = covid_admit_death/cases) %>%
#                        ungroup(),
#                      enclave.nonpara.wt %>%
#                        filter(covid_treat == 1) %>%
#                        group_by(MASS.cat, mod.vax.status) %>%
#                        summarize(covid_admit_death = sum(covid_admit_death),
#                                  cases = sum(outpt.wts),
#                                  treat.risk = covid_admit_death/cases) %>%
#                        ungroup() %>%
#                        select(treat.risk))
# 
# est.nonpara <- cbind(enclave.nonpara.wt %>%
#                 filter(covid_treat == 0) %>%
#                 group_by(MASS.cat, mod.vax.status) %>%
#                 summarize(covid_admit_death.nonpara = round(sum(covid_admit_death),0),
#                           cases.nonpara = round(sum(wt.nonpara),0),
#                           untreat.risk.nonpara = covid_admit_death.nonpara/cases.nonpara) %>%
#                     ungroup(),
#                 enclave.nonpara.wt %>%
#                 filter(covid_treat == 1) %>%
#                   group_by(MASS.cat, mod.vax.status) %>%
#                   summarize(covid_admit_death.nonpara = sum(covid_admit_death),
#                             cases.nonpara = sum(wt.nonpara),
#                             treat.risk.nonpara = covid_admit_death.nonpara/cases.nonpara) %>%
#                   ungroup() %>%
#                   select(treat.risk.nonpara))
# 
# # estimate total cases by non-parametric NUE
# boot_totalcases <- function(d, i) {
#   d2 <- d[i,]
#   return(floor(sum(d2$wt.nonpara, na.rm=TRUE)))
# }
# 
# bootcases.nonpara <- boot(enclave.nonpara.wt , 
#                   boot_totalcases, R=2099)  
# tidy(bootcases.nonpara, conf.int = TRUE,  conf.level = .95) 
# 
# 
# 
# ## bootstrap 95%CI for strata-specific estimates, use in table 2
# 
# ### No outpatient therapy for COVID
# # create function to get percent covid_admit_death
# boot_freq.severe <- function(d, i) {
#   d2 <- d[i,]
#   return(sum(d2$covid_admit_death, na.rm=TRUE)/sum(if_else(d2$outpatient.coviddx == 1, d2$wt.nonpara, 0), na.rm=TRUE))
# }
# # limited data set for untreated, outpatient covid diagnosed patients (weights to full population)
# enclave.eligible.wt.boot <- enclave.nonpara.wt%>% 
#   filter(covid_treat == 0)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat3.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat3.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# cat3.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat5.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat5.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# cat5.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat6.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat6.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# cat6.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# ihc.vax.lt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# ihc.vax.gt8mo<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# ihc.unvax<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# 
# enclave.nonpara.wt%>%
#   filter(covid_treat == 0) %>%
#   group_by(MASS.cat, mod.vax.status) %>%
#   summarise(severecovid = sum(covid_admit_death),
#             total = floor(sum(wt.nonpara))) %>%
#   mutate(risk = severecovid/total)
# 
# 
# severecovid.risk.est.untreat.all.nonpara  <- rbind(cat3.unvax, cat3.vax.gt8mo, cat3.vax.lt8mo, 
#                                            cat5.unvax, cat5.vax.gt8mo, cat5.vax.lt8mo,  
#                                            cat6.unvax, cat6.vax.gt8mo, cat6.vax.lt8mo, 
#                                            ihc.unvax, ihc.vax.gt8mo, ihc.vax.lt8mo) %>%
#   mutate(statistic = signif(100*statistic, 2),
#          conf.low = signif(100*conf.low , 2),
#          conf.high = signif(100*conf.high , 2),
#          risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
#          mod.vax.status = case_when(
#            mod.vax.status == "Unvaccinated or fewer than 3 doses" ~ "Not fully vaccinated (< 3 doses)",
#            TRUE ~ mod.vax.status
#          )) %>%
#   select(MASS.cat, mod.vax.status, risk) %>%
#   rename("MASS score" = MASS.cat, 
#          "Vaccination status" = mod.vax.status,
#          "No outpatient COVID treatment" = risk) 
# 
# # TREATED
# # create function to get percent covid_admit_death
# boot_freq.severe <- function(d, i) {
#   d2 <- d[i,]
#   return(sum(d2$covid_admit_death, na.rm=TRUE)/sum(if_else(d2$outpatient.coviddx == 1, d2$wt.nonpara, 0), na.rm=TRUE))
# }
# # limited data set for untreated, outpatient covid diagnosed patients (weights to full population)
# enclave.eligible.wt.boot <- enclave.nonpara.wt%>% 
#   filter(covid_treat == 1)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat3.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat3.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 3 or less" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# cat3.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 3 or less", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat5.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat5.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 4 and 5" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# cat5.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 4 and 5", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat6.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# cat6.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "MASS 6 or greater" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# cat6.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "MASS 6 or greater", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 < 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# ihc.vax.lt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 < 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Vaccinated, last dose or COVID-19 ≥ 8 months prior"), 
#                  boot_freq.severe, R=2099)  
# ihc.vax.gt8mo.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last dose or COVID-19 ≥ 8 months prior") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# bootfreq <- boot(enclave.eligible.wt.boot %>%
#                    filter(MASS.cat == "Severe immunocompromise" & mod.vax.status == "Unvaccinated or fewer than 3 doses"), 
#                  boot_freq.severe, R=2099)  
# ihc.unvax.tx<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
#   mutate(measure = "Covid Hosp or Death", MASS.cat = "Severe immunocompromise", mod.vax.status = "Unvaccinated or fewer than 3 doses") %>%
#   select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
# 
# 
# 
# severecovid.risk.est.treat.all.nonpara  <- rbind(cat3.unvax.tx, cat3.vax.gt8mo.tx, cat3.vax.lt8mo.tx, 
#                                          cat5.unvax.tx, cat5.vax.gt8mo.tx, cat5.vax.lt8mo.tx,  
#                                          cat6.unvax.tx, cat6.vax.gt8mo.tx, cat6.vax.lt8mo.tx, 
#                                          ihc.unvax.tx, ihc.vax.gt8mo.tx, ihc.vax.lt8mo.tx) %>%
#   mutate(statistic = signif(100*statistic, 2),
#          conf.low = signif(100*conf.low , 2),
#          conf.high = signif(100*conf.high , 2),
#          risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
#          mod.vax.status = case_when(
#            mod.vax.status == "Unvaccinated or fewer than 3 doses" ~ "Not fully vaccinated (< 3 doses)",
#            TRUE ~ mod.vax.status
#          )) %>%
#   select(MASS.cat, mod.vax.status, risk) %>%
#   rename("MASS score" = MASS.cat, 
#          "Vaccination status" = mod.vax.status,
#          "Prescribed outpatient COVID treatment" = risk) 
# 
rm(list = ls())
