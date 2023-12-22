library(tidyverse)
library(lubridate)
library(broom)
library(marginaleffects)
rm(list = ls())


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

#assess positivity assumption
quantile(predict(denom.fit, type = "response"))

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
weights <- inpatients %>%
  group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat) %>%
  summarise(admit_death = n(),
            outpt.wts = mean(outpt.wts) - 1/n() ) %>% ungroup() # mean() used to pull a single value, not the mean, includes final -1 observation for Chapman formula


## calculate projected population and incidence in pseudopopulation 
enclave.eligible.wt <- left_join(
  enclave.eligible, weights, by = c("MASS.cat", "age.cat", "mod.vax.status", "race.eth.collapsed", "highADI", "covid_treat") ) %>% 
  mutate(outpt.wts = if_else(outpatient.coviddx == 1 & covid_admit_death == 0, outpt.wts, 1)) # use weights for outpatient dx (Chapman NUE) to estimate denominator, no weights for numerator as all observed

save(enclave.eligible.wt, file="enclave.eligible.wt.Rdata")

load("enclave.eligible.wt.Rdata")
enclave.eligible.wt %>% filter(is.na(outpt.wts))

# projected cases
est.cases <-enclave.eligible.wt  %>%
  summarise(n = sum(outpt.wts), na.rm = TRUE)
  #  summarise(n = sum(if_else(outpatient.coviddx == 1, outpt.wts, 1), na.rm = TRUE))
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
            total = sum(if_else(outpatient.coviddx == 1, outpt.wts, 1), na.rm = TRUE)) %>%
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

## estimates for results section

#total projected cases
sum(enclave.eligible.wt$outpt.wts)

#total undocumented cases
sum(enclave.eligible.wt$outpt.wts)- count(enclave.eligible.wt)

#ratio of documented:projected
sum(enclave.eligible.wt$outpt.wts) / count(enclave.eligible.wt)


#estimated severe outcomes and 95%CI
avg_predictions(
  geepack::geeglm(
  covid_admit_death ~ 1,
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt )
  )

#estimated hospitalization and 95%CI
avg_predictions(
  geepack::geeglm(
    covid_admit ~ 1,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt )
)


#estimated death and 95%CI
avg_predictions(
  geepack::geeglm(
    death.4wks ~ 1,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt )
)

# risk among vaccinated with MASS 3 or less
avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 3 or less"))) 


# risk among vaccinated with MASS 4 or more
avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 4 and 5", "MASS 6 or greater", "Severe immunocompromise"))) 


# risk among vaccinated with MASS 4+
avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 4 and 5", "MASS 6 or greater", "Severe immunocompromise"),
                     mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior"))) 


# risk among vaccinated with MASS 3 or less
avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 3 or less"),
                     mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior"))) 


# assess prictors for treatment
tidy(geepack::geeglm(
  covid_treat ~ MASS.cat + age.cat + mod.vax.status  + race.eth.collapsed + highADI,
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt  ), exponentiate = TRUE)



## estimate risk by comorbidity, vaccination, and treatment strata

## Untreated
# poisson with robust standard errors, including interaction with age.cat
#using outpatient treatment weights (inference to population with covid and at risk for hosp/death)
covid.risk.untreated <- geepack::geeglm(
  covid_admit_death ~ (MASS.cat + age.cat + mod.vax.status  + race.eth.collapsed + highADI)* age.cat,
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt %>% filter(covid_treat == 0) )


#obtaining est (95%CI) from g-estimation within strata of MASS.cat, vaccination, treatment (age, ADI, race/eth at dataset mean)
risks.untreated <-tibble(predictions(covid.risk.untreated, 
            newdata = datagrid(MASS.cat = c("MASS 3 or less", "MASS 4 and 5", "MASS 6 or greater", "Severe immunocompromise"),
                               mod.vax.status = c("Unvaccinated or fewer than 3 doses", "Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior"))) %>% arrange(MASS.cat, mod.vax.status)) %>%
    mutate(statistic = signif(100*estimate, 2),
           conf.low = signif(100*conf.low, 1),
           conf.high = signif(100*conf.high , 1),
           risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
  select(MASS.cat, mod.vax.status, risk) 



## Treated
# poisson with robust standard errors, including interaction with age.cat
#using outpatient treatment weights (inference to population with covid and at risk for hosp/death)

covid.risk.treated <- geepack::geeglm(
  covid_admit_death ~ (MASS.cat + age.cat + mod.vax.status  + race.eth.collapsed + highADI)* age.cat,
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt %>% filter(covid_treat == 1) )


#obtaining est (95%CI) from g-estimation within strata of MASS.cat, vaccination, treatment (age, ADI, race/eth at dataset mean)
risks.treated <-tibble(predictions(covid.risk.treated, 
                           newdata = datagrid(MASS.cat = c("MASS 3 or less", "MASS 4 and 5", "MASS 6 or greater", "Severe immunocompromise"),
                                              mod.vax.status = c("Unvaccinated or fewer than 3 doses", "Vaccinated, last dose or COVID-19 ≥ 8 months prior", 
                                                                 "Vaccinated, last dose or COVID-19 < 8 months prior"))) %>% arrange(MASS.cat, mod.vax.status)) %>%
  mutate(statistic = signif(100*estimate, 2),
         conf.low = signif(100*conf.low, 1),
         conf.high = signif(100*conf.high , 1),
         risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
  select(MASS.cat, mod.vax.status, risk)




################################
## RD/NNT within MASS.cat strata
################################


#MASS 3 or less
#ipw for covid_treat
MASS3orless <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 3 or less")
denom.fit <- glm( covid_treat ~ as.factor(age.cat) + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed)   + as.factor(highADI),
                  family = binomial(),
                  data = MASS3orless)

# denominator of weights
pd.treat <- predict(denom.fit, type = "response")

# numerator of stabilized weights, overall probabilty of outpatient treatment
numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS3orless)
pn.treat <- predict(numer.fit, type = "response")

# calculate weights
MASS3orless$sw.treat <-if_else(MASS3orless$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                               (pn.treat / pd.treat))
MASS3orless$rd.weight <- MASS3orless$sw.treat * MASS3orless$outpt.wts

quantile(MASS3orless$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS3orless$outpt.wts)

quantile(MASS3orless$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS3orless$sw.treat)

quantile(MASS3orless$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS3orless$rd.weight)

# estimate risk difference among MASS3orless
treat.MASS3orless <- geepack::geeglm(
  covid_admit_death ~ covid_treat,
  family = poisson(link=log),
  weights = rd.weight,
  id = seqn,
  corstr = "independence",
  data = MASS3orless )

rd.MASS3orless <- tibble(avg_comparisons(treat.MASS3orless)) %>% 
  mutate(MASS.cat = "MASS 3 or less",
         NNT = signif(-1/estimate, 2),
         NNT.low = signif(-1/conf.low,2),
         NNT.high = if_else(conf.high >= 0, "Inf", as.character(signif(-1/conf.high,2))),
         NNT.CI = paste0(NNT, " (", NNT.low, " to ", NNT.high, ")" ),
         statistic = signif(100*estimate, 2),
                conf.low = signif(100*conf.low, 1),
                conf.high = signif(100*conf.high , 1),
                rd = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
         rd.nnt = paste0("RD ", rd, "\nNNTB ", NNT.CI  )) %>%
           select(MASS.cat, rd.nnt)

#MASS 4 and 5 
#ipw for covid_treat
MASS4and5 <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 4 and 5")
denom.fit <- glm( covid_treat ~ as.factor(age.cat) + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed)   + as.factor(highADI),
                  family = binomial(),
                  data = MASS4and5)

# denominator of weights
pd.treat <- predict(denom.fit, type = "response")

# numerator of stabilized weights, overall probabilty of outpatient treatment
numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS4and5)
pn.treat <- predict(numer.fit, type = "response")

# calculate weights
MASS4and5$sw.treat <-if_else(MASS4and5$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                               (pn.treat / pd.treat))
MASS4and5$rd.weight <- MASS4and5$sw.treat * MASS4and5$outpt.wts

quantile(MASS4and5$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS4and5$outpt.wts)

quantile(MASS4and5$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS4and5$sw.treat)

quantile(MASS4and5$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS4and5$rd.weight)

# estimate risk difference among MASS4and5
treat.MASS4and5 <- geepack::geeglm(
  covid_admit_death ~ covid_treat,
  family = poisson(link=log),
  weights = rd.weight,
  id = seqn,
  corstr = "independence",
  data = MASS4and5 )

rd.MASS4and5 <- tibble(avg_comparisons(treat.MASS4and5)) %>% 
  mutate(MASS.cat = "MASS 4 and 5",
         NNT = signif(-1/estimate, 2),
         NNT.low = signif(-1/conf.low,2),
         NNT.high = if_else(conf.high >= 0, "Inf", as.character(signif(-1/conf.high,2))),
         NNT.CI = paste0(NNT, " (", NNT.low, " to ", NNT.high, ")" ),
         statistic = signif(100*estimate, 2),
         conf.low = signif(100*conf.low, 1),
         conf.high = signif(100*conf.high , 1),
         rd = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
         rd.nnt = paste0("RD ", rd, "\nNNTB ", NNT.CI  )) %>%
  select(MASS.cat, rd.nnt)

#MASS 3 or less
#ipw for covid_treat
MASS3orless <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 3 or less")
denom.fit <- glm( covid_treat ~ as.factor(age.cat) + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed)   + as.factor(highADI),
                  family = binomial(),
                  data = MASS3orless)

# denominator of weights
pd.treat <- predict(denom.fit, type = "response")

# numerator of stabilized weights, overall probabilty of outpatient treatment
numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS3orless)
pn.treat <- predict(numer.fit, type = "response")

# calculate weights
MASS3orless$sw.treat <-if_else(MASS3orless$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                               (pn.treat / pd.treat))
MASS3orless$rd.weight <- MASS3orless$sw.treat * MASS3orless$outpt.wts

quantile(MASS3orless$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS3orless$outpt.wts)

quantile(MASS3orless$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS3orless$sw.treat)

quantile(MASS3orless$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS3orless$rd.weight)

# estimate risk difference among MASS3orless
treat.MASS3orless <- geepack::geeglm(
  covid_admit_death ~ covid_treat,
  family = poisson(link=log),
  weights = rd.weight,
  id = seqn,
  corstr = "independence",
  data = MASS3orless )

rd.MASS3orless <- tibble(avg_comparisons(treat.MASS3orless)) %>% 
  mutate(MASS.cat = "MASS 3 or less",
         NNT = signif(-1/estimate, 2),
         NNT.low = signif(-1/conf.low,2),
         NNT.high = if_else(conf.high >= 0, "Inf", as.character(signif(-1/conf.high,2))),
         NNT.CI = paste0(NNT, " (", NNT.low, " to ", NNT.high, ")" ),
         statistic = signif(100*estimate, 2),
                conf.low = signif(100*conf.low, 1),
                conf.high = signif(100*conf.high , 1),
                rd = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
         rd.nnt = paste0("RD ", rd, "\nNNTB ", NNT.CI  )) %>%
           select(MASS.cat, rd.nnt)

#MASS 4 and 5 
#ipw for covid_treat
MASS4and5 <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 4 and 5")
denom.fit <- glm( covid_treat ~ as.factor(age.cat) + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed)   + as.factor(highADI),
                  family = binomial(),
                  data = MASS4and5)

# denominator of weights
pd.treat <- predict(denom.fit, type = "response")

# numerator of stabilized weights, overall probabilty of outpatient treatment
numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS4and5)
pn.treat <- predict(numer.fit, type = "response")

# calculate weights
MASS4and5$sw.treat <-if_else(MASS4and5$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                               (pn.treat / pd.treat))
MASS4and5$rd.weight <- MASS4and5$sw.treat * MASS4and5$outpt.wts

quantile(MASS4and5$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS4and5$outpt.wts)

quantile(MASS4and5$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS4and5$sw.treat)

quantile(MASS4and5$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS4and5$rd.weight)

# estimate risk difference among MASS4and5
treat.MASS4and5 <- geepack::geeglm(
  covid_admit_death ~ covid_treat,
  family = poisson(link=log),
  weights = rd.weight,
  id = seqn,
  corstr = "independence",
  data = MASS4and5 )

rd.MASS4and5 <- tibble(avg_comparisons(treat.MASS4and5)) %>% 
  mutate(MASS.cat = "MASS 4 and 5",
         NNT = signif(-1/estimate, 2),
         NNT.low = signif(-1/conf.low,2),
         NNT.high = if_else(conf.high >= 0, "Inf", as.character(signif(-1/conf.high,2))),
         NNT.CI = paste0(NNT, " (", NNT.low, " to ", NNT.high, ")" ),
         statistic = signif(100*estimate, 2),
         conf.low = signif(100*conf.low, 1),
         conf.high = signif(100*conf.high , 1),
         rd = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
         rd.nnt = paste0("RD ", rd, "\nNNTB ", NNT.CI  )) %>%
  select(MASS.cat, rd.nnt)



#MASS 6 or greater
#ipw for covid_treat
MASS6greater <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 6 or greater")
denom.fit <- glm( covid_treat ~ as.factor(age.cat) + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed)   + as.factor(highADI),
                  family = binomial(),
                  data = MASS6greater)

# denominator of weights
pd.treat <- predict(denom.fit, type = "response")

# numerator of stabilized weights, overall probabilty of outpatient treatment
numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS6greater)
pn.treat <- predict(numer.fit, type = "response")

# calculate weights
MASS6greater$sw.treat <-if_else(MASS6greater$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                             (pn.treat / pd.treat))
MASS6greater$rd.weight <- MASS6greater$sw.treat * MASS6greater$outpt.wts

quantile(MASS6greater$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS6greater$outpt.wts)

quantile(MASS6greater$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS6greater$sw.treat)

quantile(MASS6greater$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
mean(MASS6greater$rd.weight)

# estimate risk difference among MASS6greater
treat.MASS6greater <- geepack::geeglm(
  covid_admit_death ~ covid_treat,
  family = poisson(link=log),
  weights = rd.weight,
  id = seqn,
  corstr = "independence",
  data = MASS6greater )

rd.MASS6greater <- tibble(avg_comparisons(treat.MASS6greater)) %>% 
  mutate(MASS.cat = "MASS 6 or greater",
         NNT = signif(-1/estimate, 2),
         NNT.low = signif(-1/conf.low,2),
         NNT.high = if_else(conf.high >= 0, "Inf", as.character(signif(-1/conf.high,2))),
         NNT.CI = paste0(NNT, " (", NNT.low, " to ", NNT.high, ")" ),
         statistic = signif(100*estimate, 2),
         conf.low = signif(100*conf.low, 1),
         conf.high = signif(100*conf.high , 1),
         rd = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
         rd.nnt = paste0("RD ", rd, "\nNNTB ", NNT.CI  )) %>%
  select(MASS.cat, rd.nnt)


#MASS 6 or greater
#ipw for covid_treat
severeIC <- enclave.eligible.wt %>% filter(MASS.cat == "Severe immunocompromise")
denom.fit <- glm( covid_treat ~ as.factor(age.cat) + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed)   + as.factor(highADI),
                  family = binomial(),
                  data = severeIC)

# denominator of weights
pd.treat <- predict(denom.fit, type = "response")

# numerator of stabilized weights, overall probabilty of outpatient treatment
numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = severeIC)
pn.treat <- predict(numer.fit, type = "response")

# calculate weights
severeIC$sw.treat <-if_else(severeIC$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                                (pn.treat / pd.treat))
severeIC$rd.weight <- severeIC$sw.treat * severeIC$outpt.wts

quantile(severeIC$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
mean(severeIC$outpt.wts)

quantile(severeIC$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
mean(severeIC$sw.treat)

quantile(severeIC$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
mean(severeIC$rd.weight)

# estimate risk difference among severeIC
treat.severeIC <- geepack::geeglm(
  covid_admit_death ~ covid_treat,
  family = poisson(link=log),
  weights = rd.weight,
  id = seqn,
  corstr = "independence",
  data = severeIC )

rd.severeIC <- tibble(avg_comparisons(treat.severeIC)) %>% 
  mutate(MASS.cat = "Severe immunocompromise",
         NNT = signif(-1/estimate, 2),
         NNT.low = signif(-1/conf.low,2),
         NNT.high = if_else(conf.high >= 0, "Inf", as.character(signif(-1/conf.high,2))),
         NNT.CI = paste0(NNT, " (", NNT.low, " to ", NNT.high, ")" ),
         statistic = signif(100*estimate, 2),
         conf.low = signif(100*conf.low, 1),
         conf.high = signif(100*conf.high , 1),
         rd = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)"),
         rd.nnt = paste0("RD ", rd, "\nNNTB ", NNT.CI  )) %>%
  select(MASS.cat, rd.nnt)


risk.table <- 
  left_join(
    left_join(risks.untreated, 
              risks.treated, 
              by = c("MASS.cat", "mod.vax.status")),
    rbind( rd.MASS3orless, rd.MASS4and5, rd.MASS6greater, rd.severeIC),
              by = "MASS.cat")%>%
  rename(risks.untreated = risk.x, risk.treated= risk.y)
save(risk.table , file="risk.table.Rdata")


rm(list = ls())
#Prior code
# 
# risks.untreated <-tibble(predictions(covid.risk.msm, 
#                                      newdata = datagrid(MASS.cat = c("MASS 3 or less", "MASS 4 and 5", "MASS 6 or greater", "Severe immunocompromise"),
#                                                         mod.vax.status = c("Unvaccinated or fewer than 3 doses", "Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior"))) %>% arrange(MASS.cat, mod.vax.status)) %>%
#   mutate(statistic = signif(100*estimate, 2),
#          conf.low = signif(100*conf.low, 2),
#          conf.high = signif(100*conf.high , 2),
#          risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
#   select(MASS.cat, mod.vax.status, risk)
# 
# 
# 
# 
# 
# 
# ## bootstrap 95% CI for total number of increased risk COVID cases
# 
# # create function
# boot_totalcases <- function(d, i) {
#   d2 <- d[i,]
#   return(floor(sum(d2$outpt.wts, na.rm=TRUE)))
# }
# 
# bootcases <- boot(enclave.eligible.wt,
#                   boot_totalcases, R=2099)
# tidy(bootcases, conf.int = TRUE,  conf.level = .95)
# 
# 
# 
# ## bootstrap 95%CI for strata-specific estimates, used in table 2
# 
# ### No outpatient therapy for COVID
# # create function to get percent covid_admit_death
# boot_freq.severe <- function(d, i) {
#   d2 <- d[i,]
#   return(sum(d2$covid_admit_death , na.rm=TRUE)/sum(if_else(d2$outpatient.coviddx == 1, d2$outpt.wts, 0), na.rm=TRUE))
# }
# # limited data set for untreated, outpatient covid diagnosed patients (weights to full population)
# enclave.eligible.wt.boot <- enclave.eligible.wt %>%
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
# enclave.eligible.wt %>%
#   filter(covid_treat == 0) %>%
#   group_by(MASS.cat, mod.vax.status) %>%
#   summarise(severecovid = sum(covid_admit_death),
#             total = floor(sum(outpt.wts))) %>%
#   mutate(risk = severecovid/total)
# 
# 
# severecovid.risk.est.untreat.all  <- rbind(cat3.unvax, cat3.vax.gt8mo, cat3.vax.lt8mo,
#                                            cat5.unvax, cat5.vax.gt8mo, cat5.vax.lt8mo,
#                                            cat6.unvax, cat6.vax.gt8mo, cat6.vax.lt8mo,
#                                            ihc.unvax, ihc.vax.gt8mo, ihc.vax.lt8mo) %>%
#   mutate(statistic = signif(100*statistic, 2),
#          conf.low = signif(100*conf.low , 2),
#          conf.high = signif(100*conf.high , 2),
#          risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
#   select(MASS.cat, mod.vax.status, risk) %>%
#   rename("MASS score" = MASS.cat,
#          "Vaccination status" = mod.vax.status,
#          "No outpatient COVID treatment" = risk)
# save(severecovid.risk.est.untreat.all, file="severecovid.risk.est.untreat.all.Rdata")
# 
# rm(cat3.unvax, cat3.vax.gt8mo, cat3.vax.lt8mo,
#    cat5.unvax, cat5.vax.gt8mo, cat5.vax.lt8mo,
#    cat6.unvax, cat6.vax.gt8mo, cat6.vax.lt8mo,
#    ihc.unvax, ihc.vax.gt8mo, ihc.vax.lt8mo)
# 
# # TREATED
# # create function to get percent covid_admit_death
# boot_freq.severe <- function(d, i) {
#   d2 <- d[i,]
#   return(sum(d2$covid_admit_death, na.rm=TRUE)/sum(if_else(d2$outpatient.coviddx == 1, d2$outpt.wts, 0), na.rm=TRUE))
# }
# # limited data set for untreated, outpatient covid diagnosed patients (weights to full population)
# enclave.eligible.wt.boot <- enclave.eligible.wt %>%
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
# severecovid.risk.est.treat.all  <- rbind(cat3.unvax.tx, cat3.vax.gt8mo.tx, cat3.vax.lt8mo.tx,
#                                          cat5.unvax.tx, cat5.vax.gt8mo.tx, cat5.vax.lt8mo.tx,
#                                          cat6.unvax.tx, cat6.vax.gt8mo.tx, cat6.vax.lt8mo.tx,
#                                          ihc.unvax.tx, ihc.vax.gt8mo.tx, ihc.vax.lt8mo.tx) %>%
#   mutate(statistic = signif(100*statistic, 2),
#          conf.low = signif(100*conf.low , 2),
#          conf.high = signif(100*conf.high , 2),
#          risk = paste0(statistic, "% (", conf.low, "% to ", conf.high, "%)")) %>%
#   select(MASS.cat, mod.vax.status, risk) %>%
#   rename("MASS score" = MASS.cat,
#          "Vaccination status" = mod.vax.status,
#          "Prescribed outpatient COVID treatment" = risk)
# save(severecovid.risk.est.treat.all, file="severecovid.risk.est.treat.all.Rdata")
# 
# rm(cat3.unvax.tx, cat3.vax.gt8mo.tx, cat3.vax.lt8mo.tx,
#    cat5.unvax.tx, cat5.vax.gt8mo.tx, cat5.vax.lt8mo.tx,
#    cat6.unvax.tx, cat6.vax.gt8mo.tx, cat6.vax.lt8mo.tx,
#    ihc.unvax.tx, ihc.vax.gt8mo.tx, ihc.vax.lt8mo.tx)
# 
# ### sensitivity analyses
# 
# ## non-parametric estimation of weights
# nonpara.wts<- enclave.eligible %>%
#   group_by(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI,
#            covid_treat) %>%
#   summarize(a = sum(if_else(outpatient.coviddx == 1 & covid_admit_death == 1, 1, 0)),
#             b = sum(if_else(outpatient.coviddx == 1 & covid_admit_death == 0, 1, 0)),
#             c = sum(if_else(outpatient.coviddx == 0 & covid_admit_death == 1, 1, 0))) %>%
#   ungroup() %>%
#   mutate(NUE1 = ((b + 1)*(c + 1)/(a + 1))-1,
#          NUE2 = a + b + c + ((b*c)/(a+1)),
#          MLE1 = if_else(a != 0, a + b + c + (b*c)/a, 0),
#          MLE2 = if_else(a != 0, (a + b)*(a + c)/a, 0),
#          NUEd = (b*c)/(a+1),
#          MLEd = if_else(a != 0,(b*c)/(a), 0),
#          wt.nonpara = 1/((a + 1) / (a + c + 2))) %>% summarise(NUE1 = sum(NUE1), NUE2 = sum(NUE2), MLE1 = sum(MLE1), MLE2 = sum(MLE2), NUEd = sum(NUEd), MLEd = sum(MLEd)) %>%
#   select(MASS.cat, age.cat, mod.vax.status, race.eth.collapsed, highADI, covid_treat, wt.nonpara)
# 
# enclave.nonpara.wt <- left_join(
#   enclave.eligible, nonpara.wts, by = c("MASS.cat", "age.cat", "mod.vax.status", "race.eth.collapsed", "highADI", "covid_treat") ) %>%
#   mutate(outpt.wts = if_else(outpatient.coviddx == 1 & covid_admit_death == 0, wt.nonpara, 1)) %>% # use weights for outpatient dx (Chapman NUE) to estimate denominator, no weights for numerator as all observed
#   mutate(mod.vax.status = fct_relevel(mod.vax.status, "Unvaccinated or fewer than 3 doses", "Vaccinated, last dose or COVID-19 ≥ 8 months prior"))
# 
# # projected cases
# est.cases.nonpara <- enclave.nonpara.wt  %>% filter(outpatient.coviddx == 1) %>%
#   summarise(n = sum(outpt.wts))
# 
# pct.nonpara.vs.para<- est.cases.nonpara/est.cases
# 
# rm(est.cases.nonpara, bootcases, denom.fit, bootfreq, inpatients)

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
