library(tidyverse)
library(lubridate)
library(broom)
library(boot)

# NNT

load("enclave.eligible.wt.Rdata")
## propensity of treatment within strata, by age, race/eth, ADI

# overall population to assess bias in prescribing patterns


denom.fit <- glm(
  covid_treat ~ as.factor(MASS.cat) + as.factor(age.cat) 
  + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed) 
  + as.factor(highADI),
  family = binomial(),
  data = enclave.eligible.wt )

tidy(denom.fit, conf.int = TRUE, exponentiate = TRUE)


## logistic models to generate denominator if IP weights, probability of treatment within MASS/IC strata

# calculate weights for treatment and then obtain product
#MASS 3 or less
MASS3orless <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 3 or less")
    denom.fit <- glm(
      covid_treat ~ as.factor(age.cat) 
      + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed) 
      + as.factor(highADI),
      family = binomial(),
      data = MASS3orless)
    
    tidy(denom.fit, conf.int = TRUE, exponentiate = TRUE)
    
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


    

##MASS 4 or 5
MASS4or5 <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 4 and 5" ) 
    denom.fit <- glm(
      covid_treat ~ as.factor(age.cat) 
      + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed) 
      + as.factor(highADI),
      family = binomial(),
      data = MASS4or5)
    
    tidy(denom.fit, conf.int = TRUE, exponentiate = TRUE)
    
    # denominator of weights
    pd.treat <- predict(denom.fit, type = "response")
    
    # numerator of stabilized weights, overall probabilty of outpatient treatment
    numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS4or5)
    pn.treat <- predict(numer.fit, type = "response")
    
    ## calculate weights
    MASS4or5$sw.treat <-if_else(MASS4or5$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                                   (pn.treat / pd.treat))
    
    MASS4or5$rd.weight <- MASS4or5$sw.treat * MASS4or5$outpt.wts

    
    quantile(MASS4or5$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
    mean(MASS4or5$outpt.wts)
    
    quantile(MASS4or5$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
    mean(MASS4or5$sw.treat)
    
    quantile(MASS4or5$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
    mean(MASS4or5$rd.weight)
    

##MASS 6 or greater
MASS6ormore <- enclave.eligible.wt %>% filter(MASS.cat == "MASS 6 or greater" ) 
    denom.fit <- glm(
      covid_treat ~ as.factor(age.cat) 
      + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed) 
      + as.factor(highADI),
      family = binomial(),
      data = MASS6ormore)
    
    tidy(denom.fit, conf.int = TRUE, exponentiate = TRUE)
    
    # denominator of weights
    pd.treat <- predict(denom.fit, type = "response")
    
    # numerator of stabilized weights, overall probabilty of outpatient treatment
    numer.fit <- glm(covid_treat ~ 1, family = binomial(),  data = MASS6ormore)
    pn.treat <- predict(numer.fit, type = "response")
    
    ## calculate weights
    MASS6ormore$sw.treat <-if_else(MASS6ormore$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                                (pn.treat / pd.treat))
    MASS6ormore$rd.weight <- MASS6ormore$sw.treat * MASS6ormore$outpt.wts
    
    
    quantile(MASS6ormore$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
    mean(MASS6ormore$outpt.wts)
    
    quantile(MASS6ormore$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
    mean(MASS6ormore$sw.treat)
    
    quantile(MASS6ormore$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
    mean(MASS6ormore$rd.weight)

##severe immunocompromise
severeIC <- enclave.eligible.wt %>% filter(MASS.cat == "Severe immunocompromise" ) 
    denom.fit <- glm(
      covid_treat ~ as.factor(age.cat) 
      + as.factor(mod.vax.status)  +  as.factor(race.eth.collapsed) 
      + as.factor(highADI),
      family = quasibinomial(),
      weights = outpt.wts,
      data = severeIC)
    
    tidy(denom.fit, conf.int = TRUE, exponentiate = TRUE)
    
    # denominator of weights
    pd.treat <- predict(denom.fit, type = "response")
    
    # numerator of stabilized weights, overall probabilty of outpatient treatment
    numer.fit <- glm(covid_treat ~ 1, family = quasibinomial(), weights = outpt.wts,  data = severeIC)
    pn.treat <- predict(numer.fit, type = "response")
    
    ## calculate weights
    severeIC$sw.treat <-if_else(severeIC$covid_treat == 0, ((1- pn.treat) / (1-pd.treat)),
                                   (pn.treat / pd.treat))
    severeIC$rd.weight <- severeIC$sw.treat * severeIC$outpt.wts
    
    
    quantile(severeIC$outpt.wts, c(0, 0.025, 0.5, 0.975, 1))
    mean(severeIC$outpt.wts)
    
    quantile(severeIC$sw.treat, c(0, 0.025, 0.5, 0.975, 1))
    mean(severeIC$sw.treat)
    
    quantile(severeIC$rd.weight, c(0, 0.025, 0.5, 0.975, 1))
    mean(severeIC$rd.weight)
    
## re-combine dataset
    treatment.NNT.wt <- rbind(MASS3orless, MASS4or5, MASS6ormore, severeIC)

## point estimate
    
    # descriptive statistics, overall
    treatment.NNT.wt %>%
      group_by(MASS.cat, covid_treat) %>%
      summarize(covid_admit_death = sum(covid_admit_death),
                death.4wks = sum(death.4wks),
                covid_admit = sum(covid_admit),
                total = sum(if_else(outpatient.coviddx == 1, rd.weight, 0))) %>%
      ungroup() %>% 
      mutate(death.pct = death.4wks/total,
             hosp.pct = covid_admit/total,
             severe.covid.pct = covid_admit_death/total) %>%
      select(MASS.cat, covid_treat, severe.covid.pct) %>%
      pivot_wider(names_from = covid_treat, values_from = severe.covid.pct) %>%
      rename(treat = `1`, untreat = `0`) %>%
      mutate(rd = untreat - treat,
             NNT = 1/rd)
    
    


# bootstrap 95% CI
    # 1/(
    #   #calculate risk difference
    #   ## risk of untreated
    #   (
    #     sum( (subset(treatment.NNT.wt, covid_treat == 0))$covid_admit_death, na.rm = TRUE) /
    #       sum(if_else((subset(treatment.NNT.wt, covid_treat == 0))$outpatient.coviddx == 1, (subset(treatment.NNT.wt, covid_treat == 0))$rd.weight, 0)) 
    #   ) -
    #   ## subrct risk of treated  
    #   (
    #     sum( (subset(treatment.NNT.wt, covid_treat == 1))$covid_admit_death, na.rm = TRUE) /
    #       sum(if_else((subset(treatment.NNT.wt, covid_treat == 1))$outpatient.coviddx == 1, (subset(treatment.NNT.wt, covid_treat == 1))$rd.weight, 0)) 
    #   )
    # )
    
    # create function to get NNT
    boot_NNT <- function(d, i) {
      d2 <- d[i,]
      return(
        1/(
        #calculate risk difference
        ## risk of untreated
          (
        sum( (subset(d2, covid_treat == 0))$covid_admit_death, na.rm = TRUE) /
        sum(if_else((subset(d2, covid_treat == 0))$outpatient.coviddx == 1, (subset(d2, covid_treat == 0))$rd.weight, 0)) 
        ) -
        ## subrct risk of treated  
          (
          sum( (subset(d2, covid_treat == 1))$covid_admit_death, na.rm = TRUE) /
          sum(if_else((subset(d2, covid_treat == 1))$outpatient.coviddx == 1, (subset(d2, covid_treat == 1))$rd.weight, 0)) 
          ))
      )
    }
    
    set.seed(4454) # reproducible results
    
    # MASS 3 or less
    bootfreq <- boot(treatment.NNT.wt %>%
                       filter(MASS.cat == "MASS 3 or less" ), 
                     boot_NNT, R=2099)  
    MASS3orless.nnt <- tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
      mutate(measure = "Number Needed to Treat", MASS.cat = "MASS 3 or less", mod.vax.status = "Vaccinated, last booster or COVID < 8 months prior") %>%
      select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
    
    
    # MASS 4 and 5
    bootfreq <- boot(treatment.NNT.wt %>%
                       filter(MASS.cat == "MASS 4 and 5" ), 
                     boot_NNT, R=2099)  
    MASS4and5.nnt<-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
      mutate(measure = "Number Needed to Treat", MASS.cat = "MASS 4 and 5", mod.vax.status = "Vaccinated, last booster or COVID < 8 months prior") %>%
      select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
    
    
    # MASS 6 or greater
    bootfreq <- boot(treatment.NNT.wt %>%
                       filter(MASS.cat == "MASS 6 or greater" ), 
                     boot_NNT, R=2099)  
    MASS6orgreater.nnt <-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
      mutate(measure = "Number Needed to Treat", MASS.cat = "MASS 6 or greater", mod.vax.status = "Vaccinated, last booster or COVID < 8 months prior") %>%
      select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
    
    # SevereIC 
    bootfreq <- boot(treatment.NNT.wt %>%
                       filter(MASS.cat == "Severe immunocompromise" ), 
                     boot_NNT, R=2099)  
    severeIC.nnt <-tidy(bootfreq, conf.int = TRUE,  conf.level = .95) %>% 
      mutate(measure = "Number Needed to Treat", MASS.cat = "Severe immunocompromise", mod.vax.status = "Vaccinated, last booster or COVID < 8 months prior") %>%
      select(measure, MASS.cat, mod.vax.status, statistic, conf.low, conf.high)
    
    
    nnt.all  <- rbind( MASS3orless.nnt, MASS4and5.nnt, MASS6orgreater.nnt, severeIC.nnt) %>%
      #risk difference and 95%CI
      mutate(rd = round(100/-statistic, 1),
              rd.high = round(100/(-conf.low), 1),
              rd.low = round(100/(-conf.high), 1),
              rd = paste0(rd, "% (", rd.low, "% to ", rd.high, "%)")
              
      ) %>%
      #NNT and 95%CI
      mutate(statistic = round(statistic, 0),
             conf.low = round(conf.low , 0),
             conf.high = round(conf.high , 0),
             nnt = paste0(statistic, " (", conf.low, " to ", conf.high, ")")) %>%
      # combine
      mutate(rd.nnt = paste0("RD ", rd, "\nNNTB ", nnt )) %>%
      select(MASS.cat, rd.nnt) %>%
      rename("MASS score" = MASS.cat) 
    
save(nnt.all, file = "nnt.all.Rdata")
    
    rm(list = ls())
    