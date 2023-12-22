library(tidyverse)
library(marginaleffects)
library(boot)
library(flextable)
library(officer)
rm(list = ls())
load("enclave.eligible.wt.Rdata")

## NIH tier
nih.stats <- enclave.eligible.wt %>% filter(covid_treat == 0) %>%
  group_by(NIHtier) %>%
  summarize(covid_admit_death = sum(covid_admit_death*outpt.wts, na.rm = TRUE),
            death.4wks = sum(death.4wks*outpt.wts, na.rm = TRUE),
            covid_admit = sum(covid_admit*outpt.wts, na.rm = TRUE),
            total = sum(outpt.wts)) %>%
  ungroup() %>%
  mutate(death.pct = death.4wks/total,
         hosp.pct = covid_admit/total,
         severe.covid.pct = covid_admit_death/total) %>%
  mutate(NIHtier = case_when(
    NIHtier == "1" ~ "Tier 1", 
    NIHtier == "2" ~ "Tier 2", 
    NIHtier == "3" ~ "Tier 3" 
  )) %>% select(-total)

#obtain CI projected total
#svy object
nih.svy<- survey::svydesign(id = ~seqn, weights = ~outpt.wts, data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
                              mutate(NIHtier = case_when(
                                NIHtier == "1" ~ "Tier 1", 
                                NIHtier == "2" ~ "Tier 2", 
                                NIHtier == "3" ~ "Tier 3" )))

nih.stats <-
  left_join(nih.stats,
            enframe(survey::svytotal(~NIHtier, nih.svy)) %>%
              mutate(estimate = as.numeric(value), 
                     pct = estimate/sum(estimate),
                     lci = prettyNum(round(estimate - 1.96*survey::SE(value),0),  big.mark = ","),
                     uci = prettyNum(round(estimate + 1.96*survey::SE(value),0),  big.mark = ","),
                     estimate = paste0( prettyNum(round(estimate,0), big.mark = ","), " (", round(pct*100,0), "%)"),
                     total = paste0(estimate, "\n[", lci,  " to ", uci, "]"),
                     NIHtier = str_sub(name, 8)) %>%
              select(NIHtier, total), 
            by = "NIHtier")


## NIH Untreated


covid.risk.untreated.nih <- geepack::geeglm(
  covid_admit_death ~ NIHtier*(race.eth.collapsed + highADI),
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
    mutate(NIHtier = case_when(
    NIHtier == "1" ~ "Tier 1", 
    NIHtier == "2" ~ "Tier 2", 
    NIHtier == "3" ~ "Tier 3" 
  )))
summary(covid.risk.untreated.nih)


#obtaining est (95%CI) from g-estimation within strata of MASS.cat, vaccination, treatment (age, ADI, race/eth at dataset mean)
risks.untreated.nih <- 
    left_join(
      tibble(predictions(covid.risk.untreated.nih, 
                                       newdata = datagrid(NIHtier = c("Tier 1", "Tier 2", "Tier 3")))) %>%
              mutate(statistic = table1::signif_pad(100*estimate, digits = 2),
              conf.low = table1::signif_pad(100*conf.low, digits= 2),
              conf.high = table1::signif_pad(100*conf.high ,digits = 2),
              risk = paste0(statistic, "%\n[", conf.low, " to ", conf.high, "%]")) %>% select(NIHtier, risk),
    nih.stats, by = "NIHtier") %>% select(NIHtier, total, covid_admit, death.4wks, risk) 
## table creation

footer <- as_paragraph("Abbreviation: CI, confidence interval\n",
                       as_sup("a"), " Hospitalization within 14 days or death within 28 days of COVID-19 diagnosis.")
header <- "Table 3. Risk of severe COVID-19 by NIH prioritization tier for increased risk patients not accessing antiviral treatment— Mass General Brigham, June to December 2022"
tier1 <- "Individuals with moderate or severe immunocompromise, unvaccinated and 75 or older, or unvaccinated and 65 and or older with additional risk factors"
tier2 <- "Unvaccinated individuals not included in Tier 1 who are 65 or older or have clinical risk factors"
tier3 <- "Vaccinated individuals not included in Tier 1 who are 65 or older or have clinical risk factors"


flextable_nih <- 
  flextable(risks.untreated.nih) %>%
  # add footnotes
  mk_par(
    part = "header", j = 1, i = 1,
    value = as_paragraph(as_b("Risk group"))) %>%
  mk_par(
    part = "header", j = 2, i = 1,
    value = as_paragraph(as_b("Projected cases\n(95% CI)"))) %>%
  mk_par(
    part = "header", j = 3, i = 1,
    value = as_paragraph(as_b("Hospitalizations"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 4, i = 1,
    value = as_paragraph(as_b("Deaths"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 5, i = 1,
    value = as_paragraph(as_b("Covid Hosp or Death\n(95% CI)"))) %>%
  mk_par(
    part = "body", j = 1, i = 1,
    value = as_paragraph(as_b("Tier 1:\n"), tier1)) %>%
  mk_par(
    part = "body", j = 1, i = 2,
    value = as_paragraph(as_b("Tier 2:\n"), tier2)) %>%
  mk_par(
    part = "body", j = 1, i = 3,
    value = as_paragraph(as_b("Tier 3:\n"), tier3)) %>%
  bg(part = "body", bg = "#ebf4f1", i = c(2))%>%
  add_header_lines(header) %>%
  add_footer_lines(footer) %>%
  bold(i = 1, part = "header") %>%
  hline_top(part = "header",
            border = fp_border(color = "#037377",
                               width = 3,
                               style = "solid")) %>%
  hline(i = 1,
        part = "header",
        border = fp_border(color = "#037377",
                           width = 0.25,
                           style = "solid")) %>%
  hline_top(part = "body",
            border = fp_border(color = "#037377",
                               width = 0.25,
                               style = "solid")) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "#037377",
                                  width = 0.25,
                                  style = "solid")) %>%
  border_inner_h(part = "body",
                 border = fp_border(color = "#037377",
                                    width = 0.25,
                                    style = "dotted")) %>%
  autofit(part = "body") %>%
  align(part = "all", align = "center") %>%
  align(j = 1, part = "all", align = "left") %>%
  set_table_properties(layout = "autofit")

  



sect_properties <- prop_section(
  page_size = page_size(orient = "landscape", 
  width = 8.5, height = 11.7),
  page_margins = page_mar()
)

# Save as word .docx
save_as_docx(flextable_nih, path = "Table 3 NIH tier.docx",
             pr_section = sect_properties)


### risk by age

age.stats <- enclave.eligible.wt %>% filter(covid_treat == 0) %>%
  mutate(age.decade = case_when(
    age.at.covid < 30 ~ "18 to 29",
    age.at.covid < 50~ "30 to 49",
    age.at.covid < 65~ "50 to 64",
    age.at.covid < 75~ "65 to 74",
    age.at.covid < 85~ "75 to 84",
    TRUE ~ "85 and older"))  %>%
  group_by(age.decade) %>%
  summarize(covid_admit_death = sum(covid_admit_death*outpt.wts, na.rm = TRUE),
            death.4wks = sum(death.4wks*outpt.wts, na.rm = TRUE),
            covid_admit = sum(covid_admit*outpt.wts, na.rm = TRUE),
            total = sum(outpt.wts)) %>%
  ungroup() %>%
  mutate(death.pct = death.4wks/total,
         hosp.pct = covid_admit/total,
         severe.covid.pct = covid_admit_death/total) %>% select(-total)

#obtain CI projected total
#svy object
age.svy<- survey::svydesign(id = ~seqn, weights = ~outpt.wts, data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
                                mutate(age.decade = case_when(
                                  age.at.covid < 30 ~ "18 to 29",
                                  age.at.covid < 50~ "30 to 49",
                                  age.at.covid < 65~ "50 to 64",
                                  age.at.covid < 75~ "65 to 74",
                                  age.at.covid < 85~ "75 to 84",
                                  TRUE ~ "85 and older")) )

age.stats <-
  left_join(age.stats,
            enframe(survey::svytotal(~age.decade, age.svy)) %>%
              mutate(estimate = as.numeric(value), 
                     pct = estimate/sum(estimate),
                     lci = prettyNum(round(estimate - 1.96*survey::SE(value),0),  big.mark = ","),
                     uci = prettyNum(round(estimate + 1.96*survey::SE(value),0),  big.mark = ","),
                     estimate = paste0( prettyNum(round(estimate,0), big.mark = ","), " (", round(pct*100,0), "%)"),
                     total = paste0(estimate, "\n[", lci,  " to ", uci, "]"),
                     age.decade = str_sub(name,11)) %>%
              select(age.decade, total), 
            by = "age.decade")
## working to get percentages/formats to match NIH

covid.risk.untreated.age.decade <- geepack::geeglm(
  covid_admit_death ~ age.decade*(race.eth.collapsed + highADI + mod.vax.status),
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
    mutate(age.decade = case_when(
      age.at.covid < 30 ~ "18 to 29",
      age.at.covid < 50~ "30 to 49",
      age.at.covid < 65~ "50 to 64",
      age.at.covid < 75~ "65 to 74",
      age.at.covid < 85~ "75 to 84",
      TRUE ~ "85 and older")) )
summary(covid.risk.untreated.age.decade )


#obtaining est (95%CI) from g-estimation within strata of MASS.cat, vaccination, treatment (age, ADI, race/eth at dataset mean)
risks.untreated.age.decade <- 
  left_join(
    tibble(predictions(covid.risk.untreated.age.decade, 
                       newdata = datagrid(age.decade = c("85 and older", "75 to 84", "65 to 74", "50 to 64", "30 to 49", "18 to 29"   )))) %>%
      mutate(statistic = table1::signif_pad(100*estimate, digits = 2),
             conf.low = table1::signif_pad(100*conf.low, digits= 2),
             conf.high = table1::signif_pad(100*conf.high ,digits = 2),
             risk = paste0(statistic, "%\n[", conf.low, "% to ", conf.high, "%]")) %>% select(age.decade , risk),
    age.stats, by = "age.decade") %>% select(age.decade , total, covid_admit, death.4wks, risk)
## table creation


footer <- as_paragraph("Abbreviation: CI, confidence interval\n",
                       as_sup("a"), " Hospitalization within 14 days or death within 28 days of COVID-19 diagnosis.")
header <- "Table 5. Risk of severe COVID-19 by age group for increased risk patients not accessing antiviral treatment— Mass General Brigham, June to December 2022"

flextable_age <- 
  flextable(risks.untreated.age.decade) %>%
  # add footnotes
  mk_par(
    part = "header", j = 1, i = 1,
    value = as_paragraph(as_b("Age group"))) %>%
  mk_par(
    part = "header", j = 2, i = 1,
    value = as_paragraph(as_b("Projected cases\n(95% CI)"))) %>%
  mk_par(
    part = "header", j = 3, i = 1,
    value = as_paragraph(as_b("Hospitalizations"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 4, i = 1,
    value = as_paragraph(as_b("Deaths"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 5, i = 1,
    value = as_paragraph(as_b("Covid Hosp or Death\n[95% CI]"))) %>%
  bg(part = "body", bg = "#ebf4f1", i = c(2,4))%>%
  add_header_lines(header) %>%
  add_footer_lines(footer) %>%
  bold(i = 1, part = "header") %>%
  hline_top(part = "header",
            border = fp_border(color = "#037377",
                               width = 3,
                               style = "solid")) %>%
  hline(i = 1,
        part = "header",
        border = fp_border(color = "#037377",
                           width = 0.25,
                           style = "solid")) %>%
  hline_top(part = "body",
            border = fp_border(color = "#037377",
                               width = 0.25,
                               style = "solid")) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "#037377",
                                  width = 0.25,
                                  style = "solid")) %>%
  border_inner_h(part = "body",
                 border = fp_border(color = "#037377",
                                    width = 0.25,
                                    style = "dotted")) %>%
  autofit(part = "body") %>%
  align(part = "all", align = "center") %>%
  align(j = 1, part = "all", align = "left") %>%
  set_table_properties(layout = "autofit")


sect_properties <- prop_section(
  page_size = page_size(orient = "landscape", 
                        width = 8.5, height = 11.7),
  page_margins = page_mar()
)

# Save as word .docx
save_as_docx(flextable_age, path = "Table 5 Age group.docx",
             pr_section = sect_properties)



### WHO groupings
who.stats <- enclave.eligible.wt %>% filter(covid_treat == 0) %>%
  group_by(WHOrisk) %>%
  summarize(covid_admit_death = sum(covid_admit_death*outpt.wts, na.rm = TRUE),
            death.4wks = sum(death.4wks*outpt.wts, na.rm = TRUE),
            covid_admit = sum(covid_admit*outpt.wts, na.rm = TRUE),
            total = sum(outpt.wts)) %>%
  ungroup() %>%
  mutate(death.pct = death.4wks/total,
         hosp.pct = covid_admit/total,
         severe.covid.pct = covid_admit_death/total) %>%
  mutate(WHOrisk = case_when(
    WHOrisk == "1" ~ "High risk", 
    WHOrisk == "2" ~ "Moderate risk", 
    WHOrisk == "3" ~ "Low risk" 
  )) %>% select(-total)

#obtain CI projected total
#svy object
who.svy<- survey::svydesign(id = ~seqn, weights = ~outpt.wts, data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
                              mutate(WHOrisk = case_when(
                                WHOrisk == "1" ~ "High risk", 
                                WHOrisk == "2" ~ "Moderate risk", 
                                WHOrisk == "3" ~ "Low risk" )))

who.stats <-
  left_join(who.stats,
enframe(survey::svytotal(~WHOrisk, who.svy)) %>%
  mutate(estimate = as.numeric(value), 
         pct = estimate/sum(estimate),
         lci = prettyNum(round(estimate - 1.96*survey::SE(value),0),  big.mark = ","),
         uci = prettyNum(round(estimate + 1.96*survey::SE(value),0),  big.mark = ","),
         estimate = paste0( prettyNum(round(estimate,0), big.mark = ","), " (", round(pct*100,0), "%)"),
         total = paste0(estimate, "\n[", lci,  " to ", uci, "]"),
         WHOrisk = str_sub(name, 8)) %>%
  select(WHOrisk, total), 
by = "WHOrisk")


## WHO Untreated


covid.risk.untreated.who <- geepack::geeglm(
  covid_admit_death ~ WHOrisk*(race.eth.collapsed + highADI + mod.vax.status),
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
    mutate(WHOrisk = case_when(
      WHOrisk == "1" ~ "High risk", 
      WHOrisk == "2" ~ "Moderate risk", 
      WHOrisk == "3" ~ "Low risk" 
    )))
summary(covid.risk.untreated.who)
tidy(covid.risk.untreated.who)



#obtaining est (95%CI) from g-estimation within strata of MASS.cat, vaccination, treatment (age, ADI, race/eth at dataset mean)
risks.untreated.who <- 
  left_join(
    tibble(predictions(covid.risk.untreated.who, 
                       newdata = datagrid(WHOrisk = c("High risk", "Moderate risk",  "Low risk")))) %>%
      mutate(statistic = table1::signif_pad(100*estimate, digits = 2),
             conf.low = table1::signif_pad(100*conf.low, digits= 2),
             conf.high = table1::signif_pad(100*conf.high ,digits = 2),
             risk = paste0(statistic, "%\n(", conf.low, "% to ", conf.high, "%)")) %>% select(WHOrisk, risk),
    who.stats, by = "WHOrisk") %>% select(WHOrisk, total, covid_admit, death.4wks, risk) 

## table creation


footer <- as_paragraph("Abbreviation: CI, confidence interval\n",
                       as_sup("a"), " Hospitalization within 14 days or death within 28 days of COVID-19 diagnosis.\n",
                       as_sup("b"), " Estimates only includes individuals with at least one increased risk condition. Consequently, estimated risk likely higher than if individuals without increased risk conditions were included")
header <- "Table 4. Risk of severe COVID-19 by WHO risk groups for increased risk patients not accessing antiviral treatment— Mass General Brigham, June to December 2022"

who1 <- "People who have been diagnosed with immunodeficiency syndromes, undergone a solid organ transplant and receiving immunosuppressants, or have autoimmune disease and are receiving immunosuppressants"
who2 <- "People over 65 years old or with certain health conditions: obesity, diabetes, chronic cardiopulmonary disease, chronic kidney disease, chronic liver disease, active cancer, disability, and comorbidities of chronic disease"
who3 <- "People who are not high or moderate risk"


flextable_who <- 
  flextable(risks.untreated.who) %>%
  # add footnotes
  mk_par(
    part = "header", j = 1, i = 1,
    value = as_paragraph(as_b("WHO risk group"))) %>%
  mk_par(
    part = "header", j = 2, i = 1,
    value = as_paragraph(as_b("Projected cases\n(95% CI)"))) %>%
  mk_par(
    part = "header", j = 3, i = 1,
    value = as_paragraph(as_b("Hospitalizations"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 4, i = 1,
    value = as_paragraph(as_b("Deaths"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 5, i = 1,
    value = as_paragraph(as_b("Covid Hosp or Death\n(95% CI)"))) %>%
  mk_par(
    part = "body", j = 1, i = 1,
    value = as_paragraph(as_b("High risk:\n"), who1)) %>%
  mk_par(
    part = "body", j = 1, i = 2,
    value = as_paragraph(as_b("Moderate risk:\n"), who2)) %>%
  mk_par(
    part = "body", j = 1, i = 3,
    value = as_paragraph(as_b("Low risk:\n"), who3)) %>%
  mk_par(
    part = "body", j = 5, i = 3,
    value = as_paragraph(as.character(risks.untreated.who[3,5]), as_sup("b"))) %>%
  bg(part = "body", bg = "#ebf4f1", i = c(2))%>%
  add_header_lines(header) %>%
  add_footer_lines(footer) %>%
  bold(i = 1, part = "header") %>%
  hline_top(part = "header",
            border = fp_border(color = "#037377",
                               width = 3,
                               style = "solid")) %>%
  hline(i = 1,
        part = "header",
        border = fp_border(color = "#037377",
                           width = 0.25,
                           style = "solid")) %>%
  hline_top(part = "body",
            border = fp_border(color = "#037377",
                               width = 0.25,
                               style = "solid")) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "#037377",
                                  width = 0.25,
                                  style = "solid")) %>%
  border_inner_h(part = "body",
                 border = fp_border(color = "#037377",
                                    width = 0.25,
                                    style = "dotted")) %>%
  autofit(part = "body") %>%
  align(part = "all", align = "center") %>%
  align(j = 1, part = "all", align = "left") %>%
  set_table_properties(layout = "autofit")


sect_properties <- prop_section(
  page_size = page_size(orient = "landscape", 
                        width = 8.5, height = 11.7),
  page_margins = page_mar()
)

# Save as word .docx
save_as_docx(flextable_who, path = "Table 4 WHO groupings.docx",
             pr_section = sect_properties)


## assess mod vs severe immunocompromise

covid.risk.untreated.ic <- geepack::geeglm(
  covid_admit_death ~ IC.type*(race.eth.collapsed + highADI + mod.vax.status + age.cat),
  family = poisson(link=log),
  weights = outpt.wts,
  id = seqn,
  corstr = "independence",
  data = enclave.eligible.wt %>% filter(covid_treat == 0) %>%
    mutate(IC.type = case_when(
      severe.immunosuppression == "1" ~ "Severe immunocompromise", 
      MASS.ic == "3" ~ "Moderate immunocompromise",
      TRUE ~ "No recorded immunocompromise"
    )) )
summary(covid.risk.untreated.ic)


#obtaining est (95%CI) from g-estimation within strata of MASS.cat, vaccination, treatment (age, ADI, race/eth at dataset mean)
risks.untreated.ic <- 
    tibble(predictions(covid.risk.untreated.ic, 
                       newdata = datagrid(IC.type = c("Severe immunocompromise", "Moderate immunocompromise", "No recorded immunocompromise")))) %>%
      mutate(statistic = table1::signif_pad(100*estimate, digits = 2),
             conf.low = table1::signif_pad(100*conf.low, digits= 2),
             conf.high = table1::signif_pad(100*conf.high ,digits = 2),
             risk = paste0(statistic, "%\n[", conf.low, " to ", conf.high, "%]")) %>% select(IC.type, risk)
risks.untreated.ic


