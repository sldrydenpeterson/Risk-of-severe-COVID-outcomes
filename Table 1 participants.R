rm(list = ls())
library(tidyverse)
library(janitor)
library(flextable)
library(officer)
library(gtsummary)
library(survey)
library(table1)

fixzip <- function(x){ifelse(nchar(x)<5, paste0(0,x), x)}

load("enclave.eligible.wt.Rdata")

enclave.eligible.tab1 <- enclave.eligible.wt %>%
  mutate(
    outpatient.coviddx = fct_relevel(as.factor(outpatient.coviddx), "0"),
    death.4wks = fct_relevel(as.factor(death.4wks), "0"),
    covid_admit = fct_relevel(as.factor(covid_admit), "0"),
    severe.immunosuppression = fct_relevel(as.factor(severe.immunosuppression ), "0"),
    vaxstatus = as.factor(vaxstatus),
    vaxstatus = fct_relevel(vaxstatus, "Vaccinated and boosted", "Vaccinated"),
    lastvaxgt8mo  = as.factor(lastvaxgt8mo),
    race.eth =as.factor(race.eth),
    highADI = fct_relevel(as.factor(highADI), "0"),
    MASS.ic = as.factor(if_else(MASS.ic == 3, 1, 0)),
    CVD = as.factor(CVD),
    diabetes= as.factor(if_else(MASS.dm == 2, 1, 0)),
    pulm.dis = as.factor(pulm.dis),
    inflam.dis = as.factor(inflam.dis),
    age.cat = as.factor(age.cat),
    sex = as.factor(sex),
    solid.malig = as.factor(solid.malig),
    heme.malig = as.factor(heme.malig),
    MASS.cat = as.factor(MASS.cat),
    psych.non.unipolar = as.factor(psych.non.unipolar),
    depr.anxiety = as.factor(depr.anxiety),
    covid_treat = as.factor(covid_treat),
    seqid = row_number()) %>% 
  select(seqid, age.cat, sex, race.eth, highADI, mod.vax.status,   MASS.cat,   
         obesity.cat, MASS.ic, severe.immunosuppression, diabetes, CVD, pulm.dis, psych.non.unipolar, depr.anxiety,
         heme.malig, solid.malig, inflam.dis, covid_treat, covid_admit, death.4wks, outpt.wts, covid_admit_death, outpatient.coviddx) 



#unweighted (observed)
tab1.raw<- tbl_summary(enclave.eligible.tab1 %>% select(-c(covid_admit_death, covid_admit, 
                                                           death.4wks, outpt.wts, seqid, outpatient.coviddx)), by = NULL,
                       type = list(c(highADI,  MASS.ic, diabetes, CVD, pulm.dis, inflam.dis, solid.malig, heme.malig,
                                     depr.anxiety, psych.non.unipolar, severe.immunosuppression, covid_treat)~  "dichotomous"),
                       label = list( mod.vax.status ~ "Vaccination status", 
                                     covid_treat ~ "Received outpatient COVID antiviral therapy",
                                     sex ~ "Sex",
                                     highADI ~ "Increased neighborhood disadvantage (ADI)",
                                     race.eth ~ "Race and ethnicity",
                                     obesity.cat ~ "Body mass index (BMI)", 
                                     age.cat ~ "Age", 
                                     MASS.cat ~"Comorbidity score", 
                                     MASS.ic ~ "Immunocompromise", 
                                     diabetes ~ "Diabetes", 
                                     CVD ~ "Heart disease or stroke",
                                     pulm.dis ~"Pulmonary disease",
                                     inflam.dis ~ "Rheumatologic or inflammatory bowel disease", 
                                     solid.malig ~"Solid tumor malignancy", 
                                     heme.malig ~ "Hematologic malignancy", 
                                     depr.anxiety ~ "Depression and anxiety", 
                                     psych.non.unipolar ~ "Bipolar, schizophrenia, and other disorders", 
                                     severe.immunosuppression ~ "Severe immunocompromise"))


total.raw<- tab1.raw$table_styling$header$modify_stat_N[1]

#create tibble of gt table to format effectively
tab1.raw<- as_tibble(tab1.raw) 

#rename columns
colnames(tab1.raw)[1] <- "Characteristic"
colnames(tab1.raw)[2] <- "raw.counts"

tab1.raw<- tab1.raw %>%
  mutate(row.n = str_split(raw.counts, " ")) %>%
  unnest_wider(row.n, names_sep = "-") %>%
  mutate(row.n.raw = as.numeric(gsub(",", "", `row.n-1`))) %>%
  select(Characteristic, raw.counts, row.n.raw) %>%
  add_row(.before = 1,
          Characteristic = "Total",
          row.n.raw = as.numeric(total.raw),
          raw.counts = paste0(prettyNum(total.raw, big.mark = ","), " (100%)"))





## weighted (projected)
enclave.tab1.svy <- svydesign(id=~seqid,  weights=~outpt.wts, 
                              data= enclave.eligible.tab1 %>%     
                                select(-c(outpatient.coviddx, covid_admit, death.4wks, covid_admit_death)) )

tab1.wt <- tbl_svysummary(enclave.tab1.svy, by = NULL,
                     type = list(c(highADI,  MASS.ic, diabetes, CVD, pulm.dis, inflam.dis, solid.malig, heme.malig,
                                   depr.anxiety, psych.non.unipolar, severe.immunosuppression, covid_treat) ~  "dichotomous"),
                     label = list( mod.vax.status ~ "Vaccination status", 
                                   covid_treat ~ "Received outpatient COVID antiviral therapy",
                                   sex ~ "Sex",
                                   highADI ~ "Increased neighborhood disadvantage (ADI)",
                                   race.eth ~ "Race and ethnicity",
                                   obesity.cat ~ "Body mass index (BMI)", 
                                   age.cat ~ "Age", 
                                   MASS.cat ~"Comorbidity score", 
                                   MASS.ic ~ "Immunocompromise", 
                                   diabetes ~ "Diabetes", 
                                   CVD ~ "Heart disease or stroke",
                                   pulm.dis ~"Pulmonary disease",
                                   inflam.dis ~ "Rheumatologic or inflammatory bowel disease", 
                                   solid.malig ~"Solid tumor malignancy", 
                                   heme.malig ~ "Hematologic malignancy", 
                                   depr.anxiety ~ "Depression and anxiety", 
                                   psych.non.unipolar ~ "Bipolar, schizophrenia, and other disorders", 
                                   severe.immunosuppression ~ "Severe immunocompromise")) %>%
  add_ci(statistic = list(all_categorical() ~ "{conf.low},{conf.high}"),
         style_fun = all_categorical() ~ purrr::partial(style_sigfig, digits = 9)) 

total.wt<- tab1.wt$table_styling$header$modify_stat_N[1]

N.total.wt<- svytotal(~marker, svydesign(id=~seqid,  weights=~outpt.wts, 
                                         data= enclave.eligible.tab1 %>% mutate(marker =1)))

#create tibble of gt table to format effectively
tab1.wt<- as_tibble(tab1.wt) 

#rename columns
colnames(tab1.wt)[1] <- "Characteristic"
colnames(tab1.wt)[2] <- "counts"
colnames(tab1.wt)[3] <- "proportions"

tab1.wt <- tab1.wt %>%
  filter(Characteristic !="seqid" & Characteristic != "outpt.wts") %>%
  mutate(p_weight = str_split(proportions, ",")) %>%
  unnest_wider(p_weight, names_sep = "-") %>%
  rename(pct_lci = `p_weight-1`,
         pct_uci = `p_weight-2`,
         p_counts = counts) %>%
  mutate(pct_lci = as.numeric(pct_lci),
         pct_uci = as.numeric(pct_uci)) %>%
  mutate(lci = if_else(!is.na(pct_lci), round(pct_lci*total.wt, 0), as.numeric(NA)),
         uci = if_else(!is.na(pct_uci), round(pct_uci*total.wt, 0), as.numeric(NA)),
         uci = prettyNum(uci, big.mark = ",", scientific=FALSE),
         lci = prettyNum(lci, big.mark = ",", scientific=FALSE),
         N.wt = if_else(!is.na(p_counts), paste0(p_counts, " [", lci, " to ", uci, "]" ), as.character(NA)),
         row.n = str_split(p_counts, " ")) %>%
  unnest_wider(row.n, names_sep = "-") %>%
  mutate(row.n = as.numeric(gsub(",", "", `row.n-1`))) %>%
  select(Characteristic, N.wt, p_counts, row.n ) %>%
  add_row(.before = 1, 
          Characteristic = "Total", 
          row.n = as.numeric(N.total.wt),
          N.wt = paste0(prettyNum(as.numeric(N.total.wt), big.mark = ","), " (100%) [", 
                        prettyNum(as.numeric(N.total.wt)-1.96*as.numeric(SE(N.total.wt)), big.mark = ","), " to ",
                        prettyNum(as.numeric(N.total.wt)+1.96*as.numeric(SE(N.total.wt)), big.mark = ","), "]")
  )




#obtain strata counts of severe covid
tab1.raw.severe<- as_tibble(tbl_summary(enclave.eligible.tab1 %>% 
            filter(covid_admit_death == 1) %>%
            select(-c(covid_admit_death, covid_admit, 
                      death.4wks, outpt.wts, seqid, outpatient.coviddx)), by = NULL,
          type = list(c(highADI,  MASS.ic, diabetes, CVD, pulm.dis, inflam.dis, solid.malig, heme.malig,
                        depr.anxiety, psych.non.unipolar, severe.immunosuppression, covid_treat)~  "dichotomous"),
          label = list( mod.vax.status ~ "Vaccination status", 
                        covid_treat ~ "Received outpatient COVID antiviral therapy",
                        sex ~ "Sex",
                        highADI ~ "Increased neighborhood disadvantage (ADI)",
                        race.eth ~ "Race and ethnicity",
                        obesity.cat ~ "Body mass index (BMI)", 
                        age.cat ~ "Age", 
                        MASS.cat ~"Comorbidity score", 
                        MASS.ic ~ "Immunocompromise", 
                        diabetes ~ "Diabetes", 
                        CVD ~ "Heart disease or stroke",
                        pulm.dis ~"Pulmonary disease",
                        inflam.dis ~ "Rheumatologic or inflammatory bowel disease", 
                        solid.malig ~"Solid tumor malignancy", 
                        heme.malig ~ "Hematologic malignancy", 
                        depr.anxiety ~ "Depression and anxiety", 
                        psych.non.unipolar ~ "Bipolar, schizophrenia, and other disorders", 
                        severe.immunosuppression ~ "Severe immunocompromise")))
#rename columns
colnames(tab1.raw.severe)[1] <- "Characteristic"
colnames(tab1.raw.severe)[2] <- "severe"

tab1.raw.severe<-
  tab1.raw.severe %>%
  mutate(severe = str_split(severe, " ")) %>%
  unnest_wider(severe, names_sep = "-") %>%
  rename(severe = `severe-1`) %>%
  select(Characteristic, severe) %>%
  add_row(.before = 1, Characteristic = "Total", 
          severe = prettyNum(sum(enclave.eligible.tab1$covid_admit_death), big.mark = ","))




#######################
#constuct final dataframe

dichotomous<- c("Immunocompromise",
                "Severe immunocompromise",
                "Diabetes",
                "Heart disease or stroke",
                "Pulmonary disease",
                "Bipolar, schizophrenia, and other disorders",
                "Depression and anxiety",
                "Hematologic malignancy",
                "Solid tumor malignancy",
                "Rheumatologic or inflammatory bowel disease",
                "Received outpatient COVID antiviral therapy",
                "Total")

tab1 <-  cbind(tab1.raw, tab1.wt %>% select(-Characteristic), tab1.raw.severe %>% select(-Characteristic)) %>%
  mutate(ratio = table1::signif_pad(row.n / row.n.raw, digits = 3 ),
         Characteristic = case_when(
           is.na(raw.counts) ~ Characteristic,
           Characteristic %in% dichotomous ~ Characteristic,
           !is.na(raw.counts) ~ paste0("   ", Characteristic)
         )) %>%
  select(Characteristic, raw.counts, N.wt, ratio, severe)

tab1[25,1] <- paste0("   ", "Severe immmunocompromise") #include needed spacing, removed above as same name later in table





#######################
#table formatting
tab1_defaults <- function(){
  set_flextable_defaults(font.family = "Calibri",
                         font.size = 10,
                         border.color = "black")
}
custom_tab <- function(df, header, footer){
  flextable(df) %>%
    theme_zebra(
      odd_header = "transparent",
      odd_body = "#ebf4f1",
      even_header = "transparent",
      even_body = "white"
    )%>%
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
   #bg(part = "body", bg = "#ebf4f1") %>%
    align(part = "all", align = "center") %>%
    align(j = 1, part = "all", align = "left") 
}
# Set Table header

header <- str_squish("Table 1. Characteristics of increased risk individuals with COVID-19— recorded and estimated total cases\nMass General Brigham, June to December 2022")
footer <- as_paragraph("Abbreviations: ADI, Area Deprivation Index, CI, confidence interval\n",
                       as_sup("a"), " Cases estimated through accounting for selection bias using the strata-specific, inverse probability of a preceding recorded positive COVID-19 test prior to the onset of a severe outcome as weights.\n")

# Set custom_tab() defaults
tab1_defaults()

# Create the flextable object
flextable_1 <- custom_tab(tab1, header, footer) %>%
  mk_par(
    part = "header", j = 1, i = 2,
    value = as_paragraph(as_b("Characteristic"))) %>%
  mk_par(
    part = "header", j = 2, i = 2,
    value = as_paragraph(as_b("Recorded cases"))) %>%
  mk_par(
    part = "header", j = 3, i = 2,
    value = as_paragraph(as_b("Estimated total cases"), as_sup("a"), "\n[95% CI]")) %>%
  mk_par(
    part = "header", j = 4, i = 2,
    value = as_paragraph(as_b("Ratio of estimated\ntotal to recorded\ncases"))) %>%
  mk_par(
    part = "header", j = 5, i = 2,
    value = as_paragraph(as_b("Severe outcomes\nof COVID-19"))) %>%
  set_table_properties(layout = "autofit")

# Save as word .docx
save_as_docx(flextable_1, path = "Table 1 participants.docx",
             pr_section =
               prop_section(page_size = page_size(orient = "landscape"),
                            type = "continuous"))


rm(list = ls())