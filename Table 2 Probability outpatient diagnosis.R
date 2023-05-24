library(flextable)
library(officer)


# Set Table header/footer
header <- str_squish(str_remove("Table 2. Odds ratios of preceding outpatient diagnosis of COVID-19* among patients with severe outcomes by selected characteristics— Mass General Brigham, June to December 2022", "\n"))
footer <- paste0("* Recorded positive COVID-19 test result two or more calendar days prior to COVID hospitalization or death.\n", 
                 "† Adjusted odds ratios from multivariable logistic model of apriori-selected characteristics are shown. In subsequent univariate analysis, all characteristics were significantly associated with preceding outpatient diagnosis (all p < 0.01), except neighborhood disadvantage (p = 0.8).")

denom.fit.rename <- glm(
  outpatient.coviddx ~ as.factor(MASS.cat) + as.factor(age.cat) + as.factor(mod.vax.status)  
  + as.factor(race.eth.collapsed)  + as.factor(highADI) + as.factor(covid_treat) ,
  family = binomial(),
  data = inpatients %>% 
    mutate(highADI = if_else(highADI == 0, "1st to 74th percentile", "75th or higher percentile"),
                               covid_treat = if_else(covid_treat == 0, "No treatment", "Received treatment")))

predict.outpatientdx <- gtsummary::tbl_regression(denom.fit.rename, exponentiate = TRUE, label = list(`as.factor(MASS.cat)` ~ "Comorbidity score", `as.factor(age.cat)` ~ "Age",
                              `as.factor(mod.vax.status)` ~ "Vaccination status", `as.factor(race.eth.collapsed)` ~ "Race and Ethnicity",
                              `as.factor(highADI)` ~ "Neighborhood disadvantage (ADI)", `as.factor(covid_treat)` ~ "Outpatient COVID-19 treatment")) %>%
  gtsummary::add_global_p() %>% # type 3
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>% # remove footnotes
  gtsummary::modify_header(estimate = "**Odds Ratio†**") %>%
  gtsummary::as_flex_table() %>%
  add_header_lines(header) %>%
  add_footer_lines(footer) %>%
  bold(i = 1, part = "header") %>% 
  hline_top(part = "header", 
            border = fp_border(color = "red", 
                               width = 3, 
                               style = "solid")) %>% 
  hline(i = 1, 
        part = "header", 
        border = fp_border(color = "black", 
                           width = 0.25, 
                           style = "solid")) %>% 
  hline_top(part = "body", 
            border = fp_border(color = "black", 
                               width = 0.25, 
                               style = "solid")) %>% 
  hline_bottom(part = "body", 
               border = fp_border(color = "black", 
                                  width = 0.25, 
                                  style = "solid")) %>% 
  hline_bottom(part = "footer", 
               border = fp_border(color = "black", 
                                  width = 0.25, 
                                  style = "solid")) %>% 
  border_inner_h(part = "body", 
                 border = fp_border(color = "black", 
                                    width = 0.25, 
                                    style = "dotted")) %>% 
  autofit(part = "body") %>% 
  bg(part = "body", bg = "#f5f5f5") %>% 
  align(part = "all", align = "center") %>% 
  align(j = 1, part = "all", align = "left")



# Save as word .docx
save_as_docx(predict.outpatientdx, path = "Table 2 probability outpatient diagnosis.docx",
             pr_section =
               prop_section(page_size = page_size(orient = "portrait"), 
                            type = "continuous"))




