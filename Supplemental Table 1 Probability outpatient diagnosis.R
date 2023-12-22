library(flextable)
library(officer)


load("inpatients.Rdata")

# Set Table header/footer
header <- as_paragraph("Supplemental Table 1. Odds ratios of preceding outpatient diagnosis of COVID-19 among patients with severe outcomes",
                        as_sup("a"), " by selected characteristics— Mass General Brigham, June to December 2022", "\n")
# footer <- paste0("* Recorded positive COVID-19 test result two or more calendar days prior to COVID hospitalization or death.\n", 
#                  "† Adjusted odds ratios from multivariable logistic model of apriori-selected characteristics are shown.")

footer <- as_paragraph("Abbreviations: ADI, Area Deprivation Index; CI, confidence interval\n",
                       as_sup("a"), " Recorded positive COVID-19 test result two or more calendar days prior to COVID hospitalization or death.\n",
                       as_sup("b"), " Adjusted odds ratios from multivariable logistic model of apriori-selected characteristics are shown.\n")



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
  gtsummary::as_flex_table() %>%
  add_header_lines(header) %>%
  add_footer_lines(footer) %>%
  bg(part = "body", bg = "#ebf4f1", i = c(2,4,6,8, 10,12,14,16, 18, 20, 22))%>%
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
  mk_par(
    part = "header", j = 2, i = 2,
    value = as_paragraph(as_b("Odds Ratio"), as_sup("b"))) 



# Save as .docx
sect_properties <- prop_section(
  page_size = page_size(orient = "portrait",
                        width = 9.5, height = 13),  #allows to fit on one page with output
  type = "continuous",
  page_margins = page_mar()
)

save_as_docx(predict.outpatientdx,
             path =  "Supplmental Table 1 probability outpatient diagnosis.docx", 
             pr_section = sect_properties)

rm(list = ls())
