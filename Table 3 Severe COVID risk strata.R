library(tableone)
library(flextable)
library(officer)


# Set Table header/footer
header <- "Table 3. Risk of severe COVID-19 by comorbidity, immunity, and outpatient treatment status— Mass General Brigham, June to December 2022"
# footer <- paste0(" Hospitalization within 14 days or death within 28 days of COVID diagnosis. The number of individuals with COVID-19 at risk, including those without a recorded positive test, estimated via modified capture-recapture methodology utilizing strata-specific weights. The methodology used does not provide adjustment for possible differential risk within strata between individuals utilizing and not utilizing Confidence limits estimated via bootstrapping with 2000 replicates.\n",
#           "† Monoclonal antibody screening score (MASS) calculated as age 65 years and older (2 points), BMI 35 and higher (2), diabetes (2), chronic kidney disease (3), cardiovascular disease in a patient 55 years and older (2), chronic respiratory disease in a patient 55 years and older (3), hypertension in a patient 55 years and older (1), and immunocompromised status (3). Risk for individuals with severe immunocompromise— organ or stem cell transplantation, hematologic malignancy, or receiving mTOR inhibitors, cyclosporine, mycophenolate, or anti-CD20 therapy— estimated separately.")

footer <- as_paragraph("Abbreviations: CI, confidence interval; NNTB, number needed to treat to benefit; RD, risk difference\n",
                       as_sup("a"), " Hospitalization within 14 days or death within 28 days of COVID-19 diagnosis.\n",
                       as_sup("b"), " Monoclonal antibody screening score (MASS) calculated as age 65 years and older (2 points), BMI 35 and higher (2), diabetes (2), chronic kidney disease (3), cardiovascular disease in a patient 55 years and older (2), chronic respiratory disease in a patient 55 years and older (3), hypertension in a patient 55 years and older (1), and immunocompromised status (3). Risk for individuals with severe immunocompromise— organ or stem cell transplantation, hematologic malignancy, or receiving mTOR inhibitors, cyclosporine, mycophenolate, or anti-CD20 therapy— estimated separately.\n",
                       as_sup("c"), " Outpatient treatments during the study period included oral nirmatrelvir-ritonavir (90%), oral molnupiravir (2%), intravenous bebtelovimab (5%), and intravenous remdesivir (3%).\n",
                       as_sup("d"), " Strata-specific estimates calculated through marginal inverse-probability weighted models to reduce expected bias in access and acceptance to outpatient treatment by vaccination status, race and ethnicity, neighborhood disadvantage, and age.")
spanner <- as_paragraph(as_chunk(c("", "Hospitalization and/or death risk (95% CI)", "")))


load("severecovid.risk.est.untreat.all.Rdata")
load("severecovid.risk.est.treat.all.Rdata")
load("nnt.all.Rdata")

flextable_3 <- flextable(
                left_join(
                left_join(
                    severecovid.risk.est.untreat.all, severecovid.risk.est.treat.all, by= c("MASS score", "Vaccination status")),
                    nnt.all, by = "MASS score") %>%
                mutate(
                  'Prescribed outpatient COVID treatment' = if_else(`Prescribed outpatient COVID treatment` == "0% (0% to 0%)", "0% (unable to estimate)", 
                                                                    `Prescribed outpatient COVID treatment`),
                  `MASS score` = case_when(
                    `MASS score` == "MASS 3 or less" ~ "MASS 3 or less",
                    `MASS score` == "MASS 4 and 5" ~ "MASS 4 and 5",
                    `MASS score` == "MASS 6 or greater" ~ "MASS 6 or greater",
                    `MASS score` == "Severe immunocompromise" ~ "Severe\nimmunocompromise",
                    TRUE ~ "test"
                    
                  )
                ) ) %>%
  merge_v( j = c("MASS score", "rd.nnt")) %>%            
  add_header_row(colwidths = c(2, 2, 1), values= c("", "Hospitalization and/or death risk (95% CI)", "")) %>%
  # add footnotes
  mk_par(
    part = "header", j = 3, i = 1,
    value = as_paragraph(as_b("Hospitalization and/or death risk (95% CI)"), as_sup("a"))) %>%
  mk_par(
    part = "header", j = 1, i = 2,
    value = as_paragraph(as_b("MASS score"), as_sup("b"))) %>%
  mk_par(
    part = "header", j = 2, i = 2,
    value = as_paragraph(as_b("Vaccination status"))) %>%
  mk_par(
    part = "header", j = 3, i = 2,
    value = as_paragraph(as_b("No outpatient COVID-19 antiviral"), as_sup("c"))) %>%
  mk_par(
    part = "header", j = 4, i = 2,
    value = as_paragraph(as_b("Prescribed outpatient COVID-19 antiviral"))) %>%
  mk_par(
    part = "header", j = 5, i = 2,
    value = as_paragraph(as_b("Adjusted risk difference and number needed to treat to benefit\n(95% CI)"), as_sup("d"))) %>%
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
save_as_docx(flextable_3, path = "Table 3 Severe COVID outcomes risk strata.docx",
             
             pr_section =
               prop_section(page_size = page_size(orient = "landscape"), 
                            type = "continuous"))

sect_properties <- prop_section(
  page_size = page_size(orient = "landscape",
                        width = 9.5, height = 11.8),  #allows to fit on one page with output
  type = "continuous",
  page_margins = page_mar()
)

save_as_docx(flextable_3,
             path =  "Table 3 Severe COVID outcomes risk strata.docx", 
             pr_section = sect_properties)

rm(list = ls())