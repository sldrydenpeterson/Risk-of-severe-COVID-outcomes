library(tableone)
library(flextable)
library(officer)


# Set Table header/footer
header <- "Table 3. Risk of COVID-19 hospitalization or death by comorbidity and immunity status— Mass General Brigham, June to December 2022"
footer <- paste0("* Hospitalization within 14 days or death within 28 days of COVID diagnosis. The number of individuals with COVID-19 at risk, including those without a recorded positive test, estimated via modified capture-recapture methodology utilizing strata-specific weights. The methodology used does not provide adjustment for possible differential risk within strata between individuals utilizing and not utilizing outpatient COVID treatment (i.e., nirmatrelvir-ritonavir, molnupiravir, bebtelovimab, and remdesivir). Confidence limits estimated via bootstrapping with 2000 replicates.\n",
          "† Monoclonal antibody screening score (MASS) calculated as age 65 years and older (2 points), BMI 35 and higher (2), diabetes (2), chronic kidney disease (3), cardiovascular disease in a patient 55 years and older (2), chronic respiratory disease in a patient 55 years and older (3), hypertension in a patient 55 years and older (1), and immunocompromised status (3). Risk for individuals with severe immunocompromise— organ or stem cell transplantation, hematologic malignancy, or receiving mTOR inhibitors, cyclosporine, mycophenolate, or anti-CD20 therapy— estimated separately.")

flextable_1 <- flextable(left_join(severecovid.risk.est.untreat.all, severecovid.risk.est.treat.all, by= c("MASS score", "Vaccination status")) %>%
                           mutate('Prescribed outpatient COVID treatment' = if_else(`Prescribed outpatient COVID treatment` == "0% (0% to 0%)", "0% (unable to estimate)", `Prescribed outpatient COVID treatment`)) %>%
                           rename( "MASS score†" = "MASS score")
                              ) %>%
                  
  add_header_row(colwidths = c(2, 2), values= c("", "Hospitalization and/or death risk (95% CI)*"))  %>%
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
save_as_docx(flextable_1, path = "Table 3 Severe COVID outcomes risk strata.docx",
             pr_section =
               prop_section(page_size = page_size(orient = "landscape"), 
                            type = "continuous"))




