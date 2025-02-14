rm(list = ls())

library(tidyverse)
library(marginaleffects)
library(boot)
library(flextable)
library(officer)

load("enclave.eligible.wt.Rdata")

# risk among vaccinated with MASS 3 or less
risk.figure.dt <- rbind(
  tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 3 or less")#,
                    # mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     )) ) %>% mutate(Label =  "MASS 3 or less", predictor = "Comorbidities") %>%
  select(Label, estimate, conf.low, conf.high, predictor),
tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 4 and 5")#,
                    # mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     )) 
) %>% mutate(Label =  "MASS 4 and 5", predictor = "Comorbidities") %>%
  select(Label, estimate, conf.low, conf.high, predictor),
tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("MASS 6 or greater")#,
                    # mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     )) 
) %>% mutate(Label =  "MASS 6 or greater", predictor = "Comorbidities") %>%
  select(Label, estimate, conf.low, conf.high, predictor),
tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ MASS.cat + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0) ),
  newdata = datagrid(MASS.cat = c("Severe immunocompromise")#,
                     #mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     )) 
) %>% mutate(Label =  "Severe immunocompromise", predictor = "Comorbidities") %>%
  select(Label, estimate, conf.low, conf.high, predictor) %>%
  add_row(predictor = "Comorbidities", Label = "    "),

#age cat
tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ age.decade + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
          mutate(age.decade = case_when(
            age.at.covid < 30 ~ "18 to 29",
            age.at.covid < 50~ "30 to 49",
            age.at.covid < 65~ "50 to 64",
            age.at.covid < 75~ "65 to 74",
            age.at.covid < 85~ "75 to 84",
            TRUE ~ "85 and older"))),
  newdata = datagrid(age.decade = c("18 to 29")#, mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "18 to 29 years", predictor = "Age Groups") %>%
  select(Label, estimate, conf.low, conf.high, predictor),

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ age.decade + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(age.decade = case_when(
        age.at.covid < 30 ~ "18 to 29",
        age.at.covid < 50~ "30 to 49",
        age.at.covid < 65~ "50 to 64",
        age.at.covid < 75~ "65 to 74",
        age.at.covid < 85~ "75 to 84",
        TRUE ~ "85 and older"))),
  newdata = datagrid(age.decade = c("30 to 49")#, mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "30 to 49 years", predictor = "Age Groups") %>%
  select(Label, estimate, conf.low, conf.high, predictor),

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ age.decade + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(age.decade = case_when(
        age.at.covid < 30 ~ "18 to 29",
        age.at.covid < 50~ "30 to 49",
        age.at.covid < 65~ "50 to 64",
        age.at.covid < 75~ "65 to 74",
        age.at.covid < 85~ "75 to 84",
        TRUE ~ "85 and older"))),
  newdata = datagrid(age.decade = c("50 to 64")#,# mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "50 to 64 years", predictor = "Age Groups") %>%
  select(Label, estimate, conf.low, conf.high, predictor),

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ age.decade + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(age.decade = case_when(
        age.at.covid < 30 ~ "18 to 29",
        age.at.covid < 50~ "30 to 49",
        age.at.covid < 65~ "50 to 64",
        age.at.covid < 75~ "65 to 74",
        age.at.covid < 85~ "75 to 84",
        TRUE ~ "85 and older"))),
  newdata = datagrid(age.decade = c("65 to 74")#,
                    # mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                    ))) %>% 
  mutate(Label =  "65 to 74 years", predictor = "Age Groups") %>%
  select(Label, estimate, conf.low, conf.high, predictor) ,

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ age.decade + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(age.decade = case_when(
        age.at.covid < 30 ~ "18 to 29",
        age.at.covid < 50~ "30 to 49",
        age.at.covid < 65~ "50 to 64",
        age.at.covid < 75~ "65 to 74",
        age.at.covid < 85~ "75 to 84",
        TRUE ~ "85 and older"))),
  newdata = datagrid(age.decade = c("75 to 84")#,
                     #mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "75 to 84 years", predictor = "Age Groups") %>%
  select(Label, estimate, conf.low, conf.high, predictor),

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ age.decade + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(age.decade = case_when(
        age.at.covid < 30 ~ "18 to 29",
        age.at.covid < 50~ "30 to 49",
        age.at.covid < 65~ "50 to 64",
        age.at.covid < 75~ "65 to 74",
        age.at.covid < 85~ "75 to 84",
        TRUE ~ "85 and older"))),
  newdata = datagrid(age.decade = c("85 and older")#,
                    #mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                    ))) %>% 
  mutate(Label =  "85 years and older", predictor = "Age Groups") %>%
  select(Label, estimate, conf.low, conf.high, predictor) %>%
  add_row(predictor = "Age Groups", Label = "  "),

# WHO categories

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ WHOrisk + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(WHOrisk = case_when(
        WHOrisk == "1" ~ "High risk", 
        WHOrisk == "2" ~ "Moderate risk", 
        WHOrisk == "3" ~ "Low risk" 
      ))),
  newdata = datagrid(WHOrisk = c("Low risk")#,
                     #mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "Low risk", predictor = "WHO Risk Category") %>%
  select(Label, estimate, conf.low, conf.high, predictor),

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ WHOrisk + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(WHOrisk = case_when(
        WHOrisk == "1" ~ "High risk", 
        WHOrisk == "2" ~ "Moderate risk", 
        WHOrisk == "3" ~ "Low risk" 
      ))),
  newdata = datagrid(WHOrisk = c("Moderate risk")#,
                     #mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "Moderate risk", predictor = "WHO Risk Category") %>%
  select(Label, estimate, conf.low, conf.high, predictor),

tibble(avg_predictions(
  geepack::geeglm(
    covid_admit_death ~ WHOrisk + mod.vax.status,
    family = poisson(link=log),
    weights = outpt.wts,
    id = seqn,
    corstr = "independence",
    data = enclave.eligible.wt %>% filter(covid_treat == 0)  %>%
      mutate(
        WHOrisk = case_when(
        WHOrisk == "1" ~ "High risk", 
        WHOrisk == "2" ~ "Moderate risk", 
        WHOrisk == "3" ~ "Low risk" ))),
  newdata = datagrid(WHOrisk = c("High risk")#,
                     #mod.vax.status = c("Vaccinated, last dose or COVID-19 ≥ 8 months prior", "Vaccinated, last dose or COVID-19 < 8 months prior")
                     ))) %>% 
  mutate(Label =  "High risk", predictor = "WHO Risk Category") %>%
  select(Label, estimate, conf.low, conf.high, predictor)) %>%
  mutate(Label = fct_relevel(Label, "High risk", "Moderate risk", "Low risk", "  ", 
                             "Severe immunocompromise", 'MASS 6 or greater', "MASS 4 and 5", "MASS 3 or less", "    ", 
                             "85 years and older", "75 to 84 years", "65 to 74 years", "50 to 64 years", "30 to 49 years", "18 to 29 years"),
         estimate = estimate*100,
         conf.high = conf.high*100,
         conf.low = conf.low*100)

risk.figure.dt

fig2 <- ggplot(risk.figure.dt ) +
  # ggalt::geom_dumbbell(aes( x = conf.low, xend = conf.high, y = Label), size_x = 5, size_xend = 5, size =5,
  #                      color="#e3e2e1", 
  #               colour_x = "#bad744", colour_xend = "#bad744",
  #               dot_guide=TRUE, dot_guide_size=.25) + 
  geom_segment(aes(x = 0, xend = conf.low, y =Label), size = 0.25, linetype = 3) +
  geom_segment(aes(x= conf.low, xend = conf.high, y = Label), size = 5, color = "#e3e2e1") +
  scale_color_manual(values = c("#4393c3", "#bad744", "#f4a582"))+ 
  geom_point(aes(x = conf.low, y = Label, color = predictor), size = 5) +
  geom_point(aes(x = conf.high, y = Label, color = predictor), size = 5) +
  geom_point(aes(x = estimate, y = Label), shape = 124, size = 6) +
  geom_text(aes( x = estimate, y = Label, label = paste0(gtsummary::style_sigfig(estimate, digit = 2), "%")), nudge_x = 0.08, nudge_y = 0.3) +
  theme_classic() + 
  annotate("rect",  xmin = 5.75, xmax = 6.75, ymin = 9.25, ymax = 15.5, fill = "#4393c3") +
  annotate("text",  x = 6.25, y = 12.4, label = "Age\ngroups", angle = 270, color = "white", size =6) +
  annotate("rect",   xmin = 5.75, xmax = 6.75, ymin = 4.25, ymax = 8.75, fill = "#bad744") +
  annotate("text",  x = 6.25, y = 6.5, label = "Number/severity\nof comorbidities", angle = 270, color = "white", size =6) +
  annotate("rect",   xmin = 5.75, xmax = 6.75, ymin = 0.25, ymax = 3.75, fill = "#f4a582") +
  annotate("text",  x = 6.25, y = 2, label = "WHO risk\ncategories", angle = 270, color = "white", size =6) +
  coord_cartesian(xlim = c(0,7), ylim = c(0,15), clip = "off") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5), labels =c("0", "1%", "2%", "3%", "4%", "5%") ,minor_breaks = NULL, expand = c(0, 0)) +
  guides(color = "none") +
  labs(title="Hospitalization or death following COVID-19 infection",
       subtitle = "Increased-risk individuals without outpatient antiviral therapy— 236,757 total infections (recorded and unrecorded)") +
  bbc_style() +
  theme(plot.title = element_text(size =25), 
        plot.subtitle = element_text(size =16),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color="#cbcbcb", linewidth = 0.25))

#ggsave("Figure 2.pdf", width = 12, height = 12)
 
footer <- ""
finalise_plot(
  plot_name = fig2,
  source = "",
  save_filepath = "Figure 2.png",
  width_pixels = 900,
  height_pixels = 900
)

