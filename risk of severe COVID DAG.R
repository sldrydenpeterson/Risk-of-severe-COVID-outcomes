library(tidyverse)
library(ggdag)


severecovid_dag <- list(
  x = c(hosp.death = 1.25, covid.dx = 1, comorbidities = 0, vaccine = 0.5, ses = 0.25, race.eth = 0.25),
  y = c(hosp.death = 1, covid.dx = 1.05, comorbidities = 1, vaccine = .95, ses = 0.5, race.eth = 1.5)
)

severe.covid.risk <- 
  dagify(covid.dx ~ ses,
         covid.dx ~ race.eth,
         covid.dx ~ ses,
         covid.dx ~ vaccine,
         covid.dx ~ hosp.death,
         covid.dx ~ comorbidities,
         vaccine ~ race.eth,
         vaccine ~ ses,
         vaccine ~ comorbidities ,
         comorbidities ~ race.eth,
         ses ~ comorbidities,
         ses ~ race.eth,
         exposure = "comorbidities",
         outcome = "covid.dx",
         coords = severecovid_dag,
         
         labels = c(
           "covid.dx" = "Documented COVID\ndiagnosis",
           "hosp.death" = "Hospitalization or\ndeath",
           "comorbidities" = "Comorbidities\nand age",
           "ses" = "Socioeconomic\nposition",
           "race.eth" = "Race and\nethnicity",
           "vaccine" = "Vaccination & prior\ninfection status") ) %>%
  tidy_dagitty() %>%
  mutate(confounder = fct_relevel(if_else(name %in% c("hosp.death", "comorbidities"), "causal", "covariates" ),  "causal", "covariates"),
         conditional = if_else(name %in% c("covid.dx"), "condit", "uncondit"))

## graphics


severe.covid.risk %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend, shape = conditional, label = label, color = confounder)) +
  geom_dag_edges(aes(label = NA)) +
  geom_dag_node() +
  geom_dag_label_repel(
    ggplot2::aes( label = label, fill = confounder ), 
    col = "white", show.legend = FALSE
  ) +
  scale_shape_manual(values = c(15,19)) +
  ggsci::scale_color_uchicago() +
  ggsci::scale_fill_uchicago() +
  guides(color = "none", fill= "none", shape = "none") +
  expand_plot(expand_y = expansion(c(0.2, 0.2))) +
  theme_dag_blank()

ggsave("severe covid risk DAG.pdf", width = 8, height = 8)



# 
# 
# 
# ggdag_drelationship(severe.covid.risk, controlling_for = c("vaccine", "ses", "race.eth", "ses"),
#                     text = FALSE, use_labels = "label",
#                     stylized = TRUE, node_size = 10,   text_size = 2.8, edge_type="link_arc", collider_lines = TRUE)
# 
# ggdag_adjustment_set(severe.covid.risk, text = FALSE, use_labels = "label",
#                      shadow = TRUE, stylized = TRUE, node_size = 10,   text_size = 2.8)
# 
# 
# ggdag(severe.covid.risk, text = FALSE, use_labels = "label", edge_type = "link_arc")
# 
# 
# ggdag_paths(severe.covid.risk, text = FALSE, use_labels = "label", shadow = TRUE,
# 
#             #directed=false indicates all paths, regardless of direction
#              directed=TRUE,
# 
#             #changing node/text sizes for aesthetics
#             node_size = 8,   text_size = 2, stylized = TRUE)