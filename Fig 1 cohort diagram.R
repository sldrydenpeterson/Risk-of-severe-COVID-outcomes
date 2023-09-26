library(tidyverse)

load("enclave.analysis.Rdata")

# CONSORT type flow diagram

data <- tibble(x= 1:200, y= 1:200)
head(data)

c <-  data %>% 
  ggplot(aes(x, y)) +
  scale_x_continuous(minor_breaks = seq(10, 200, 10)) +
  scale_y_continuous(minor_breaks = seq(10, 240, 10)) +
  theme_linedraw()
c

# starting box
start.total<- enclave.analysis %>%
  count()
c<- c +
  geom_rect(xmin = 50, xmax=150, ymin=220, ymax=230, color='black',
            fill='#DEDFE4', size=0.25) +
  annotate('text', x= 100, y=225,  size=2.5, hjust=0.5, vjust=0.5,
           label= str_wrap(paste0(prettyNum(start.total, big.mark =  ","), " Patients (age 18 and older) recorded positive SARS-CoV-2 testing between June 1 and December 31, 2022"),
                           indent =0 , width =60, exdent = 0)) +
  geom_segment(
    x=100, xend=100, y=220, yend=200.5, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=100, xend=124.5, y=210, yend=210, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) 

c

# box of residence in MA/NH
residents<- enclave.analysis %>%    
  filter(state == "Massachusetts" | state == "New Hampshire") %>% 
   count()

c<- c +
  geom_rect(xmin = 50, xmax=150, ymin=190, ymax=200, color='black',
            fill='#DEDFE4', size=0.25) +
  annotate('text', x= 100, y=195,  size=2.5, hjust=0.5, vjust=0.5,
           label= str_wrap(paste0(prettyNum(residents, big.mark =  ","), " Patients resident of Massachusetts or New Hampshire"),
                           indent =0 , width =60, exdent = 0)) + 
  
  geom_rect(xmin = 125, xmax=182, ymin=205, ymax=215, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 127, y=210,  size=2.5, hjust=0, vjust=0.5,
           label= str_wrap(paste0(prettyNum(start.total-residents, big.mark =  ","), " Patients residing outside of study area"),
                           indent =0 , width = 38 , exdent = 0)) +
  geom_segment(
    x=100, xend=100, y=190, yend=170.5, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=100, xend=124.5, y=180, yend=180, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) 

c


# box of high risk
increased.risk <- enclave.analysis %>%    
  filter(state == "Massachusetts" | state == "New Hampshire") %>%
  filter(CDC.increasedrisk == 1 ) %>%
  count()

c<- c +
  geom_rect(xmin = 50, xmax=150, ymin=160, ymax=170, color='black',
            fill='#DEDFE4', size=0.25) +
  annotate('text', x= 100, y=165,  size=2.5, hjust=0.5, vjust=0.5,
           label= str_wrap(paste0(prettyNum(increased.risk, big.mark =  ","), " Patients with one or more CDC increased risk criteria"),
                           indent =0 , width =60, exdent = 0))  + 
  
  geom_rect(xmin = 125, xmax=182, ymin=175, ymax=185, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 127, y=180,  size=2.5, hjust=0, vjust=0.5,
           label= str_wrap(paste0(prettyNum(residents-increased.risk, big.mark =  ","), " Patients not known to meet CDC increased risk criteria"),
                           indent =0 , width =38, exdent = 0)) +
#   geom_segment(
#     x=100, xend=124.5, y=150, yend=150, 
#     size=0.15, linejoin = "mitre", lineend = "butt",
#     arrow = arrow(length = unit(1, "mm"), type= "closed")) 
# c


geom_segment(
  x=100, xend=100, y=160, yend=140, 
  size=0.15, linejoin = "mitre", lineend = "butt") +
  
  geom_segment(
    x=60, xend=140, y=140, yend=140, 
    size=0.15, linejoin = "mitre", lineend = "butt") +
  
  geom_segment(
    x=60, xend=60, y=140, yend=130.5, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  
  geom_segment(
    x=140, xend=140, y=140, yend=130.5, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed"))
c




# boxes of outcomes
outcome <- enclave.analysis %>%    
  filter(state == "Massachusetts" | state == "New Hampshire") %>%
  filter(CDC.increasedrisk == 1 ) %>%
  group_by(covid_admit_death) %>%
  count()

c<- c +
  geom_rect(xmin = 35, xmax=85, ymin=110, ymax=130, color='black',
            fill='#DEDFE4', size=0.25) +
  annotate('text', x= 60, y=120,  size=2.5, hjust=0.5, vjust=0.5,
           label= str_wrap(paste0(prettyNum( outcome[1,2], big.mark =  ","), " Patients without a severe outcome of COVID-19"),
                           indent =0 , width = 30, exdent = 0)) + 
  
  geom_rect(xmin = 115, xmax=165, ymin= 110, ymax=130, color='black',
            fill='#DEDFE4', size=0.25)+
  annotate('text', x= 140, y=120,  size=2.5, hjust=0.5, vjust=0.5,
           label= str_wrap(paste0(prettyNum( outcome[2,2], big.mark =  ","), " Patients with a severe outcome of COVID-19"),
                           indent =0 , width =30, exdent = 0)) + 
  
  theme_void()


c + coord_cartesian(ylim = c(100,240))





ggsave("Cohort diagram Fig 1.pdf", width = 8, height = 8, dpi = 600, bg = "white")

rm(list = ls())
