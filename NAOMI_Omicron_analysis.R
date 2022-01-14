#--------------------------------------------------------------------#
# Crick haemodialysis COVID-19 vaccination study
# Neutralisation tites against SARS-CoV-2 Variants of Concern:
# Omicron and Delta
#
# Edward Carr, Cell Biology of Infection Lab, Francis Crick Institute
# 
#
# Code adapted, from David LV Bauer's analysis of the LEGACY cohort:
# https://github.com/davidlvb/Crick-UCLH-Legacy-VOCs-2021-05/
# https://github.com/davidlvb/Crick-UCLH-Legacy-AZ-VOCs-2021-06/
#
# And our previous haemodialysis report:
# https://github.com/davidlvb/Crick-UCLH-Legacy-AZ-VOCs-2021-06
#
# This work is licensed under a CC BY 4.0 license and is free to use with attribution
#
# Last updated: 12 Jan 2021
# R 4.0.4
#--------------------------------------------------------------------#
## libaries ####
#--------------------------------------------------------------------#
library(tidyverse)
library(furniture)
library(ggpubr)
library(cowplot)

#--------------------------------------------------------------------#
## Read in data ####
#--------------------------------------------------------------------#
load(file = "NAOMI_Omicron_2022-01-12_PUBLIC.RData")

nAb <- nAb %>% 
  ungroup() %>% # ensure ungrouped

#--------------------------------------------------------------------#
## Narrative numbers for manuscript - #s for groups, days intervals ####
#--------------------------------------------------------------------#
## N of patients ##
## (each row is a serum sample in full object)
## (to count patients select a single serum sample for each patient)
nAb %>% group_by(research_identifier) %>% 
  slice_head(n=1) %>% 
  nrow()

## intervals in days after dose 2 ##
nAb %>% filter(str_detect(time_point,"mth6")) %>% 
  summarise(n = n(),
            min = fivenum(daysSinceDose2, na.rm = T)[1],
            Q1 = fivenum(daysSinceDose2, na.rm = T)[2],
            median = fivenum(daysSinceDose2, na.rm = T)[3],
            Q3 = fivenum(daysSinceDose2, na.rm = T)[4],
            max = fivenum(daysSinceDose2, na.rm = T)[5])

## intervals in days after dose 3 ##
nAb %>% filter(str_detect(time_point,"vax3")) %>% 
  summarise(n = n(),
            min = fivenum(daysSinceDose3, na.rm = T)[1],
            Q1 = fivenum(daysSinceDose3, na.rm = T)[2],
            median = fivenum(daysSinceDose3, na.rm = T)[3],
            Q3 = fivenum(daysSinceDose3, na.rm = T)[4],
            max = fivenum(daysSinceDose3, na.rm = T)[5])

## % each vaccine ##
nAb %>% filter(! duplicated(research_identifier)) %>% 
  group_by(vaccine_1) %>% tally()

## % each vaccine //time_point ##
nAb %>% 
  group_by(time_point) %>% tally()

nAb %>% 
  group_by(vaccine_1, time_point) %>% tally()

nAb %>% 
  group_by(vaccine_1, vaccine_3, time_point) %>% tally()

#--------------------------------------------------------------------#
## IC50s - 5 number summaries for 158d after dose 2 ####
#--------------------------------------------------------------------#

## NB:
#### 5 = no inhibition
#### 10 = weak inhibition
#### 5120 = complete inhibition
#### Thus:
#### Q1, median or Q3 of 5 or 10 are described as < 40 in the paper.

## day 158 after dose 2, Delta IC50s:
nAb %>%
  filter(str_detect(time_point, "mth6_bleed")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Delta)) %>%
  group_by(vaccine_1) %>%
  summarise(n = n(),
            min = fivenum(ic50_of_interest, na.rm = T)[1],
            Q1 = fivenum(ic50_of_interest, na.rm = T)[2],
            median = fivenum(ic50_of_interest, na.rm = T)[3],
            Q3 = fivenum(ic50_of_interest, na.rm = T)[4],
            max = fivenum(ic50_of_interest, na.rm = T)[5],
            .groups = "keep")

## day 158 after dose 2, Omicron IC50s:  
nAb %>%
  filter(str_detect(time_point, "mth6_bleed")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  group_by(vaccine_1) %>%
  summarise(n = n(),
            min = fivenum(ic50_of_interest, na.rm = T)[1],
            Q1 = fivenum(ic50_of_interest, na.rm = T)[2],
            median = fivenum(ic50_of_interest, na.rm = T)[3],
            Q3 = fivenum(ic50_of_interest, na.rm = T)[4],
            max = fivenum(ic50_of_interest, na.rm = T)[5],
            .groups = "keep")

#--------------------------------------------------------------------#
## IC50s - 5 number summaries for 27d after dose 3 ####
#--------------------------------------------------------------------#

## Delta IC50s:
nAb %>%
  filter(str_detect(time_point, "vax3")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Delta)) %>%
  group_by(vaccine_1) %>%
  summarise(n = n(),
            min = fivenum(ic50_of_interest, na.rm = T)[1],
            Q1 = fivenum(ic50_of_interest, na.rm = T)[2],
            median = fivenum(ic50_of_interest, na.rm = T)[3],
            Q3 = fivenum(ic50_of_interest, na.rm = T)[4],
            max = fivenum(ic50_of_interest, na.rm = T)[5],
            .groups = "keep")

## Omicron IC50s:  
nAb %>%
  filter(str_detect(time_point, "vax3")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  group_by(vaccine_1) %>%
  summarise(n = n(),
            min = fivenum(ic50_of_interest, na.rm = T)[1],
            Q1 = fivenum(ic50_of_interest, na.rm = T)[2],
            median = fivenum(ic50_of_interest, na.rm = T)[3],
            Q3 = fivenum(ic50_of_interest, na.rm = T)[4],
            max = fivenum(ic50_of_interest, na.rm = T)[5],
            .groups = "keep")

#--------------------------------------------------------------------#
## IC50s - % / N non-responders (<40) after dose 3 ####
#--------------------------------------------------------------------#

## for code, <40 called NR
## Strictly, NR lacks a formal consensus definition.
## Using NR here makes code easier to read cf TRUE/FALSE.

## For Delta
nAb %>%
  filter(str_detect(time_point, "vax3")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Delta)) %>%
  mutate(responder = ifelse(ic50_of_interest<40, 
                            "NR", 
                            "R")) %>%
  group_by(vaccine_1,ic50_of_interest<40,responder) %>%
  summarise(n = n(),
            .groups = "keep")


## Omicron
nAb %>%
  filter(str_detect(time_point, "vax3")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  mutate(responder = ifelse(ic50_of_interest<40, 
                            "NR", 
                            "R")) %>%
  group_by(vaccine_1,ic50_of_interest<40,responder) %>%
  summarise(n = n(),
            .groups = "keep")

#--------------------------------------------------------------------#
## A function to re-shape to 2x2 matrix for McNemar ####
#--------------------------------------------------------------------#
mcnemar.matrix <- function(x, ic50_selected) {
  
  x <- x %>%
    ## make a 2x2 contingency table where the sum of the table is N patients.
    ## [and not 2xN].
    ##
    ##
    ## (A) ensure ic50 numerical ##
    mutate(ic50_of_interest = as.numeric(.data[[ic50_selected]])) %>%
    filter(!is.na(ic50_of_interest)) %>%
    ## (B) define responder / detectable Ab ##
    mutate(responder = ifelse(ic50_of_interest<40, 
                              "NR", 
                              "R")) %>%
    ## (C) select pairs of pre+post ##
    group_by(research_identifier) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    filter(n == 2) %>%
    ## (D) select relevant columns ##
    select(research_identifier, responder, time_point) %>%
    ## (E) Make wider ##
    pivot_wider(id_cols = research_identifier,values_from = responder, names_from = time_point) %>%
    ## (F) Tally the before/after responses ##
    ## drop=FALSE, so empty groups retained ##
    mutate(performance = paste(mth6_bleed, vax3_28d_bleed)) %>%
    mutate(performance = 
             factor(
               performance, 
               # attention to order of levels
               # influences position in 2x2
               levels = c("R R", "NR R", "R NR", "NR NR"))) %>%
    group_by(performance, .drop=FALSE) %>%
    tally() %>%
    ## (G) This is relevant count data as a column 'n' ##
    ## Need to shape this into 2x2 ##
    mutate(preDose = str_remove(as.character(performance), " .*")) %>%
    mutate(postDose = str_remove(as.character(performance), ".* ")) %>%
    select(n, preDose, postDose) %>%
    pivot_wider(id_cols = postDose, names_from = preDose, names_prefix = "preDose_", values_from = n) %>%
    # this is now a 3x2
    as.data.frame(.)
  
  ## mcnemar test is base function - requires matrix ##
  ## data select and change class to matrix ##
  y <- x[,2:3]
  rownames(y) <- paste("postDose", x[,1], sep="_")
  y <- as.matrix(y)
  
  return(y)
}
#--------------------------------------------------------------------#
## McNemar test 1. Delta. both vaccines pre+post dose 3####
#--------------------------------------------------------------------#

# Generate the 2x2
mcn1 <- nAb %>%
  mcnemar.matrix(ic50_selected = "ic50_Delta")

# View 2x2
## R = responder / detectable
## NR = nonresponder / non-detectable
mcn1

# Check the sum of 2x2 is equal to number of paired patients (NOT samples)
sum(mcn1) == nAb %>%
  ## (A) ensure ic50 numerical ##
  mutate(ic50_of_interest = as.numeric(ic50_Delta)) %>%
  filter(!is.na(ic50_of_interest)) %>%
  ## (B) select pairs of pre+post ##
  group_by(research_identifier) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  group_by(research_identifier) %>%
  slice_head(n=1) %>%
  nrow()

# Do McNemar test:
mcnemar.test(mcn1)

#--------------------------------------------------------------------#
## McNemar test 2. Omicron. both vaccines pre+post dose 3####
#--------------------------------------------------------------------#

# Generate the 2x2
mcn2 <- nAb %>%
  mcnemar.matrix(ic50_selected = "ic50_Omicron")

# View 2x2
## R = responder / detectable
## NR = nonresponder / non-detectable
mcn2

# Check the sum of 2x2 is equal to number of paired patients NOT samples.
sum(mcn2) == nAb %>%
  ## (A) ensure ic50 numerical ##
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  filter(!is.na(ic50_of_interest)) %>%
  ## (B) select pairs of pre+post ##
  group_by(research_identifier) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  group_by(research_identifier) %>%
  slice_head(n=1) %>%
  nrow()

# Do McNemar test:
mcnemar.test(mcn2)


#--------------------------------------------------------------------#
## McNemar test 3. Delta. AZD1222 group pre+post dose 3####
#--------------------------------------------------------------------#

# Generate the 2x2
mcn3 <- nAb %>%
  filter(vaccine_1 == "AZD1222") %>%
  mcnemar.matrix(ic50_selected = "ic50_Delta")

# View 2x2
## R = responder / detectable
## NR = nonresponder / non-detectable
mcn3

# Check the sum of 2x2 is equal to number of paired patients NOT samples.
sum(mcn3) == nAb %>%
  filter(vaccine_1 == "AZD1222") %>%
  ## (A) ensure ic50 numerical ##
  mutate(ic50_of_interest = as.numeric(ic50_Delta)) %>%
  filter(!is.na(ic50_of_interest)) %>%
  ## (B) select pairs of pre+post ##
  group_by(research_identifier) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  group_by(research_identifier) %>%
  slice_head(n=1) %>%
  nrow()

# Do McNemar test:
mcnemar.test(mcn3)


#--------------------------------------------------------------------#
## McNemar test 4. Omicron. AZD1222 group pre+post dose 3####
#--------------------------------------------------------------------#

# Generate the 2x2
mcn4 <- nAb %>%
  filter(vaccine_1 == "AZD1222") %>%
  mcnemar.matrix(ic50_selected = "ic50_Omicron")


# View 2x2
## R = responder / detectable
## NR = nonresponder / non-detectable
mcn4

# Check the sum of 2x2 is equal to number of paired patients NOT samples.
sum(mcn4) == nAb %>%
  filter(vaccine_1 == "AZD1222") %>%
  ## (A) ensure ic50 numerical ##
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  filter(!is.na(ic50_of_interest)) %>%
  ## (B) select pairs of pre+post ##
  group_by(research_identifier) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  group_by(research_identifier) %>%
  slice_head(n=1) %>%
  nrow()

# Do McNemar test:
mcnemar.test(mcn4)


#--------------------------------------------------------------------#
## McNemar test 5. Delta. BNT162b2 group pre+post dose 3####
#--------------------------------------------------------------------#

# Generate the 2x2
mcn5 <- nAb %>%
  filter(vaccine_1 == "BNT162b2") %>%
  mcnemar.matrix(ic50_selected = "ic50_Delta")


# View 2x2
## R = responder / detectable
## NR = nonresponder / non-detectable
mcn5

# Check the sum of 2x2 is equal to number of paired patients NOT samples.
sum(mcn5) == nAb %>%
  filter(vaccine_1 == "BNT162b2") %>%
  ## (A) ensure ic50 numerical ##
  mutate(ic50_of_interest = as.numeric(ic50_Delta)) %>%
  filter(!is.na(ic50_of_interest)) %>%
  ## (B) select pairs of pre+post ##
  group_by(research_identifier) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  group_by(research_identifier) %>%
  slice_head(n=1) %>%
  nrow()

# Do McNemar test:
mcnemar.test(mcn5)


#--------------------------------------------------------------------#
## McNemar test 6. Omicron. BNT162b2 group pre+post dose 3####
#--------------------------------------------------------------------#

# Generate the 2x2
mcn6 <- nAb %>%
  filter(vaccine_1 == "BNT162b2") %>%
  mcnemar.matrix(ic50_selected = "ic50_Omicron")


# View 2x2
## R = responder / detectable
## NR = nonresponder / non-detectable
mcn6

# Check the sum of 2x2 is equal to number of paired patients NOT samples.
sum(mcn6) == nAb %>%
  filter(vaccine_1 == "BNT162b2") %>%
  ## (A) ensure ic50 numerical ##
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  filter(!is.na(ic50_of_interest)) %>%
  ## (B) select pairs of pre+post ##
  group_by(research_identifier) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  group_by(research_identifier) %>%
  slice_head(n=1) %>%
  nrow()

# Do McNemar test:
mcnemar.test(mcn6)






#--------------------------------------------------------------------#
## IC50s - 5 number summaries for 27d after dose 3 in immunosuppressed ####
#--------------------------------------------------------------------#

## Omicron IC50s:  
nAb %>%
  filter(immunosuppressed == "Y") %>%
  filter(str_detect(time_point, "vax3")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  group_by(vaccine_1) %>%
  summarise(n = n(),
            min = fivenum(ic50_of_interest, na.rm = T)[1],
            Q1 = fivenum(ic50_of_interest, na.rm = T)[2],
            median = fivenum(ic50_of_interest, na.rm = T)[3],
            Q3 = fivenum(ic50_of_interest, na.rm = T)[4],
            max = fivenum(ic50_of_interest, na.rm = T)[5],
            .groups = "keep")

#--------------------------------------------------------------------#
## IC50s - 5 number summaries for 27d after dose 3 in UN-immunosuppressed ####
#--------------------------------------------------------------------#

## Omicron IC50s:  
nAb %>%
  filter(immunosuppressed == "N") %>%
  filter(str_detect(time_point, "vax3")) %>%
  mutate(ic50_of_interest = as.numeric(ic50_Omicron)) %>%
  group_by(vaccine_1) %>%
  summarise(n = n(),
            min = fivenum(ic50_of_interest, na.rm = T)[1],
            Q1 = fivenum(ic50_of_interest, na.rm = T)[2],
            median = fivenum(ic50_of_interest, na.rm = T)[3],
            Q3 = fivenum(ic50_of_interest, na.rm = T)[4],
            max = fivenum(ic50_of_interest, na.rm = T)[5],
            .groups = "keep")


#--------------------------------------------------------------------#
## Table 1 - demographics ####
#--------------------------------------------------------------------#

## Age and gender are not included in the public dataset.

table1(nAb %>%
         group_by(research_identifier) %>%
         slice_head(n=1) %>%
         ungroup() %>%
         select(
           # age_in_years, 
          # gender,
           vaccine_1, 
                dialysis_centre,
                immunosuppressed),
       splitby = "vaccine_1", 
       na.rm = T,
       test = FALSE,
       # assume parametric distribution:
       param = FALSE, output = "text")


#--------------------------------------------------------------------#
## Figure 1 TOP PANEL ####
#--------------------------------------------------------------------#

## top pane, violins
violins <- nAb %>% 
  ## Long format ###
  pivot_longer(starts_with("ic50"),
               names_to = "variant",
               values_to = "ic50") %>%
  ## clean labels for 'variant' ###
  mutate(variant = str_remove(variant, "ic50_")) %>%
  mutate(variant =
           factor(variant,
                  levels =
                    rev(c("Delta", "Omicron")))) %>%
  mutate(ic50 = as.numeric(ic50)) %>%
  ## plot ##
  ggplot(
    aes(x=time_point, 
        y=ic50, 
        col=vaccine_1, 
        group=research_identifier)) + 
  geom_line(col="black", 
            alpha=0.5,
            position=position_dodge(width=0.2)) + 
  ## violin aes needs group=NULL: ##
  geom_violin(trim=T, 
              aes(group=NULL), 
              alpha=0.5) + 
  geom_point(alpha=0.3,
             size=2,
             position=position_dodge(width=0.2)) + 
  ## colours + facets ##
  scale_colour_manual(values=c("#8E0152", "#276419")) +
  facet_wrap(vaccine_1~variant) +
  ## scales ##
  scale_y_continuous(trans='log2',
                     limits=c(4, 5121),
                     breaks=c(5,10, 40, 256, 1024, 2560, 5120),
                     labels=c("no", "weak","40", "256", "1024","2560", "compl\n"),
                     minor_breaks = c(64,128,512,2048),
                     name=bquote('Virus Neutralisation,'~IC[50]~'')) + 
  scale_x_discrete(labels=c("D2+158d","D3+27d"),
    name="") +
  stat_summary(fun=median,
               inherit.aes=F,
               aes(y=ic50,
                   x=time_point,
                   group=vaccine_1),
               position=position_dodge(width=0.9),
               geom="point",
               color="black", 
               shape=5,
               size=3,
               stroke=2) +
  ## Labels, theme adjustment ##
  ggpubr::theme_pubr(base_size=12, legend="none") +
  ggpubr::rotate_x_text(angle=45) + 
  theme(strip.background=element_blank(),
        strip.text = element_text(size=12),
        panel.grid.major.y = 
          element_line(colour="lightgrey",
                       linetype = 2,
                       size = 0.5),
        panel.grid.minor.y = 
          element_line(colour="lightgrey",
                       linetype = 3,
                       size = 0.5)) +
  labs(title = 
         "Haemodialysis: Omicron nAbTs after third doses",
       subtitle = 
         "AZD1222 or BNT162b2 as doses 1 & 2, all BNT162b2 dose 3\n")


# violins

#--------------------------------------------------------------------#
## Figure 1 BOTTOM PANEL ####
#--------------------------------------------------------------------#

bars <- nAb %>% 
  ## change labels ##
  mutate(time_point = as.character(time_point)) %>%
  mutate(time_point = 
           str_replace(time_point, "mth6_bleed", "D2+158d")) %>%
  mutate(time_point = 
           str_replace(time_point, "vax3_28d_bleed", "D3+27d")) %>%
  mutate(time_point = factor(time_point)) %>%
  ## make long format ##
  pivot_longer(starts_with("ic50"),
               names_to = "variant",
               values_to = "ic50") %>%
  mutate(variant = str_remove(variant, "ic50_")) %>%
  mutate(variant =
           factor(variant,
                  levels =
                    rev(c("Delta", "Omicron")))) %>%
  mutate(ic50 = as.numeric(ic50)) %>%
  filter(!is.na(ic50)) %>%
  mutate(ic50Range = 
           cut(ic50,
               breaks=c(0,40, 256,9999),
               right=FALSE, label=c("<40", "40-256", ">256"))) %>%
  ggplot(aes(x=ic50Range, fill=vaccine_1)) + 
  geom_bar(aes(y=100*..prop.., 
               group=vaccine_1),
           alpha=0.3,
           stat="count",
           position=position_dodge(preserve="single"),
           col="black") +
  coord_flip() +
  ylim(0,100) + 
  ylab("%") + 
  xlab(bquote(
    'Virus Neutralisation,'~IC[50]~'')) +
  scale_fill_manual(values=c("#8E0152", "#276419")) +
  facet_grid(vaccine_1~variant+time_point) +
  theme_pubr(base_size=12, legend = "none") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.grid.major.x = 
          element_line(colour="lightgrey", linetype = 2,
                       size = 0.5),
        panel.grid.minor.x = 
          element_line(colour="lightgrey", linetype = 3,
                       size = 0.5)) + rotate_x_text(angle = 45)




#--------------------------------------------------------------------#
## Figure 1 rendered together, panel labels, exported ####
#--------------------------------------------------------------------#

cowplot::plot_grid(violins, bars,
                   rel_heights = c(2,1),
                   ncol = 1, labels = "AUTO", label_size = 18)

ggsave(filename = "Omicron_HD_Figure1.png", height=9, width=8,
       dpi=300)

