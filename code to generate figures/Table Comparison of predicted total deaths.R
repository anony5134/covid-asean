## Input data are provided in the Epidemia_Models/epidemia_rt_output folder, so it's not necessary to run the Epidemia models to run this code
###################################################
## calculate numbers and make table1 delta deaths with CF for 2020 and 2021
###################################################
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(jcolors)
library(kableExtra)
library(gt)
###################################################
###################################################
##### COVID-19
###################################################
df_noSI_2020 <-readRDS("/Users///COVID-Counterfactual/Top5/noSI/counterfactual_2020_top5_iter200_noSI.rds")
df_SI_2020 <-readRDS("/Users///COVID-Counterfactual/Top5/SI/counterfactual_2020_top5_iter200_SI.rds")
df_noSI_2021 <- readRDS("/Users///COVID-Counterfactual/Top5/2021_noSI/counterfactual_2021_top5_iter200_noSI.rds")
df_SI_2021 <-readRDS("/Users///COVID-Counterfactual/Top5/2021_SI/counterfactual_2021_top5_iter200_SI.rds")
#
# if df created prior to edits at line 626
#
df<- df_noSI_2020
source("/Users///COVID-Counterfactual/Top5/table1_transform_df_plotting_top5.R")
df_noSI_2020 <- df
#
df<- df_SI_2020
source("/Users///COVID-Counterfactual/Top5/table1_transform_df_plotting_top5.R")
df_SI_2020 <- df
#
df<- df_noSI_2021
source("/Users///COVID-Counterfactual/Top5/table1_transform_df_plotting_top5.R")
df_noSI_2021 <- df
#
df<- df_SI_2021
source("/Users///COVID-Counterfactual/Top5/table1_transform_df_plotting_top5.R")
df_SI_2021 <- df
remove(df)
#
#TOP 5 ("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
#
# Group by "group" and calculate the sum of "deaths" for each group
grouped_data_SI_2020 <- df_SI_2020 %>%
  group_by(group) %>%
  summarize(Total_Deaths_SI_2020 = round(sum(deaths),), Modelled_deaths_SI_2020 = round(sum(deaths_median),))
grouped_data_SI_2020 <- grouped_data_SI_2020 %>%
  mutate(
    Averted_SI_2020 = (Total_Deaths_SI_2020 - Modelled_deaths_SI_2020),
  ) %>%
rename(CounterFactual = group)
#
grouped_data_noSI_2020 <- df_noSI_2020 %>%
  group_by(group) %>%
  summarize(Total_Deaths_noSI_2020 = round(sum(deaths),), Modelled_deaths_noSI_2020 = round(sum(deaths_median),))
grouped_data_noSI_2020 <- grouped_data_noSI_2020 %>%
  mutate(
    Averted_noSI_2020 = (Total_Deaths_noSI_2020 - Modelled_deaths_noSI_2020),
  ) %>%
  rename(CounterFactual = group)
#2021
grouped_data_SI_2021 <- df_SI_2021 %>%
  group_by(group) %>%
  summarize(Total_Deaths_SI_2021 = round(sum(deaths),), Modelled_deaths_SI_2021 = round(sum(deaths_median),))
grouped_data_SI_2021 <- grouped_data_SI_2021 %>%
  mutate(
    Averted_SI_2021 = (Total_Deaths_SI_2021 - Modelled_deaths_SI_2021),
  ) %>%
  rename(CounterFactual = group)
#
grouped_data_noSI_2021 <- df_noSI_2021 %>%
  group_by(group) %>%
  summarize(Total_Deaths_noSI_2021 = round(sum(deaths),), Modelled_deaths_noSI_2021 = round(sum(deaths_median),))
grouped_data_noSI_2021 <- grouped_data_noSI_2021 %>%
  mutate(
    Averted_noSI_2021 = (Total_Deaths_noSI_2021 - Modelled_deaths_noSI_2021),
  ) %>%
  rename(CounterFactual = group)
# create single table based on counterfactuals
Table1_2020 <- merge(grouped_data_noSI_2020, grouped_data_SI_2020, by = "CounterFactual")
Table1_2020 <- Table1_2020[order(Table1_2020$CounterFactual), ]
Table1_2020 <- Table1_2020 %>%
  arrange(case_when(
    startsWith(CounterFactual, "F") ~ 1,
    startsWith(CounterFactual, "A") ~ 2,
    startsWith(CounterFactual, "R") ~ 3,
    TRUE ~ 4
  ))
#
gt_table1_2020 <- gt(Table1_2020) %>%
  data_color(
    columns = c("Averted_noSI_2020","Averted_SI_2020"),
    method = "bin",
    bins  = c(-Inf, -0.1),
    palette = c("red"), ) %>%
  cols_align(align = "center") %>%
  tab_header( title = "Predicted deaths from Epidemia Models for ASEAN countries (2020)") %>%
  tab_spanner( label = "Modelled without SI",columns = c(Total_Deaths_noSI_2020,Modelled_deaths_noSI_2020,Averted_noSI_2020	))  %>%
  cols_label("Total_Deaths_noSI_2020" = "Total Deaths") %>%
  cols_label("Modelled_deaths_noSI_2020" = "Modelled deaths") %>%
  cols_label("Averted_noSI_2020" = "Averted deaths") %>%
  tab_spanner( label = "Modelled using SI",columns = c(Total_Deaths_SI_2020,Modelled_deaths_SI_2020,Averted_SI_2020	))  %>%
  cols_label("Total_Deaths_SI_2020" = "Total Deaths") %>%
  cols_label("Modelled_deaths_SI_2020" = "Modelled deaths") %>%
  cols_label("Averted_SI_2020" = "Averted deaths") %>%
  cols_width(everything() ~ px(150))

gtsave(gt_table1_2020,"/Users///COVID-Counterfactual/Top5/table1_2020.html")

#
# create single table based on counterfactuals
Table1_2021 <- merge(grouped_data_noSI_2021, grouped_data_SI_2021, by = "CounterFactual")
Table1_2021 <- Table1_2021[order(Table1_2021$CounterFactual), ]
Table1_2021 <- Table1_2021 %>%
  arrange(case_when(
    startsWith(CounterFactual, "F") ~ 1,
    startsWith(CounterFactual, "A") ~ 2,
    startsWith(CounterFactual, "R") ~ 3,
    TRUE ~ 4
  ))
#
gt_table1_2021 <- gt(Table1_2021) %>%
  data_color(
    columns = c("Averted_noSI_2021","Averted_SI_2021"),
    method = "bin",
    bins  = c(-Inf, -0.1), 
    palette = c("red"), 
    ) %>%
  cols_align(align = "center") %>%
  tab_header( title = "Predicted deaths from Epidemia Models for ASEAN countries(2021)") %>%
  tab_spanner( label = "Modelled without SI",columns = c(Total_Deaths_noSI_2021,Modelled_deaths_noSI_2021,Averted_noSI_2021))  %>%
  cols_label("Total_Deaths_noSI_2021" = "Total Deaths") %>%
  cols_label("Modelled_deaths_noSI_2021" = "Modelled deaths") %>%
  cols_label("Averted_noSI_2021" = "Averted deaths") %>%
  tab_spanner( label = "Modelled using SI",columns = c(Total_Deaths_SI_2021,Modelled_deaths_SI_2021,Averted_SI_2021))  %>%
  cols_label("Total_Deaths_SI_2021" = "Total Deaths") %>%
  cols_label("Modelled_deaths_SI_2021" = "Modelled deaths") %>%
  cols_label("Averted_SI_2021" = "Averted deaths")



gtsave(gt_table1_2021,"/Users///COVID-Counterfactual/Top5/table1_2021.html")

