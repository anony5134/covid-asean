library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(parallel)
library(parallelly)
library(zoo)
library(lubridate)
library(tseries)
library(cowplot)
library(jcolors)
library(tidyquant)
library(timetk)
library(epidemia)
source("/Users///COVID-Counterfactual/utils.R")
source("/Users///COVID-Counterfactual/ts_block_bootstrap_sp.R")
source("/Users///COVID-Counterfactual/get_quantiles.R")
## Input data are provided in the Epidemia_Models/epidemia_rt_output folder, so it's not necessary to run the Epidemia models to run this code
###################################################
## Plot Rt and incidence (Figure 1b)
###################################################
###################################################
fit=readRDS("/Users///COVID-Counterfactual/Top5/SI/fit_2020_top5_iter200_SI.rds")
all_covid_rt = extract_Rt(fit)
all_covid_rt$organism <- "SARS-CoV-2"
head(all_covid_rt)
#
all_covid_rt_avg <- all_covid_rt %>%
  group_by(group, tag, level, organism) %>%
  mutate_at(c("lower", "upper", "median"), ~ zoo::rollmean(., k = 15, align = "center", fill = NA)) %>%
  ungroup() %>%
  filter(date < as.Date("2021-01-01")) %>%
  filter(!is.na(median))
range(all_covid_rt_avg$date)
#
head(all_covid_rt)
# "2020-02-25" "2022-06-30"
#ours: "2020-02-08" "2020-12-24"
all_covid_rt %>%
  filter(date >= as.Date("2020-01-01")) %>%
  mutate(
    lower_perc = 100 * (lower - median) / median,
    upper_perc = 100 * (upper - median) / median
  )
head(all_covid_rt)
# by March 1, 2020, lower and upper prediction intervals are 13% below and 15% above median Rt
#
# add OSI data here
#
combined_mob <- read_csv("/Users///COVID-Counterfactual/all_OSI.csv") %>% dplyr::select(location, date, stringency_index)
combined_mob$date <- dmy(combined_mob$date)
#combined_mob$epi_date <- as.Date(combined_mob$epi_date)
range(combined_mob$date)
#
head(combined_mob)
# plot of all country's OSI for 2021
#
#combine the 2 df for ccf calculations not for plots below
#
all_covid_rt_avg <- rename(all_covid_rt_avg, location = group)
# Combine dataframes 
combined_rt_osi <- merge(combined_mob, all_covid_rt_avg, by = c("date", "location"))
head(combined_rt_osi)
#
saveRDS(combined_rt_osi,"/Users///COVID-Counterfactual/Top5/SI/rt_osi_combine_top5.rds")
#
#######################################################################################
## Combine mobility and Rt data into one data frame for statistical analyses
########################################################################################
#combined_rt_osi <- read_rds("/Users///COVID-Counterfactual/2021/top3/rt_osi_combine_top3.rds")
combined_mob_avg <- combined_rt_osi %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))
#
combined <- rename(combined_mob_avg, COVID2 = median)
#create epi_date, no idea what it is needed for ts_block
combined$epi_date <- combined$date
#choose country
combined_Malaysia <- combined %>% filter(location %in% c("Malaysia"))
###########################################
## Actual CCF Loop
###########################################
metrics <- c( "stringency_index")
#this will take stringency and place it on scale where median = 0
comb_weekly <- combined_Malaysia %>%
  mutate_at(vars(all_of(metrics)), ~ scale(.x) %>% as.vector()) %>%
  filter(epi_date >= as.Date("2020-01-01") & epi_date < as.Date("2020-12-31")) %>%
  group_by(epi_date) %>%
  summarize_at(c("COVID2",metrics), ~ mean(.x, na.rm = T))
#
actual_data_and_perm <- ts_block_bootstrap_sp(df1=comb_weekly,pathogen="COVID2",window=10,l=0.077,perms=100)
#
save(actual_data_and_perm, file = "/Users///COVID-Counterfactual/Top5/SI/Malaysia_mobility_CCF_5mo_sliding_window_actual_and_null_output_osi_spearman.RData")
#
########################################################################################
# Compile 5mo block bootstrap results
# covid : need spearman OSI
load(paste0("/Users///COVID-Counterfactual/Top5/SI/Malaysia_mobility_CCF_5mo_sliding_window_actual_and_null_output_osi_spearman.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("stringency_index"))
all_covid <- bind_rows(actual_data_and_perm3)
unique(all_covid$mobility_metric)
unique(all_covid$pathogen)
########################################################################################
## combine data for all pathogens
########################################################################################
all_path_ccf <- all_covid
all_path_ccf$pathogen <- as.factor(all_path_ccf$pathogen)
levels(all_path_ccf$pathogen)
levels(all_path_ccf$mobility_metric)
unique(all_path_ccf$pathogen)
all_path_ccf$pathogen <- factor(all_path_ccf$pathogen, levels = c("COVID2"))
levels(all_path_ccf$pathogen)
levels(all_path_ccf$mobility_metric)
########################################################################################
## Import 5mo block bootstrap results
########################################################################################
#reread Rt data
all_endemic_results_avg <- combined
all_endemic_results_avg$organism <- as.factor(all_endemic_results_avg$organism)
levels(all_endemic_results_avg$organism)
rt_weekly <- all_endemic_results_avg %>%
  group_by(epi_date, tag, level, organism, location) %>%
  summarize_at(c("lower", "upper", "COVID2"), ~ mean(.x, na.rm = T)) %>%
  arrange(organism, epi_date, location) %>%
  droplevels()
#
unique(rt_weekly$organism)
levels(all_path_ccf$pathogen)
levels(rt_weekly$organism)
levels(all_path_ccf$pathogen)
levels(rt_weekly$organism) <- levels(all_path_ccf$pathogen)
##############
# PLOTS
########################################################################################
#Figure 5: Cross-correlations between Rt and mobility across SARS-CoV-2 waves
#######################################################################################
unique(all_path_ccf$mobility_metric)
#
covid_mob_red5 <- ggplot(
  all_path_ccf %>% filter(start_week >= as.Date("2020-01-01") & start_week < as.Date("2020-12-31") & pathogen == "COVID2" &
                            mobility_metric %in% c(
                              "stringency_index"
                            )),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  # facet_grid(mobility_metric ~ location) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  scale_alpha_manual(values = c(0.1, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 18) +
  ylab("Cross-Correlation Coefficient: Rt and Stringency Index") +
  # scale_fill_viridis_c() +
  scale_fill_viridis_c(begin = 0, end = 0.75) +
  xlab("5-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Weeks)") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )
save_plot(covid_mob_red5, file = "/Users///COVID-Counterfactual/Top5/SI/Malaysia_fig_5_block_bootstrap_covid_4mo_moving_window_select_mobility_indicators_spearman_neg_lags.png", base_width = 14, base_height = 14)
