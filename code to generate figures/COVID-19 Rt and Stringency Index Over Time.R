## Input data are provided in the Epidemia_Models/epidemia_rt_output folder, so it's not necessary to run the Epidemia models to run this code
###################################################
## Plot Rt and incidence (Figure 1b)
###################################################
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(jcolors)
library(tidyquant)
library(timetk)
source("/Users///COVID-Counterfactual/utils.R")
source("/Users///COVID-Counterfactual/get_quantiles.R")
###################################################
###################################################
##### COVID-19
###################################################
fit=readRDS("/Users///COVID-Counterfactual/Top5/2021_SI/fit_2021_top5_iterations200_SI.rds")

all_covid_rt = extract_Rt(fit)
all_covid_rt$organism <- "SARS-CoV-2"
#filter for country
#TOP 5 ("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
country <- "Singapore"
all_covid_rt2 <- all_covid_rt %>% filter(group == country)
#
all_covid_rt_avg <- all_covid_rt2 %>%
  group_by(group, tag, level, organism) %>%
  mutate_at(c("lower", "upper", "median"), ~ zoo::rollmean(., k = 15, align = "center", fill = NA)) %>%
  ungroup() %>%
  filter(date < as.Date("2021-12-31")) %>%
  filter(!is.na(median))
all_covid_rt2 %>%
  filter(date >= as.Date("2021-12-31")) %>%
  mutate(
    lower_perc = 100 * (lower - median) / median,
    upper_perc = 100 * (upper - median) / median
  )
# by March 1, 2020, lower and upper prediction intervals are 13% below and 15% above median Rt
#
# add OSI data here
combined_mob <- read_csv("/Users///COVID-Counterfactual/all_OSI.csv",show_col_types = FALSE) %>% dplyr::select(location, date, stringency_index)
#range(combined_mob$date)
#reformat date column
combined_mob <- combined_mob %>%
  mutate(date = as.Date(date, format = '%d/%m/%Y')) %>%
  rename(group = location)
#filter for country
combined_mob <- combined_mob %>% filter(group == country)
###################
# data plotted here
#
all_covid_plot <-  ggplot() +
#SI
    geom_line(data = combined_mob %>% filter(date >= as.Date("2021-02-01") & date <= as.Date("2021-12-31")),
            aes(x = date, y = stringency_index / 100 ), color = "blue", alpha = 0.8, lwd = 5, group = "group") +
# Rt
 geom_line(data = all_covid_rt_avg %>%
    dplyr::select(date, median) %>% distinct() %>% filter(date >= as.Date("2021-02-01") & date < as.Date("2021-12-31")),
  aes(x = date, y = median, col = "black"), lwd = 2, lty = "solid",group = "group"
) +
#confidence for Rt
geom_ribbon(
  data = all_covid_rt_avg %>%
    dplyr::select(date, level, lower, upper) %>%
    filter(level == 90) %>%
    distinct() %>% filter(date >= as.Date("2021-02-01") & date < as.Date("2021-12-31")),
  aes(x = date, ymin = lower, ymax = upper, fill = "All",group = "group"), alpha = 0.5
) +
  theme_bw(base_size = 16) +
  geom_hline(aes(yintercept = 1), lty = "dashed", color = "#004488", lwd = 1) +
scale_y_continuous(
  # Features of the first axis
  name = "Rt",
  limits = c(0, 6),
  # Add a second axis and specify its features
  sec.axis = sec_axis(~ . * 100, name = "Stringency Index"),
  expand = c(0, 0.1)
) +
scale_fill_manual(values = c("#ABB065", "#DB9D85"), name = NULL, breaks = c("All", "OSI"), labels = c("SARS-CoV-2 Rt", "Oxford Stringency Index")) +
  scale_color_manual(values = c("#ABB065", "#DB9D85"), name = NULL, breaks = c("All", "OSI"), labels = c("SARS-CoV-2 Rt", "Oxford Stringency Index")) +
  guides(fill = "none") +
  xlab("date") +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b %Y") +
  theme(legend.position = c(0.65, 0.9), legend.background = element_blank()) +
  labs(title = paste(country,": COVID-19 Rt and Stringency Index Over Time (2021)") )+
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold"))
#
save_plot(all_covid_plot, filename = (paste0("/Users///COVID-Counterfactual/Top5/2021_SI/", country, "_fig_1_rt_endemic_and_covid_combined.png")), base_width = 16, base_height = 16)
#

