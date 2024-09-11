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
library(epidemia)
source("/Users///COVID-Counterfactual/utils.R")
source("/Users///COVID-Counterfactual/get_quantiles.R")
###################################################
###################################################
##### COVID-19
###################################################
fit=readRDS("/Users///COVID-Counterfactual/2021/next5/fit_2021_next5_iter300.rds")
all_covid_rt = extract_Rt(fit)
all_covid_rt$organism <- "SARS-CoV-2"
head(all_covid_rt)
#
all_covid_rt_avg <- all_covid_rt %>%
  group_by(group, tag, level, organism) %>%
  mutate_at(c("lower", "upper", "median"), ~ zoo::rollmean(., k = 15, align = "center", fill = NA)) %>%
  ungroup() %>%
  filter(date < as.Date("2022-07-01")) %>%
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
combined_mob <- read_csv("/Users///COVID-Counterfactual/all_OSI.csv") %>% dplyr::select(location, date, stringency_index)
combined_mob$date <- dmy(combined_mob$date)
#combined_mob$epi_date <- as.Date(combined_mob$epi_date)
range(combined_mob$date)
#
head(combined_mob)
# plot of all country's OSI for 2021
#
ggplot_all_osi <- ggplot() +  
  geom_line(data = combined_mob %>% filter(date >= as.Date("2020-01-01") & date < as.Date("2021-12-31")),
            aes(x = date, y = stringency_index, color = location), alpha = 0.8, lwd = 2
  ) +
  labs(title = "Stringency Index Over Time by Location",
       x = "Date",
       y = "Stringency Index") +
  theme_bw(base_size = 16)
save_plot(ggplot_all_osi,, filename = "/Users///COVID-Counterfactual/2021/all_OSI.png", base_width = 16, base_height = 16)
#
#combine the 2 df for ccf calculations not for plots below
#
all_covid_rt_avg <- rename(all_covid_rt_avg, location = group)
# Combine dataframes 
combined_rt_osi <- merge(combined_mob, all_covid_rt_avg, by = c("date", "location"))
head(combined_rt_osi)
#
saveRDS(combined_rt_osi,"/Users///COVID-Counterfactual/2021/next5/rt_osi_combine_next5.rds")
#
# Create the facetted line plot
p <- ggplot() +
  geom_line(data = combined_rt_osi %>% filter(date >= as.Date("2020-01-01") & date < as.Date("2021-12-31")),
            aes(x = date, y = median, group = location), linewidth = 0.5
  ) +
  geom_line(data = combined_rt_osi %>% filter(date >= as.Date("2020-01-01") & date < as.Date("2021-12-31")),
            aes(x = date, y = stringency_index /50, group = location, colour = "stringency_index"), alpha = 0.8, lwd = 0.5
  ) +
#  scale_color_manual(values = c("stringency_index" = "green")) +
  xlab("Date") +
  ylab("Rt") +
  scale_x_date(expand = c(0, 0), date_breaks = "2 months", date_labels = "%b %Y") +
  facet_wrap(~ location) +  # Adjust ncol as needed
# You could set the x axis in an 90° angle to get a cleaner plot
theme(axis.text.x = element_text(angle = 90,
                                 vjust = 0.5,
                                 hjust = 1))
# Save the plot
ggsave("/Users///COVID-Counterfactual/2021/all_locations_line_plot.png", p)
#
ggplot_correlation <- ggplot() +
  geom_line(data = combined_rt_osi %>% filter(date >= as.Date("2020-01-01") & date < as.Date("2021-12-31")),
            aes(x = date, y = median, color = location), alpha = 0.8, lwd = 2
  ) +
  geom_line(data = combined_rt_osi %>% filter(date >= as.Date("2020-01-01") & date < as.Date("2021-12-31")),
            aes(x = date, y = stringency_index, color = location), alpha = 0.8, lwd = 2
  ) +
  theme_bw(base_size = 16)
save_plot(ggplot_correlation,, filename = "/Users///COVID-Counterfactual/2021/rt_OSI_correlation.png", base_width = 16, base_height = 16)
#
# change country below
#
all_covid_plot <- ggplot() +
  geom_vline(xintercept = as.Date("2021-01-01"), lty = "dashed", color = "darkgreen") +
  geom_line(
    data = combined_mob %>% filter(location=="Vietnam" & date >= as.Date("2021-01-01") & date < as.Date("2021-12-31")),
#scale OSI to get on same change secondry axes as well
    aes(x = date, y = stringency_index / 50, color = "OSI"), alpha = 0.8, lwd = 2
  ) +
  geom_line(
    data = all_covid_rt_avg %>%
      dplyr::select(date,location,median) %>% distinct() %>% filter(location == "Vietnam" & date >= as.Date("2021-01-01") & date < as.Date("2021-12-31")),
    aes(x = date, y = median, col = "All"), lwd = 1.2, lty = "solid"
  ) +
  theme_bw(base_size = 16) +
  # geom_hline(aes(yintercept = 1), lty = "dashed", color = "#004488", lwd = 1) +
   scale_y_continuous(
    # Features of the first axis
    name = "Rt",
    # Add a second axis and specify its features match above
    sec.axis = sec_axis(~ . * 50, name = "OSI"),
    #second number expands above ymax
   # expand = c(0.5, -0.5)
   limits = c(NA, 3)
  ) +
   scale_fill_manual(values = c("#ABB065", "#DB9D85"), name = NULL, breaks = c("All", "OSI"), labels = c("SARS-CoV-2 Rt", "Stringency Index")) +
   scale_color_manual(values = c("#ABB065", "#DB9D85"), name = NULL, breaks = c("All", "OSI"), labels = c("SARS-CoV-2 Rt", "Stringency Index")) +
   guides(fill = "none") +
  xlab("date") +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b %Y") +
  theme(legend.position = c(0.15, 0.9), legend.text = element_text(size = 22), legend.background = element_blank())
all_covid_plot

end_and_cov <- plot_grid(all_covid_plot + ggtitle("Vietnam SARS-CoV-2: Rt and Stringency Index") + theme(plot.title = element_text(hjust = 0.5)),
                         nrow = 1,
                         labels = "AUTO"
)
end_and_cov
#
save_plot(end_and_cov,, filename = "/Users///COVID-Counterfactual/2021/next5/Vietnam_fig_1_rt_endemic_and_covid_combined.png", base_width = 16, base_height = 16)
#
