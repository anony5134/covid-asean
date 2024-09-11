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
fit_SI=readRDS("/Users///COVID-Counterfactual/Top5/SI/fit_2020_top5_iter200_SI.rds")
fit_noSI= readRDS("/Users///COVID-Counterfactual/Top5/noSI/fit_2020_top5_iter200_noSI.rds")
# save Rt
orig_rt_SI = extract_Rt(fit_SI)
orig_rt_noSI = extract_Rt(fit_noSI)
#
#TOP 5 ("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
#
country <- "Thailand"
orig_rt_SI2 <- orig_rt_SI %>% filter(group == country)
orig_rt_noSI2 <- orig_rt_noSI %>% filter(group == country)
########
# Plot here
#
Title = paste(country,": Rt modelled with (blue) and without (red) Stringency Index Over Time (2020)")
all_rt_plot <-  ggplot() +
  geom_line(data = orig_rt_SI2 %>%
            dplyr::select(date, median) %>% distinct() %>% filter(date >= as.Date("2020-02-01") & date < as.Date("2020-12-31")),
            aes(x = date, y = median, col = "blue"), lwd = 3, lty = "solid",group = "group"
  ) +
  #confidence for Rt
  geom_ribbon(
    data = orig_rt_SI2 %>%
      dplyr::select(date, level, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% filter(date >= as.Date("2020-02-01") & date < as.Date("2021-12-31")),
    aes(x = date, ymin = lower, ymax = upper, fill = "blue",group = "group"), alpha = 0.2
  ) +
geom_line(data = orig_rt_noSI2 %>%
            dplyr::select(date, median) %>% distinct() %>% filter(date >= as.Date("2020-02-01") & date < as.Date("2020-12-31")),
          aes(x = date, y = median, col = "red"), lwd = 3, lty = "solid",group = "group"
) +
#confidence for Rt
geom_ribbon(
  data = orig_rt_noSI2 %>%
    dplyr::select(date, level, lower, upper) %>%
    filter(level == 90) %>%
    distinct() %>% filter(date >= as.Date("2020-02-01") & date < as.Date("2021-12-31")),
  aes(x = date, ymin = lower, ymax = upper, fill = "red",group = "group"), alpha = 0.2
) +
  theme_bw(base_size = 16) +
  geom_hline(aes(yintercept = 1), lty = "dashed", color = "#004488", lwd = 1) +
  xlab("date") +
  ylab("Rt") +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b %Y") +
  theme(legend.position = "none") +
  labs(title = str_wrap(Title, 60)) +
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold", lineheight = 1.2))
#
save_plot(all_rt_plot, filename = (paste0("/Users///COVID-Counterfactual/Top5/", country, "_fig_2_rt.png")), base_width = 16, base_height = 16)

