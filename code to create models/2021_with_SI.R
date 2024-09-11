#install.packages("/Users///COVID-Counterfactual/StanHeaders_2.21.0-7.tar.gz", repos = NULL, type = "source")
#install.packages("~//COVID-Counterfactual/rstan_2.21.2.tar.gz", repos = NULL, type = "source")
#install.packages("~//COVID-Counterfactual/rstanarm_2.21.1.tar.gz", repos = NULL, type = "source")
#devtools::install_github('ImperialCollegeLondon/epidemia', ref='exponential_new', force = TRUE)
library(here)
library(epidemia)
library(tidyverse)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(scales)
library(EnvStats)
library(gridExtra)
source("/Users///COVID-Counterfactual/Top5/SI/plot_with_eta.R")
library(dplyr)
library(arm)
#TOP 5 ("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
iterate = "200"
### ASEAN COV DEATH DATA
data_world_confirmed = read.csv(url("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"))
data_world_recovered = read.csv(url("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"))
data_world_deaths = read.csv(url("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"))

name_countries = unique(data_world_confirmed$Country.Region)
flag1=which(name_countries=="Cambodia")
flag2=which(name_countries=="Malaysia")
flag3=which(name_countries=="Philippines")
flag4=which(name_countries=="Singapore")
flag5=which(name_countries=="Thailand")
flag_vec=c(flag1,flag2,flag3,flag4,flag5)

country_name = unique(data_world_deaths$Country.Region)[flag1]
data_country_deaths = colSums((data_world_deaths %>% filter(Country.Region == country_name))[, 5:ncol(data_world_deaths)])
dates = seq(as.Date("2021-01-22"), as.Date("2021-01-22") + ncol(data_world_deaths)-4, by = 1) ## yyyy-mm-dd
date_initial = as.Date("2021-02-01")
date_final = as.Date("2021-12-31")
#death totals
I = abs(diff(data_country_deaths[(which(dates == date_initial)-1):(which(dates == date_final))]))
dates=as.Date(dates[(which(dates == date_initial)):(which(dates == date_final))])
data=data.frame(dates)
data_cases=data.frame(dates)
data_total_cases=data.frame(dates)
for(flag in flag_vec){
  country_name = unique(data_world_deaths$Country.Region)[flag]
  data_country_deaths = colSums((data_world_deaths %>% filter(Country.Region == country_name))[, 5:ncol(data_world_deaths)])
  data_country_cases = colSums((data_world_confirmed %>% filter(Country.Region == country_name))[, 5:ncol(data_world_confirmed)])
  dates = seq(as.Date("2021-01-22"), as.Date("2021-01-22") + ncol(data_world_deaths)-4, by = 1) ## yyyy-mm-dd
  date_initial = as.Date("2021-02-01")
  date_final = as.Date("2021-12-31")
  TC=data_country_cases[(which(dates == date_initial)):(which(dates == date_final))]
  #total deaths
  I = abs(diff(data_country_deaths[(which(dates == date_initial)-1):(which(dates == date_final))]))
  C = abs(diff(data_country_cases[(which(dates == date_initial)-1):(which(dates == date_final))]))
  dates=as.Date(dates[(which(dates == date_initial)):(which(dates == date_final))])
  #data now contains I
  data=cbind.data.frame(data,I)
  data_cases=cbind.data.frame(data_cases,C)
  data_total_cases=cbind.data.frame(data_total_cases,TC)
}
colnames(data)=c("Date","Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
rownames(data)=c(1:nrow(data))
colnames(data_cases)=c("Date","Cambodia", "Malaysia","Philippines", "Singapore", "Thailand")
rownames(data_cases)=c(1:nrow(data_cases))
colnames(data_total_cases)=c("Date", "Cambodia", "Malaysia","Philippines", "Singapore", "Thailand")
rownames(data_total_cases)=c(1:nrow(data_total_cases))
dat<-data
#save data
saveRDS(dat,"/Users///COVID-Counterfactual/Top5/2021_SI/data_2021_top5.rds")
#saveRDS(dat,paste0("/Users///COVID-Counterfactual/Top5/2021_SI/2021_SI/2021_SI/data_2021_top5_iterations",iterate,"_SI.rds"))
#
#
# onset to door time function
# onsettime as the time of the symptom onset based on the patientinterview. 
# Door time was when the patient presented to thehospital emergency department 
#
o2d<-function(){
  i2o <- EuropeCovid$obs$deaths$i2o
  shape1 <- 5.807; scale1 <- 0.948; # infection to onset https://www.acpjournals.org/doi/10.7326/M20-0504
  shape2 <- 1.454 ; scale2 <- 10.434 # using estimated of onset to death from chess data
  x1 <- rgamma(1e6,shape=shape1,scale=scale1) # infection-to-onset distribution
  x2 <- rgamma(1e6,shape=shape2,scale=scale2) # infection-to-onset distribution
  f_cached <- ecdf(x1+x2) # empirical cumulative distribtion function
  convolution = function(u) f_cached(u)
  f = rep(NA,length(EuropeCovid$obs$deaths$i2o)) # f is the probability of dying on day i given infection
  f[1] = (convolution(1.5) - convolution(0)) # first entry
  for(i in 2:length(EuropeCovid$obs$deaths$i2o)) { # all other entries
    f[i] = (convolution(i+.5) - convolution(i-.5))
  }
  return(f)
}
#
dat <- pivot_longer(dat, cols = -Date, names_to = 'country', values_to = 'deaths')
dat$deaths[is.na(dat$deaths)]=0 # to have all time series starting on the same day 30th Jan 2020
dat <- drop_na(dat) # to drop NAs just incase
dat <- mutate(dat, Date = as.Date(Date, format='%d/%m/%Y'))
dat <- rename(dat, date=Date)
dat <- filter(dat, date <= as.Date("2021-12-31")) %>% arrange(country, date) # only use data to July
#
#Add oxford stringency index data SI data to this frame
#
SI <- read_csv("/Users///COVID-Counterfactual/all_OSI.csv",show_col_types = FALSE) %>% dplyr::select(location, date, stringency_index)
# Combine dataframes 
colnames(SI) <- c("country", "date","SI")
#reformat date column
SI <- SI %>%
  mutate(date = as.Date(date, format = '%d/%m/%Y')) %>%
  mutate(date = format(date, "%Y-%m-%d"))
SI <- SI %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))
# Merge and create new column
dat_osi <- dat %>%
  left_join(SI, by = c("date", "country")) %>%
  mutate(SInd = coalesce(SI, 0)) 
dat2 <- dat_osi
#
i2o <- o2d() # note choice here to model observations using o2d
# does not sum to 1 ??
# deaths are related to infections by using infection to death distribution and ifr
#In order to infer the effects of control measures on transmission, we
# must fit the model to data. Here, daily deaths are used. In theory,
# additional types of data can be included in the model, but such
# extension are not considered here. A simple intercept model is used for
# the infection fatality rate (IFR). 
#This makes the assumption that the
# IFR is constant over time. The model can be written as follows.
#' By using `link = scaled_logit(0.02)`, we let the IFR range between
#' $0\%$ and $2\%$. 'link' must be one of logit, probit, cauchit, cloglog, identity
#' 	link: A string representing the link function used to transform the linear predictor
#' In conjunction with the symmetric prior on the
#' intercept, this gives the IFR a prior mean of $1\%$.
#' `EuropeCovid2$inf2death` is a numeric vector giving the same
#' distribution for the time from infection to death as that used in
#' @Flaxman2020.
# deaths <- epiobs(
#   formula = deaths(country, date) ~ 1,
#   prior_intercept = rstanarm::normal(log(2) ,scale = 0.03),
#   prior_aux = rstanarm::normal(10), scale=2,
#   link="logit",
#   i2o=i2o * 0.01
# )
# deaths are related to infections by using infection to death distribution and ifr
deaths <- epiobs(
  formula = deaths(country, date) ~ 1,
  prior_intercept = rstanarm::normal(location=0.05,scale = 0.03),
  prior_aux = rstanarm::normal(location=10, scale=2),
  link="identity",
  i2o=i2o * 0.01,
)
#'family' must be one of poisson, neg_binom, quasi_poisson
# sampling parameters of chains, iterations, seed and samling related parameters
args<-NULL
args$data <- dat2
args$algorithm <- "sampling"
args$obs <- list(deaths=deaths)
#
#initial run??
args$init_run = T
#iter=1500 original takes hours, 50 and 100 failed
#
#, delta 0.95-> 0.99
args$sampling_args <- list(iter=200, seed=12345, control=list(adapt_delta=0.99,max_treedepth=15))
# for rt just using a weekly random walk for each country with a separate R0
args$rt <- epirt(
  #  formula = R(country, date) ~ 0 + country + rw(time = week, gr=country, prior_scale = 0.1),
  formula = R(country, date) ~ 0 + country + (SInd || country) + rw(time = week, gr=country, prior_scale = 0.1),
  prior = shifted_gamma(shape = 1/2, scale = 1, shift = log(1.05)/2),
  prior_intercept = rstanarm::normal(scale=0.9)
)
args$prior_tau = rstanarm::exponential(rate = 1)
# make sure for peridod before week of 13th March has same weekly phiex (same Rt)
w=rep(c(1:(round(nrow(data))/7+1)),each=7)
w=w[c(1:nrow(data))]
#
#change this to match number of countries
#
args$data$week <- as.character(rep(w,times=5))
nrow(data)
# start random walk on a given date
#args$data <- mutate(
#args$data,
#week = replace(week, which(week <= 11), NA)
#)
country=c("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
pop=c(1.677e7, 3.394e7, 1.1556e8, 5.64e6, 7.170e7)
popu=cbind.data.frame(country,pop)
args$pop=popu
args$pop_adjust <- F
args$si=EuropeCovid$si
#
#takes ALOT OF TIME:
#
#fit <- do.call(epim, args)
#
#
#saveRDS(fit,paste0("/Users///COVID-Counterfactual/Top5/2021_SI/fit_2021_top5_iterations",iterate,"_SI.rds"))
#plot_obs(fit, type = "deaths", levels = c(50, 95))
fit=readRDS("/Users///COVID-Counterfactual/Top5/2021_SI/fit_2021_top5_iterations200_SI.rds")
# original fit
fit$pop_adjust <- TRUE
fit_orig <- fit
# save Rt
orig_rt = plot_rt(fit,date_breaks = "1 month")
ggsave(here("/Users///COVID-Counterfactual/Top5/2021_SI/SI_original-rt.pdf"),orig_rt,width=11,height=4)
# date from which we will change Rt
changeDate <- as.Date('2021-03-15')
# get phiex for random walks
nms <- colnames(as.matrix(fit_orig))
idx_cam<- grep("^R\\|rw.*Cambodia", nms)
idx_mal <- grep("^R\\|rw.*Malaysia", nms)
idx_phi <- grep("^R\\|rw.*Philippines", nms)
idx_sin <- grep("^R\\|rw.*Singapore", nms)
idx_tha <- grep("^R\\|rw.*Thailand", nms)
#get phiex for R0
idx_cam_R0 <- grep("^R\\|countryCambodia", nms)
idx_mal_R0 <- grep("^R\\|countryMalaysia", nms)
idx_phi_R0 <- grep("^R\\|countryPhilippines", nms)
idx_sin_R0 <- grep("^R\\|countrySingapore", nms)
idx_tha_R0 <- grep("^R\\|countryThailand", nms)
#get phiex for seeds
idx_cam_seeds <- grep("seeds\\[Cambodia\\]", nms)
idx_mal_seeds <- grep("seeds\\[Malaysia\\]", nms)
idx_phi_seeds <- grep("seeds\\[Philippines\\]", nms)
idx_sin_seeds <- grep("seeds\\[Singapore\\]", nms)
idx_tha_seeds <- grep("seeds\\[Thailand\\]", nms)
#
nchains <- length(fit_orig$stanfit@sim$samples)
#
mat <- as.matrix(fit)
# compute orderings for the draws
order_cam <- order(mat[, idx_cam_R0])
order_mal <- order(mat[, idx_mal_R0])
order_phi <- order(mat[, idx_phi_R0])
order_sin <- order(mat[, idx_sin_R0])
order_tha <- order(mat[, idx_tha_R0])
#
camSeeds = mat[,idx_cam_seeds]
malSeeds = mat[,idx_mal_seeds]
phiSeeds = mat[,idx_phi_seeds]
sinSeeds = mat[,idx_sin_seeds]
thaSeeds = mat[,idx_tha_seeds]
#
camR=exp(mat[,idx_cam_R0])
malR=exp(mat[,idx_mal_R0])
phiR=exp(mat[,idx_phi_R0])
sinR=exp(mat[,idx_sin_R0])
thaR=exp(mat[,idx_tha_R0])
#Calculating the summary statistics
camsum <- summary(camR,digits=2)
#####################################################################
pdf(here('/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/R0_ratios.pdf'))
par(mfrow=c(3,3))
hist(camR,main="R0 Cambodia")
hist(malR,main="R0 Malaysia")
hist(phiR,main="R0 Philippines")
hist(sinR,main="R0 Singapore")
hist(thaR,main="R0 Thailand")
#
hist(camR/malR,main="R0 ratio Cambodia/Malaysia")
hist(malR/camR,main="R0 ratio Malaysia/Cambodia")
hist(camR/phiR,main="R0 ratio Cambodia/Philippines")
hist(phiR/camR,main="R0 ratio Philippines/Cambodia")
hist(camR/sinR,main="R0 ratio Cambodia/Singapore")
hist(sinR/camR,main="R0 ratio Singapore/Cambodia")
hist(camR/thaR,main="R0 ratio Cambodia/Thailand")
hist(thaR/camR,main="R0 ratio Thailand/Cambodia")
#
hist(malR/phiR,main="R0 ratio Malaysia/Philippines")
hist(phiR/malR,main="R0 ratio Philippines/Malaysia")
hist(malR/sinR,main="R0 ratio Malaysia/Singapore")
hist(sinR/malR,main="R0 ratio Singapore/Malaysia")
hist(malR/thaR,main="R0 ratio Malaysia/Thailand")
hist(thaR/malR,main="R0 ratio Thailand/Malaysia")
#
hist(phiR/sinR,main="R0 ratio Philippines/Singapore")
hist(sinR/phiR,main="R0 ratio Singapore/Philippines")
hist(phiR/thaR,main="R0 ratio Philippines/Thailand")
hist(thaR/phiR,main="R0 ratio Thailand/Philippines")
#
hist(sinR/thaR,main="R0 ratio Singapore/Thailand")
hist(thaR/sinR,main="R0 ratio Thailand/Singapore")
#
dev.off()
#####################################################################
pdf(here('/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/boxplotR0_top5.pdf'))
bp <- boxplot(camR,malR, phiR, sinR, thaR, col = c("gold", "darkgreen"), xaxt = "n", ylab="R0 (average number of secondary infections)" )
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
z=c("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
text(tick, par("usr")[3] - 0.1, z, srt = 45, xpd = TRUE)
dev.off()
#####################################################################
pdf(here('/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/PairPlot_R0_seeds.pdf'))
par(mfrow=c(1,3))
plot(camR,camSeeds,pch=16,main="Cambodia",xlab="R0",ylab="Seeds")
plot(malR,malSeeds,pch=16,main="Malaysia",xlab="R0",ylab="Seeds")
plot(phiR,malSeeds,pch=16,main="Philippines",xlab="R0",ylab="Seeds")
plot(sinR,sinSeeds,pch=16,main="Singapore",xlab="R0",ylab="Seeds")
plot(thaR,thaSeeds,pch=16,main="Thailand",xlab="R0",ylab="Seeds")

dev.off()
#
#####################################################################
pdf(here('/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/R0_prior.pdf'))
hist(exp(rnorm(1e6,log(3.5),0.1)),100,main="Prior distribution on R0",xlab="R0",ylab="Density",col='red')
dev.off()

##########################
# generate counterfactuals
#z=c("Cambodia", "Malaysia", "Philippines", "Singapore", "Thailand")
# generate counterfactuals

e_orig <- posterior_linpred(fit)
e_orig$group
## absolute
e1 <- e_orig
w <- e_orig$time >= changeDate

e1$draws[order_cam, w & (e1$group == "Cambodia")] <- e_orig$draws[order_mal, w & (e_orig$group == "Malaysia")]
e1$draws[order_mal, w & (e1$group == "Malaysia")] <- e_orig$draws[order_phi, w & (e_orig$group == "Philippines")]
e1$draws[order_phi, w & (e1$group == "Philippines")] <- e_orig$draws[order_sin, w & (e_orig$group == "Singapore")]
e1$draws[order_sin, w & (e1$group == "Singapore")] <- e_orig$draws[order_tha, w & (e_orig$group == "Thailand")]
e1$draws[order_tha, w & (e1$group == "Thailand")] <- e_orig$draws[order_cam, w & (e_orig$group == "Cambodia")]

## absolute
e2 <- e_orig
w <- e_orig$time >= changeDate

e2$draws[order_mal, w & (e2$group == "Malaysia")] <- e_orig$draws[order_cam, w & (e_orig$group == "Cambodia")]
e2$draws[order_phi, w & (e2$group == "Philippines")] <- e_orig$draws[order_mal, w & (e_orig$group == "Malaysia")]
e2$draws[order_sin, w & (e2$group == "Singapore")] <- e_orig$draws[order_phi, w & (e_orig$group == "Philippines")]
e2$draws[order_tha, w & (e2$group == "Thailand")] <- e_orig$draws[order_sin, w & (e_orig$group == "Singapore")]
e2$draws[order_cam, w & (e2$group == "Cambodia")] <- e_orig$draws[order_tha, w & (e_orig$group == "Thailand")]

## absolute
e3 <- e_orig
w <- e_orig$time >= changeDate

e3$draws[order_cam, w & (e3$group == "Cambodia")] <- e_orig$draws[order_phi, w & (e_orig$group == "Philippines")]
e3$draws[order_phi, w & (e3$group == "Philippines")] <- e_orig$draws[order_tha, w & (e_orig$group == "Thailand")]
e3$draws[order_tha, w & (e3$group == "Thailand")] <- e_orig$draws[order_mal, w & (e_orig$group == "Malaysia")]
e3$draws[order_mal, w & (e3$group == "Malaysia")] <- e_orig$draws[order_sin, w & (e_orig$group == "Singapore")]
e3$draws[order_sin, w & (e3$group == "Singapore")] <- e_orig$draws[order_cam, w & (e_orig$group == "Cambodia")]

## absolute
e4 <- e_orig
w <- e_orig$time >= changeDate

e4$draws[order_phi, w & (e4$group == "Philippines")] <- e_orig$draws[order_cam, w & (e_orig$group == "Cambodia")]
e4$draws[order_cam, w & (e4$group == "Cambodia")] <- e_orig$draws[order_sin, w & (e_orig$group == "Singapore")]
e4$draws[order_sin, w & (e4$group == "Singapore")] <- e_orig$draws[order_mal, w & (e_orig$group == "Malaysia")]
e4$draws[order_mal, w & (e4$group == "Malaysia")] <- e_orig$draws[order_tha, w & (e_orig$group == "Thailand")]
e4$draws[order_tha, w & (e4$group == "Thailand")] <- e_orig$draws[order_phi, w & (e_orig$group == "Philippines")]

## relative
fit5=fit_orig
warmup <- args$sampling_args$iter/2
for (chain in 1:nchains) {
  order_cam <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_cam_R0]][-(1:warmup)])
  order_mal<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_mal_R0]][-(1:warmup)])
  order_phi <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_phi_R0]][-(1:warmup)])
  order_sin<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_sin_R0]][-(1:warmup)])
  order_tha <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_tha_R0]][-(1:warmup)])
  for (i in 1:length(idx_cam)) {
    fit5$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam] <- fit_orig$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal]
    fit5$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal] <- fit_orig$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi]
    fit5$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi] <- fit_orig$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin]
    fit5$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin] <- fit_orig$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha]
    fit5$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha] <- fit_orig$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam]
  }
}


#par(mfrow=c(2,2))
#plot(phiR,banR,pch=16,xlab='R0 Cambodia', ylab='R0 Malaysia', main="Posterior Cambodia R0 against Malaysia R0")
#temp_samples=as.matrix(fit_orig)
#plot(temp_samples[,idx_cam_R0],temp_samples[,idx_cam[1]],pch=16,xlab="R0 Cambodia",ylab="Cambodia first week random walk",main='Original samples')
#temp_samples=as.matrix(fit_corr_plot)
#plot(temp_samples[,idx_cam_R0],temp_samples[,idx_cam[1]],pch=16,xlab="R0 Cambodia",ylab="Cambodia first week random walk",main='Naive relative approach')
#temp_samples=as.matrix(fit3)
#plot(temp_samples[,idx_cam_R0],temp_samples[,idx_cam[1]],pch=16,xlab="R0 Cambodia",ylab="Cambodia first week random walk",main='Ordered relative approach')
#dev.off()

#relative
fit6=fit_orig
warmup <- args$sampling_args$iter/2
for (chain in 1:nchains) {
  order_cam <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_cam_R0]][-(1:warmup)])
  order_mal<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_mal_R0]][-(1:warmup)])
  order_phi <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_phi_R0]][-(1:warmup)])
  order_sin<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_sin_R0]][-(1:warmup)])
  order_tha <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_tha_R0]][-(1:warmup)])
  for (i in 1:length(idx_cam)) {
    fit6$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam] <- fit_orig$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha]
    fit6$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha] <- fit_orig$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin]
    fit6$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal] <- fit_orig$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam]
    fit6$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi] <- fit_orig$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal]
    fit6$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin] <- fit_orig$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi]
  }
}

#relative
fit7=fit_orig
warmup <- args$sampling_args$iter/2
for (chain in 1:nchains) {
  order_cam <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_cam_R0]][-(1:warmup)])
  order_mal<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_mal_R0]][-(1:warmup)])
  order_phi <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_phi_R0]][-(1:warmup)])
  order_sin<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_sin_R0]][-(1:warmup)])
  order_tha <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_tha_R0]][-(1:warmup)])
  for (i in 1:length(idx_cam)) {
    fit7$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam] <- fit_orig$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi]
    fit7$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi] <- fit_orig$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha]
    fit7$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha] <- fit_orig$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal]
    fit7$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal] <- fit_orig$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin]
    fit7$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin] <- fit_orig$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam]
  }
}

#relative
fit8=fit_orig
warmup <- args$sampling_args$iter/2
for (chain in 1:nchains) {
  order_cam <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_cam_R0]][-(1:warmup)])
  order_mal<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_mal_R0]][-(1:warmup)])
  order_phi <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_phi_R0]][-(1:warmup)])
  order_sin<- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_sin_R0]][-(1:warmup)])
  order_tha <- warmup + order(fit_orig$stanfit@sim$samples[[chain]][[idx_tha_R0]][-(1:warmup)])
  for (i in 1:length(idx_cam)) {
    fit8$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi] <- fit_orig$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam]
    fit8$stanfit@sim$samples[[chain]][[idx_cam[i]]][order_cam] <- fit_orig$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin]
    fit8$stanfit@sim$samples[[chain]][[idx_sin[i]]][order_sin] <- fit_orig$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal]
    fit8$stanfit@sim$samples[[chain]][[idx_mal[i]]][order_mal] <- fit_orig$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha]
    fit8$stanfit@sim$samples[[chain]][[idx_tha[i]]][order_tha] <- fit_orig$stanfit@sim$samples[[chain]][[idx_phi[i]]][order_phi]
  }
}

## plotting rt
# get posterior infections
seed = 100
rt <- posterior_rt(fit_orig,seed=seed)
rt_cf1 <- posterior_rt_(fit_orig, eta=e1$draws, seed=seed)
rt_cf2 <- posterior_rt_(fit_orig, eta=e2$draws, seed=seed)
rt_cf3 <- posterior_rt_(fit_orig, eta=e3$draws, seed=seed)
rt_cf4 <- posterior_rt_(fit_orig, eta=e4$draws, seed=seed)
rt_cf5 <- posterior_rt(fit5, seed=seed)
rt_cf6 <- posterior_rt(fit6, seed=seed)
rt_cf7 <- posterior_rt(fit7, seed=seed)
rt_cf8 <- posterior_rt(fit8, seed=seed)

deaths_obs <-args$data$deaths
deaths_fit <- posterior_predict(fit_orig, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf1 <- posterior_predict_(fit_orig, eta=e1$draws, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf2 <- posterior_predict_(fit_orig, eta=e2$draws, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf3 <- posterior_predict_(fit_orig, eta=e3$draws, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf4 <- posterior_predict_(fit_orig, eta=e4$draws, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf5 <- posterior_predict(fit5, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf6 <- posterior_predict(fit6, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf7 <- posterior_predict(fit7, type="deaths",posterior_mean=TRUE, seed=seed)
deaths_cf8 <- posterior_predict(fit8, type="deaths",posterior_mean=TRUE, seed=seed)

infections <- posterior_infections(fit_orig,poster_mean=TRUE, seed=seed)
infections_cf1 <- posterior_infections_(fit_orig, eta=e1$draws, posterior_mean=TRUE, seed=seed)
infections_cf2 <- posterior_infections_(fit_orig, eta=e2$draws, posterior_mean=TRUE, seed=seed)
infections_cf3 <- posterior_infections_(fit_orig, eta=e3$draws, posterior_mean=TRUE, seed=seed)
infections_cf4 <- posterior_infections_(fit_orig, eta=e4$draws, posterior_mean=TRUE, seed=seed)
infections_cf5 <- posterior_infections(fit5, posterior_mean=TRUE, seed=seed)
infections_cf6 <- posterior_infections(fit6,posterior_mean=TRUE, seed=seed)
infections_cf7 <- posterior_infections(fit7, posterior_mean=TRUE, seed=seed)
infections_cf8 <- posterior_infections(fit8,posterior_mean=TRUE, seed=seed)
#
rt$group <- as.character(rt$group)
rt_cf1$group <- as.character(rt_cf1$group)
rt_cf2$group <- as.character(rt_cf2$group)
rt_cf3$group <- as.character(rt_cf3$group)
rt_cf4$group <- as.character(rt_cf4$group)
rt_cf5$group <- as.character(rt_cf5$group)
rt_cf6$group <- as.character(rt_cf6$group)
rt_cf7$group <- as.character(rt_cf7$group)
rt_cf8$group <- as.character(rt_cf8$group)

df<- data.frame(
  date = rt$time,
  median = apply(rt$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_fit$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_fit$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_fit$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections$draws, 2, function(x) quantile(x, 0.975)),
  group = rt$group
)
df_cf1 <- data.frame(
  date = rt_cf1$time,
  median = apply(rt_cf1$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf1$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf1$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf1$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf1$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf1$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf1$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf1$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf1$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf1$group
)
df_cf2 <- data.frame(
  date = rt_cf2$time,
  median = apply(rt_cf2$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf2$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf2$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf2$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf2$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf2$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf2$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf2$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf2$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf2$group
)
df_cf3 <- data.frame(
  date = rt_cf3$time,
  median = apply(rt_cf3$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf3$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf3$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf3$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf3$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf3$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf3$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf3$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf3$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf3$group
)
df_cf4 <- data.frame(
  date = rt_cf4$time,
  median = apply(rt_cf4$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf4$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf4$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf4$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf4$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf4$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf4$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf4$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf4$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf4$group
)
df_cf5 <- data.frame(
  date = rt_cf5$time,
  median = apply(rt_cf5$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf5$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf5$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf5$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf5$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf5$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf5$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf5$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf5$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf5$group
)

df_cf6 <- data.frame(
  date = rt_cf6$time,
  median = apply(rt_cf6$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf6$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf6$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf6$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf6$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf6$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf6$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf6$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf6$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf6$group
)
df_cf7 <- data.frame(
  date = rt_cf7$time,
  median = apply(rt_cf7$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf7$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf7$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf7$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf7$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf7$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf7$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf7$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf7$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf7$group
)
df_cf8 <- data.frame(
  date = rt_cf8$time,
  median = apply(rt_cf8$draws, 2, function(x) quantile(x, 0.5)),
  median_li = apply(rt_cf8$draws, 2, function(x) quantile(x, 0.025)),
  median_ui = apply(rt_cf8$draws, 2, function(x) quantile(x, 0.975)),
  deaths_median = apply(deaths_cf8$draws, 2, function(x) quantile(x, 0.5)),
  deaths_li = apply(deaths_cf8$draws, 2, function(x) quantile(x, 0.025)),
  deaths_ui = apply(deaths_cf8$draws, 2, function(x) quantile(x, 0.975)),
  deaths=deaths_obs,
  infections_median = apply(infections_cf8$draws, 2, function(x) quantile(x, 0.5)),
  infections_li= apply(infections_cf8$draws, 2, function(x) quantile(x, 0.025)),
  infections_ui = apply(infections_cf8$draws, 2, function(x) quantile(x, 0.975)),
  group = rt_cf8$group
)

df$colour <- "Original"
df_cf1$colour <- "Scenario 1"
df_cf2$colour <- "Scenario 2"
df_cf3$colour <- "Scenario 3"
df_cf4$colour <- "Scenario 4"
df_cf5$colour <- "Scenario 5"
df_cf6$colour <- "Scenario 6"
df_cf7$colour <- "Scenario 7"
df_cf8$colour <- "Scenario 8"
df <- rbind(df, df_cf1, df_cf2, df_cf3, df_cf4,df_cf5, df_cf6, df_cf7, df_cf8)
df=df[df$date>=changeDate,]
#
saveRDS(df,"/Users///COVID-Counterfactual/Top5/2021_SI/counterfactual_2021_top5_iter200_SI.rds")
# making sure we have renamed everything as per scenarios in papers, which is
# abs/real: donor->recpient
#
source("/Users///COVID-Counterfactual/transform_df_plotting_top5.R")
#############################################################################
#plots
#
# colours thanks for Michael Betancourts aesthetics
ci <- c("#C79999")
mn <- c("#7C0000")
#
date_breaks <- "1 month"
base <- ggplot2::ggplot() +
  ggplot2::xlab("") +
  ggplot2::scale_x_date(
    date_breaks = date_breaks,
    labels = scales::date_format("%B")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 12)
  ) +
  ggplot2::theme(legend.position = "right")
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median, color = colour),
    data = df[df$country%in%c("Cambodia")& !df$group%in%c("Fit Cam"),],
    size = 1,color="red",
  ) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median),
    data = df_originals[df_originals$country%in%c("Cambodia"),],
    size = 1,color='black'
  ) +
  xlim(as.Date('2021-02-01'), as.Date('2022-01-01'))
p1 <- p1 + ggplot2::labs(y = "Median Rt") + ggplot2::facet_wrap(~group,scales="free_y") + scale_color_brewer(palette="Paired") + ggtitle("Effective R tectory for Cambodia (Recipient)")

ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-rt_cam_1.pdf",p1,width=9,height=6)

p2 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median, color = colour),
    data = df[df$country%in%c("Malaysia")& !df$group%in%c("Fit Mal"),],
    size = 1,color="red",
  ) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median),
    data = df_originals[df_originals$country%in%c("Malaysia") ,],
    size = 1,color='black') + xlim(as.Date('2021-02-01'), as.Date('2022-01-01'))
p2 <- p2 + ggplot2::labs(y = "Median Rt") + ggplot2::facet_wrap(~group,,scales="free_y") + scale_color_brewer(palette="Paired")  + ggtitle("Effective R tectory for Malaysia (Recipient)")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-rt_mal_1.pdf",p2,width=9,height=6)


p3 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median, color = colour),
    data = df[df$country%in%c("Philippines") & !df$group%in%c("Fit Phi"),],
    size = 1,color="red",
  ) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median),
    data = df_originals[df_originals$country%in%c("Philippines"),],
    size = 1,color='black'
  ) + xlim(as.Date('2021-02-01'), as.Date('2022-01-01'))
p3 <- p3 + ggplot2::labs(y = "Median Rt") + ggplot2::facet_wrap(~group,,scales="free_y") + scale_color_brewer(palette="Paired")  + ggtitle("Effective R tectory for Philippines (Recipient)")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-rt_phi_1.pdf",p3,width=9,height=6)


p4 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median, color = colour),
    data = df[df$country%in%c("Singapore") & !df$group%in%c("Fit Sin"),],
    size = 1,color="red",
  ) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median),
    data = df_originals[df_originals$country%in%c("Singapore"),],
    size = 1,color='black'
  ) + xlim(as.Date('2021-02-01'), as.Date('2022-01-01'))
p4 <- p4 + ggplot2::labs(y = "Median Rt") + ggplot2::facet_wrap(~group,,scales="free_y") + scale_color_brewer(palette="Paired")  + ggtitle("Effective R tectory for Singapore (Recipient)")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-rt_sin_1.pdf",p4,width=9,height=6)


p5 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median, color = colour),
    data = df[df$country%in%c("Thailand") & !df$group%in%c("Fit Tha"),],
    size = 1,color="red",
  ) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = median),
    data = df_originals[df_originals$country%in%c("Thailand"),],
    size = 1,color='black'
  ) + xlim(as.Date('2021-02-01'), as.Date('2022-01-01'))
p5 <- p5 + ggplot2::labs(y = "Median Rt") + ggplot2::facet_wrap(~group,,scales="free_y") + scale_color_brewer(palette="Paired")  + ggtitle("Effective R tectory for Thailand (Recipient)")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-rt_tha_1.pdf",p5,width=9,height=6)


margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
p1=p1+margin
p2=p2+margin
p3=p3+margin
p4=p4+margin
p5=p5+margin
library(cowplot)
g_rt <- plot_grid(p1, p2, p3, align = "v", nrow = 3, rel_heights = c(1/4, 1/4, 1/2))
g_rt <- grid.arrange(p1,p2,p3,p4,p5,ncol=1)  # default settings
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-rt_1.pdf",g_rt,width=11,height=11)
######################################################################################################################################################
# Deaths
######################################################################################################################################################
#remove df$group is NA
dfsub <- na.omit(df)
# removes all rows from dfsub where the value in the group column falls within the range of 6 to 25 
#(considering only the unique groups)
dfsub=dfsub[-which(dfsub$group %in% unique(dfsub$group)[c(6:25)]),]
#create new subgroup which has fit and rel
dfsub$newgroup <- rep(NA, nrow(dfsub))
dfsub$newgroup[which(dfsub$group %in% unique(dfsub$group)[c(1:5)])]="Fit"
dfsub$newgroup[-which(dfsub$group %in% unique(dfsub$group)[c(1:5)])]="Rel"
##
x=dfsub$median[which(dfsub$group%in%c("Fit Cam","Fit Phi","Fit Tha","Fit Mal","Fit Sin"))]
dfsub$rtfit <- rep(NA, nrow(dfsub))
dfsub$rtfit=rep(x,times=5)
dfrt1=data.frame("Group"=dfsub$group,"Value"=dfsub$median,"phi"=rep("Counterfactual",times=length(dfsub$median)),"Date"=dfsub$date)
dfrt1=dfrt1[-which(dfrt1$Group %in%c("Fit Philippines","Fit Cambodia","Fit Malaysia","Fit Singapore","Fit Thailand")),]
dfrt2=data.frame("Group"=dfsub$group,"Value"=dfsub$rtfit,"phi"=rep("Actual",times=length(dfsub$rtfit)),"Date"=dfsub$date)
#
dfrt=rbind.data.frame(dfrt1,dfrt2)
#
date_breaks <- "1 month"
#
p1 <-ggplot(dfrt,aes(x = Date, y = Value,color=phi)) + geom_line(size=1) +
  scale_color_manual(values = c("black","red")) +
  xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) +
  scale_x_date(
    date_breaks = date_breaks,
    labels = scales::date_format("%b")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 12)
  ) +
  ggplot2::theme(legend.position = "bottom")

p1$labels$colour="Line Type"
p1 <- p1 + ggplot2::labs(y = "Median Rt") + ggplot2::facet_wrap(~Group,scales="free_y") + 
  ggtitle("Effective R tectory March 15- Dec 31 2020") + theme(legend.title=element_text(size=12))+
  theme(legend.text=element_text(size=12)) + xlab("")

ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/next5rt2020two.pdf",p1,width=13,height=8)

######################################################################################################################################################
# Deaths
######################################################################################################################################################
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = deaths_median, color = colour),
    data = df,
    size = 1,color=mn,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=deaths), stat="identity",
    data = df,
    width = 1.5,fill='skyblue4',alpha=0.7,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = deaths_li,ymax=deaths_ui),
    data = df,
    size = 1,fill=ci,alpha=0.5,
  ) +
  xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) +  scale_x_date(
    date_breaks = date_breaks,
    labels = scales::date_format("%b")
  )
p1 <- p1 + ggplot2::labs(y = "Deaths") + ggplot2::facet_wrap(~group,scales="free_y") + 
  scale_color_brewer(palette="Paired") + ggtitle("Deaths March 15- Dec 31 2020")+xlab("")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/death2020two.pdf",p1,width=13,height=8)
#########################
# country deaths
#########################
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = deaths_median, color = colour),
    data = df[df$country%in%c("Malaysia"),],
    size = 1,color=mn,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = deaths_li,ymax=deaths_ui),stat="identity",
    data = df[df$country%in%c("Malaysia"),],
    size = 1,fill=ci,alpha=0.5,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=deaths,group = 1), stat="identity",
    data = df[df$country%in%c("Malaysia"),],
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  xlim(as.Date('2021-03-01'), as.Date('2022-01-01')) 
p1 <- p1 + ggplot2::labs(y = "Deaths") + ggplot2::facet_wrap(~group, scales="free_y") + scale_color_brewer(palette="Paired") + ggtitle("Malaysia")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-deaths_mal.pdf",p1,width=9,height=6)
#
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = deaths_median, color = colour),
    data = df[df$country%in%c("Cambodia"),],
    size = 1,color=mn,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = deaths_li,ymax=deaths_ui),stat="identity",
    data = df[df$country%in%c("Cambodia"),],
    size = 1,fill=ci,alpha=0.5,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=deaths,group = 1), stat="identity",
    data = df[df$country%in%c("Cambodia"),],
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  xlim(as.Date('2021-03-01'), as.Date('2022-01-01')) 
p1 <- p1 + ggplot2::labs(y = "Deaths") + ggplot2::facet_wrap(~group, scales="free_y") + scale_color_brewer(palette="Paired") + ggtitle("Cambodia")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-deaths_cam.pdf",p1,width=9,height=6)
#
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = deaths_median, color = colour),
    data = df[df$country%in%c("Philippines"),],
    size = 1,color=mn,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = deaths_li,ymax=deaths_ui),stat="identity",
    data = df[df$country%in%c("Philippines"),],
    size = 1,fill=ci,alpha=0.5,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=deaths,group = 1), stat="identity",
    data = df[df$country%in%c("Philippines"),],
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  xlim(as.Date('2021-03-01'), as.Date('2022-01-01')) 
p1 <- p1 + ggplot2::labs(y = "Deaths") + ggplot2::facet_wrap(~group, scales="free_y") + scale_color_brewer(palette="Paired") + ggtitle("Philippines")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-deaths_phi.pdf",p1,width=9,height=6)
#
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = deaths_median, color = colour),
    data = df[df$country%in%c("Singapore"),],
    size = 1,color=mn,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = deaths_li,ymax=deaths_ui),stat="identity",
    data = df[df$country%in%c("Singapore"),],
    size = 1,fill=ci,alpha=0.5,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=deaths,group = 1), stat="identity",
    data = df[df$country%in%c("Singapore"),],
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  xlim(as.Date('2021-03-01'), as.Date('2022-01-01')) 
p1 <- p1 + ggplot2::labs(y = "Deaths") + ggplot2::facet_wrap(~group, scales="free_y") + scale_color_brewer(palette="Paired") + ggtitle("Singapore")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-deaths_sin.pdf",p1,width=9,height=6)
#
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = deaths_median, color = colour),
    data = df[df$country%in%c("Thailand"),],
    size = 1,color=mn,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = deaths_li,ymax=deaths_ui),stat="identity",
    data = df[df$country%in%c("Thailand"),],
    size = 1,fill=ci,alpha=0.5,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=deaths,group = 1), stat="identity",
    data = df[df$country%in%c("Thailand"),],
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  xlim(as.Date('2021-03-01'), as.Date('2022-01-01')) 
p1 <- p1 + ggplot2::labs(y = "Deaths") + ggplot2::facet_wrap(~group, scales="free_y") + scale_color_brewer(palette="Paired") + ggtitle("Thailand")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-deaths_tha.pdf",p1,width=9,height=6)
#
#######################################################################################################################################################
# Infections
#######################################################################################################################################################
p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = infections_median, color = colour),
    data = df[df$country%in%c("Cambodia"),],
    size = 1,color=mn,
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = infections_li,ymax=infections_ui),
    data = df[df$country%in%c("Cambodia"),],
    size = 1,fill=ci,alpha=0.5,
  )  +
  xlim(as.Date('2021-02-01'), as.Date('2021-02-01'))
p1 <- p1 + 
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = Date, y= Cambodia), stat="identity",
    data = data_cases,
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  ggplot2::labs(y = "Infections") + ggplot2::facet_wrap(~group) + scale_color_brewer(palette="Paired") + ggtitle("Cambodia")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-infections_cam.pdf",p1,width=9,height=6)
#
p2 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = infections_median, color = colour),
    data = df[df$country%in%c("Malaysia"),],
    size = 1,color=mn,
  )+ xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = infections_li,ymax=infections_ui),
    data = df[df$country%in%c("Malaysia"),],
    size = 11,fill=ci,alpha=0.5,
  )
p2 <- p2 + 
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = Date, y= Malaysia), stat="identity",
    data = data_cases,
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  ggplot2::labs(y = "Infections") + ggplot2::facet_wrap(~group) + scale_color_brewer(palette="Paired") + ggtitle("Malaysia")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-infections_Mal.pdf",p2,width=9,height=6)
#
p3 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = infections_median, color = colour),
    data = df[df$country%in%c("Philippines"),],
    size = 1,color=mn,
  )+ xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) 
# ggplot2::geom_ribbon(
#   mapping = ggplot2::aes(x = date, ymin = infections_li, ymax = infections_ui),
#   data = df[df$country%in%c("Philippines"),],
#   size = 1,fill=ci,alpha=0.5,
# )
p3 <- p3 + 
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = Date, y= Philippines), stat="identity",
    data = data_cases,
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  ggplot2::labs(y = "Infections") + ggplot2::facet_wrap(~group) + scale_color_brewer(palette="Paired") + ggtitle("Phillipines")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-infections_Phil.pdf",p3,width=9,height=6)
#
p4 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = infections_median, color = colour),
    data = df[df$country%in%c("Singapore"),],
    size = 1,color=mn,
  )+ xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = infections_li,ymax=infections_ui),
    data = df[df$country%in%c("Singapore")&df$group%in%c("Fit Sin"),],
    size = 1,fill=ci,alpha=0.5,
  )
p4 <- p4 + 
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = Date, y= Singapore), stat="identity",
    data = data_cases,
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  ggplot2::labs(y = "Infections") + ggplot2::facet_wrap(~group) + scale_color_brewer(palette="Paired") + ggtitle("Singapore")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-infections_sin.pdf",p4,width=9,height=6)
#
p5 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = infections_median, color = colour),
    data = df[df$country%in%c("Thailand"),],
    size = 1,color=mn,
  )+ xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = infections_li,ymax=infections_ui),
    data = df[df$country%in%c("Thailand")&df$group%in%c("Fit Tha"),],
    size = 1,fill=ci,alpha=0.5,
  )
p5 <- p5 + 
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = Date, y= Thailand), stat="identity",
    data = data_cases,
    width = 0.5,fill='skyblue4',alpha=0.8
  ) +
  ggplot2::labs(y = "Infections") + ggplot2::facet_wrap(~group) + scale_color_brewer(palette="Paired") + ggtitle("Thailand")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-infections_tha.pdf",p5,width=9,height=6)
#
margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
p1=p1+margin
p2=p2+margin
p3=p3+margin
p4=p4+margin
p5=p5+margin
g_infections <- grid.arrange(p1,p2,p3,ncol=1)  # default settings
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/scenario-infections_all.pdf",g_infections,width=11,height=11)

#######################################################
# dfsub ? 
#########

p1 <- base +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = date, y = infections_median, color = colour),
    data = dfsub,
    size = 1,color=mn,
  ) +
  ggplot2::geom_bar(
    mapping = ggplot2::aes(x = date, y=infections), stat="identity",
    data = dfsub,
    width = 0.5,fill='skyblue4',alpha=0.7,
  ) +
  xlim(as.Date('2021-02-01'), as.Date('2021-12-31')) +  scale_x_date(
    date_breaks = date_breaks,
    labels = scales::date_format("%b")
  )
p1 <- p1 + ggplot2::labs(y = "Infections") + ggplot2::facet_wrap(~group,scales="free_y") + 
  scale_color_brewer(palette="Paired") + ggtitle("Infections March 15- Dec 31 2020")+xlab("")
ggsave("/Users///COVID-Counterfactual/Top5/2021_SI/CF_SI_figures/infections2020two.pdf",p1,width=13,height=8)
#
