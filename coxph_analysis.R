##### Purpose #####
#Cox regression analysis - both univariate and multivariate

#Note: To successfully use this without modifying the code, your datatable may have:
# Column 1: pin
# Column 2 to n: survival variables and events
# Column n+: covariates

#### Libraries ####
library(dplyr)
library(tidyverse)
library(lubridate)
library(survival)
library(writexl)

#### Set File Path ####
setwd("E:/file/path")

#### Data ####
## Please change these as needed!
data_full <- read.csv('datatable.csv')
# If you have multiple survival events that you would like to analyze separately, this is possible as seen below
survivaltime <- c('yeareventA_dead', 'yeareventA_censor', 'yeareventB_dead', 'yeareventB_censor')
survivalvars <- c('eventA', 'eventA', 'eventB', 'eventB')

# You will have to change "Cov1" to whatever covariate is left most in your dataframe
covariates <- colnames(data_full %>%
                         select(which(colnames(data_full)=="Cov1"):length(colnames(data_full))))


#### Functions ####
CoxPHCustom_Univ <- function(survivalvars, survivaltime, covariates, data_full) {
  univ_surv <- c()
  for (i in 1:length(survivalvars)) {
    for (j in 1:length(covariates)) {
      univ_surv <- c(univ_surv,
                     paste('Surv(', survivaltime[i], ', ', survivalvars[i], ', type = "right") ~', covariates[j],
                           collapse = ', ', sep = ''))
    }
  }
  univ_formulas <- sapply(univ_surv,
                          function(x) as.formula(paste(x)))
  univ_models <- lapply(univ_formulas,
                        function(x){coxph(x, data = data_full, ties = "breslow", tt = 365.25, id = pin)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           HR_p <- as.data.frame(x$coefficients)
                           HR_p <- subset(HR_p, select = -`exp(coef)`)
                           HR_p['covariate'] <- row.names(HR_p)
                           HR_p['n'] <- x$n
                           HR_p['n_event'] <- x$nevent
                           rownames(HR_p) <- 1:nrow(HR_p)
                           conf_limits <- as.data.frame(x$conf.int)
                           rownames(conf_limits) <- 1:nrow(conf_limits)
                         })
  univ_results2<-bind_rows(univ_results)
  univ_results2['event'] <- rep(survivalvars, each = length(univ_results2$covariate)/length(survivalvars))
  univ_results2 <- univ_results2[, c("event", "covariate", "n", "n_event",
                                     "coef", "se(coef)", "z", "Pr(>|z|)",     
                                     "exp(coef)", "exp(-coef)", "lower .95", "upper .95")]
  return(univ_results2)
}


## Set Variables to Factors.
data_full2 <- data_full %>%
  mutate_at(vars(which(colnames(data_full)=="age"):length(colnames(data_full))),
            as.factor)

## Univariate
univ_results1 <- CoxPHCustom_Univ(survivalvars, survivaltime, covariates, data_full)
univ_results2 <- univ_results1 %>%
  mutate(`p_value` = round(as.numeric(`Pr(>|z|)`),7),
         HR = round(`exp(coef)`, 7),
         lower_95 = round(`lower .95`, 7),
         upper_95 = round(`upper .95`, 7)) %>%
  left_join(covariate_labels,
            by = 'covariate') %>%
  select(event, covariate_label, HR, lower_95, upper_95, p_value)

## Multivariate
univ_results2_sig <- univ_results1 %>%
  filter(`Pr(>|z|)` <= 0.1) %>%
  mutate(multi_covariates = str_sub(covariate, 1, -2),
         time = if_else(str_detect(event, "eventA"),
                        "yeareventA_dead",
                        "yeareventA_censor")) %>%
  group_by(time, event) %>%
  distinct(multi_covariates) %>%
  mutate(multi_covariates_comb = str_c(multi_covariates, collapse = "+")) %>%
  distinct(multi_covariates_comb)

multi_surv <- mapply(function(x, y, z) paste('Surv(', y, ',', x, ')~', paste(z, collapse = "+"), sep = ""),
                     x = univ_results2_sig$event,
                     y = univ_results2_sig$time,
                     z = univ_results2_sig$multi_covariates_comb)

multi_formulas <- sapply(multi_surv,
                         function(x) as.formula(paste(x)))

multi_models <- lapply(multi_formulas,
                       function(x){coxph(x, data = data_full, ties = "breslow")})

multi_results <- lapply(multi_models,
                              function(x){ 
                                x <- summary(x)
                                HR_p <- as.data.frame(x$coefficients)
                                HR_p <- subset(HR_p, select = -`exp(coef)`)
                                HR_p['covariate'] <- row.names(HR_p)
                                HR_p['n'] <- x$n
                                HR_p['n_event'] <- x$nevent
                                rownames(HR_p) <- 1:nrow(HR_p)
                                conf_limits <- as.data.frame(x$conf.int)
                                rownames(conf_limits) <- 1:nrow(conf_limits)
                              })

multi_results2 <-bind_rows(multi_results, .id = "id")
multi_results3 <- multi_results2 %>%
  mutate(p_value = round(as.numeric(`Pr(>|z|)`),7),
         HR = round(`exp(coef)`, 7),
         lower_95 = round(`lower .95`, 7),
         upper_95 = round(`upper .95`, 7),
         event = id) %>%
  select(event, covariate, HR, lower_95, upper_95, p_value)








