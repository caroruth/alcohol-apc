library(tidyverse)
library(janitor)
library(bamp)
library(reshape2)


acute_alc <- read_csv("/Users/carolinerutherford/OneDrive - cumc.columbia.edu/APC Methods - acute alcohol/CDC Wonder data/gordon alcohol deaths 020822 clean.csv") %>%
  clean_names() %>%
  mutate(crude_rate = 100000*deaths/population) %>%
  filter(single_year_age <= 74)

case_matrix <- acast(acute_alc, year ~ single_year_age, value.var = "deaths")
pop_matrix <- acast(acute_alc,year ~ single_year_age, value.var = "population")


# rw1 = random walk of first order (prior)
bapc <- bamp(case_matrix, pop_matrix,  age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 1, verbose = T, overdisp = T)

plot(bapc)

quantiles<-c(0.05,0.5,0.95)
age<-as.array(bapc$samples$age)
age<-apply(age,2,quantile,quantiles)
age <- tibble("ll" = age[1,], "coef" = age[2,], "ul" = age[3,], level = c(16:74)) %>%
  mutate(parameter = "age")

period<-as.array(bapc$samples$period)
period<-apply(period,2,quantile,quantiles)
period <- tibble("ll" = period[1,], "coef" = period[2,], "ul" = period[3,], level = c(1999:2020)) %>%
  mutate(parameter = "period")

cohort<-as.array(bapc$samples$cohort)
cohort<-apply(cohort,2,quantile,quantiles)
cohort <- tibble("ll" = cohort[1,], "coef" = cohort[2,], "ul" = cohort[3,], level = c(1925:2004)) %>%
  mutate(parameter = "cohort")

bapc_dat <- bind_rows(age, period, cohort) %>%
  mutate(parameter = factor(parameter, levels = c("age", "period", "cohort")))

ggplot(bapc_dat, aes(level, coef)) +
  geom_line() +
  geom_line(aes(y = ll), linetype = "dashed") +
  geom_line(aes(y = ul), linetype = "dashed") +
  facet_wrap(~parameter, scales = "free_x") +
  labs(y = "Effect", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
