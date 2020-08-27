#------ load packages ------
library(readxl)
setwd("../")
#------ load data ------
wuhan<- read_excel("Wuhan_parameter.xlsx", sheet = 1)
age_profile_new <- read_excel("Wuhan_parameter.xlsx", sheet = 2)
age_profile_new_clinical <- read_excel("Wuhan_parameter.xlsx", sheet = 3)
tele <- read_excel("Wuhan_parameter.xlsx", sheet = 4)
severity <- read_excel("Wuhan_parameter.xlsx", sheet = 5)
pcr <- read_excel("Wuhan_parameter.xlsx", sheet = 6)

#------ set par ------
set.seed(123)
runtime = 10000

# screen number
screen <- c(16568, 83, 552)

# wuhan data
dat_confirm <- wuhan$confirmed
dat_clinical <- wuhan$clinical
# sum(dat_confirm) + sum(dat_clinical)
# dat_clinical <- rep(0, 4)
mild <- severity$mild * 0.01
moderate <- severity$moderate * 0.01
moderate_cli <- c(0.50025018,0.683266932,0.757075472,0)
sc_cli <- c(0.49974982,0.316733068,0.242924528,0)

# pcr par
pcr_nume <- c(pcr$num[2], pcr$num[1], pcr$num[1], pcr$num[2])
pcr_denom <- c(pcr$denomin[2], pcr$denomin[1], pcr$denomin[1], pcr$denomin[2])

# survey par
survey_nume <- tele$med[1:4]
survey_denom <- tele$fever[1:4]

# -sensitivity analysis ---#
# lower limit
# survey_nume <- c(survey_denom[1] * 0.5357143, survey_denom[2] * 0.2111801, survey_denom[3] * 0.2474227, survey_denom[4] * 0.1891892)
# upper limit
# survey_nume <- c(survey_denom[1] * 0.7857143 , survey_denom[2] * 0.3478261, survey_denom[3] * 0.4329897, survey_denom[4] * 0.5135135)


# age profile
t(age_profile_new[, -1]) -> age_con
t(age_profile_new_clinical[, -1]) -> age_cli
age_s1 <- c(0.01077, 0.17873, 0.38067, 0.42981)
age_s2 <- c(0.03330, 0.19958, 0.36801, 0.39909)

#------ pcr rate ------
pcr_rate <- c()
for (i in 1:length(pcr_nume)) {
  rbinom(n = runtime, size = pcr_denom[i], prob = pcr_nume[i]/pcr_denom[i])/pcr_denom[i] -> tmp
  pcr_rate <- cbind(pcr_rate, tmp)
}

#------ pcr adjust ------
t(dat_confirm/t(pcr_rate)) -> dat_confirm_pcr

#------confirm mild -----
t(t(dat_confirm_pcr)*mild) -> dat_confirm_pcr_mild

#------confirm moderate ------
t(t(dat_confirm_pcr)*moderate) -> dat_confirm_pcr_moderate

#------confirm s&c ------
t(t(dat_confirm_pcr)*(1 - mild - moderate)) -> dat_confirm_pcr_sc

#----- confirmed age profile ------
# mild age
dat_confirm_pcr_mild_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_mild[i, ]*t(age_con)) -> dat_confirm_pcr_mild_age[[i]]
}

# moderate age
dat_confirm_pcr_moderate_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_moderate[i, ]*t(age_con)) -> dat_confirm_pcr_moderate_age[[i]]
}

# s&c age
dat_confirm_pcr_sc_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_sc[i, ]*t(age_con)) -> dat_confirm_pcr_sc_age[[i]]
}

#------clinical moderate -----
t(dat_clinical*t(moderate_cli)) -> dat_clinical_moderate

# x <- c(1,2)
# y <- matrix(1:12, nrow = 2)
# t(x * t(y))
# matrix(rep(1:2,6), nrow = 2, ncol = 6, byrow = T)

#------clinical sc ------
t(dat_clinical*t(sc_cli)) -> dat_clinical_sc

#------ clinical age profile ------
# moderate
t(dat_clinical_moderate[,1]*t(age_cli)) -> dat_clinical_moderate_age
# sc
t(dat_clinical_sc[,1]*t(age_cli)) -> dat_clinical_sc_age

#------ clinical age profile ------
t(dat_clinical*t(age_cli)) -> dat_clinical_age

#------ survey rate ------
s_nume <- c()
survey_rate <- c()
for (i in 1:length(survey_nume)) {
  rbinom(n = runtime, size = survey_denom[i], prob = survey_nume[i]/survey_denom[i]) -> tmp1
  s_nume <- cbind(s_nume, tmp1)
  tmp1 / survey_denom[i] -> tmp
  survey_rate <- cbind(survey_rate, tmp)
}
survey_rate_mean <- apply(s_nume, 1, sum) / sum(survey_denom)
# quantile(survey_rate[,1], 0.975)
# quantile(survey_rate[,1], 0.025)
# quantile(survey_rate[,2], 0.975)
# quantile(survey_rate[,2], 0.025)
# quantile(survey_rate[,3], 0.975)
# quantile(survey_rate[,3], 0.025)
# quantile(survey_rate[,4], 0.975)
# quantile(survey_rate[,4], 0.025)
#------ screen ------   # 0.474 is mild proportion between period 2-4???0.3 is moderate proportion between period 4-5
screen[1]*0.474*(1 - survey_rate_mean)/pcr_rate[, 1] -> s1
screen[2]*(1 - survey_rate_mean)/pcr_rate[, 2] -> s2

screen[1]*0.474*(1 - survey_rate_mean)/pcr_rate[, 1] +
  screen[1]*0.30*(1 - survey_rate_mean)/pcr_rate[, 1] -> s1_alt
screen[2]*(1 - survey_rate_mean)/pcr_rate[, 2] +
  screen[3]*(1 - survey_rate_mean)/pcr_rate[, 2] -> s2_alt

#------ medical consult ------
dat_confirm_pcr_mild + dat_confirm_pcr_moderate + dat_confirm_pcr_sc -> dat_confirm_pcr_all

dat_confirm_pcr_all_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_all[i, ]*t(age_con)) -> dat_confirm_pcr_all_age[[i]]
}

dat_fin <- c()
for (i in 1:runtime) {
  # i = 1
  tmp <- dat_confirm_pcr_all_age[[i]]
  rbind(dat_fin, apply(tmp, 1, sum)) -> dat_fin
}

t(apply(dat_clinical_age, 1, sum) + t(dat_fin)) -> dat_fin

screen_age <- c()
for (i in 1:runtime) {
  rbind(screen_age, s1[i]*age_s1 + s2[i]*age_s2) -> screen_age
}

dat_fin - screen_age -> dat_fin

r1 <- apply(dat_fin, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = T)

#------med sensitivity ------
dat_confirm_pcr_mild + dat_confirm_pcr_moderate + dat_confirm_pcr_sc -> dat_confirm_pcr_all

dat_confirm_pcr_all_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_all[i, ]*t(age_con)) -> dat_confirm_pcr_all_age[[i]]
}

dat_fin <- c()
for (i in 1:runtime) {
  # i = 1
  tmp <- dat_confirm_pcr_all_age[[i]]
  rbind(dat_fin, apply(tmp, 1, sum)) -> dat_fin
}

t(apply(dat_clinical_age, 1, sum) + t(dat_fin)) -> dat_fin


screen_age_alt <- c()
for (i in 1:runtime) {
  rbind(screen_age_alt, s1_alt[i]*age_s1 + s2_alt[i]*age_s2) -> screen_age_alt
}

dat_fin - screen_age_alt -> dat_fin

r1_alt <- apply(dat_fin, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = T)


#------ hospitalized ------
dat_confirm_pcr_moderate + dat_confirm_pcr_sc -> dat_confirm_pcr_all

dat_confirm_pcr_all_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_all[i, ]*t(age_con)) -> dat_confirm_pcr_all_age[[i]]
}

dat_fin <- c()
for (i in 1:runtime) {
  # i = 1
  tmp <- dat_confirm_pcr_all_age[[i]]
  rbind(dat_fin, apply(tmp, 1, sum)) -> dat_fin
}

t(apply(dat_clinical_age, 1, sum) + t(dat_fin)) -> dat_fin

r2 <- apply(dat_fin, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = T)

#------hos sensitivity ------
dat_confirm_pcr_moderate * 0.75 + dat_confirm_pcr_sc -> dat_confirm_pcr_all

dat_confirm_pcr_all_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_all[i, ]*t(age_con)) -> dat_confirm_pcr_all_age[[i]]
}

dat_fin <- c()
for (i in 1:runtime) {
  # i = 1
  tmp <- dat_confirm_pcr_all_age[[i]]
  rbind(dat_fin, apply(tmp, 1, sum)) -> dat_fin
}

t(apply(dat_clinical_moderate_age * 0.75, 1, sum) + apply(dat_clinical_sc_age, 1, sum) + t(dat_fin)) -> dat_fin

r2_alt <- apply(dat_fin, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = T)




#------ symptomatic ------
screen_age <- c()
for (i in 1:runtime) {
  rbind(screen_age, s1[i]*age_s1 + s2[i]*age_s2) -> screen_age
}

screen_age_alt <- c()
for (i in 1:runtime) {
  rbind(screen_age_alt, s1_alt[i]*age_s1 + s2_alt[i]*age_s2) -> screen_age_alt
}

dat_mild_pcr_age_fin <- c()
for (i in 1:runtime) {
  dat_confirm_pcr_mild_age[[i]] -> tmp
  rbind(dat_mild_pcr_age_fin, apply(tmp, 1, sum)) -> dat_mild_pcr_age_fin
}

dat_moderate_pcr_age_fin <- c()
for (i in 1:runtime) {
  dat_confirm_pcr_moderate_age[[i]] -> tmp
  rbind(dat_moderate_pcr_age_fin, apply(tmp, 1, sum)) -> dat_moderate_pcr_age_fin
}

dat_confirm_pcr_moderate + dat_confirm_pcr_sc -> dat_confirm_pcr_all

dat_confirm_pcr_all_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_all[i, ]*t(age_con)) -> dat_confirm_pcr_all_age[[i]]
}

dat_fin <- c()
for (i in 1:runtime) {
  # i = 1
  tmp <- dat_confirm_pcr_all_age[[i]]
  rbind(dat_fin, apply(tmp, 1, sum)) -> dat_fin
}

(dat_mild_pcr_age_fin - screen_age)/survey_rate + screen_age + dat_fin -> dat_fin

t(apply(dat_clinical_age, 1, sum) + t(dat_fin)) -> dat_fin

r3 <- apply(dat_fin, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = T)

# sym alternative sensitivity
dat_confirm_pcr_sc -> dat_confirm_pcr_all

dat_confirm_pcr_all_age <- list()
for (i in 1:runtime) {
  t(dat_confirm_pcr_all[i, ]*t(age_con)) -> dat_confirm_pcr_all_age[[i]]
}

dat_fin <- c()
for (i in 1:runtime) {
  # i = 1
  tmp <- dat_confirm_pcr_all_age[[i]]
  rbind(dat_fin, apply(tmp, 1, sum)) -> dat_fin
}

(dat_mild_pcr_age_fin + dat_moderate_pcr_age_fin - screen_age_alt)/survey_rate + screen_age_alt +  dat_fin -> dat_fin 

t(apply(dat_clinical_age, 1, sum) + t(dat_fin)) -> dat_fin

r3_alt <- apply(dat_fin, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = T)

#------ save result ------
r1 <- as.data.frame(r1)
r1_alt <- as.data.frame(r1_alt)
r2 <- as.data.frame(r2)
r2_alt <- as.data.frame(r2_alt)
r3 <- as.data.frame(r3)
r3_alt <- as.data.frame(r3_alt)

r1 <- round(r1, digits = 0)
r1$all <- apply(r1, 1, sum)
r1_alt <- round(r1_alt, digits = 0)
r1_alt$all <- apply(r1_alt, 1, sum)
r2 <- round(r2, digits = 0)
r2$all <- apply(r2, 1, sum)
r2_alt <- round(r2_alt, digits = 0)
r2_alt$all <- apply(r2_alt, 1, sum)
r3 <- round(r3, digits = 0)
r3$all <- apply(r3, 1, sum)
r3_alt <- round(r3_alt, digits = 0)
r3_alt$all <- apply(r3_alt, 1, sum)



r1$p <- "med"
r1_alt$p <- "med_alt"
r2$p <- "hos"
r2_alt$p <- "hos_alt"
r3$p <- "sym"
r3_alt$p <- "sym_alt"

rbind(r3, r1, r2) -> fin_result
rbind( r3_alt, r1_alt, r2_alt) -> fin_result_sensi
fin_result <- t(fin_result)
fin_result <- fin_result[, c("50%", "2.5%", "97.5%", "50%1", "2.5%1", "97.5%1", "50%2", "2.5%2", "97.5%2")]
fin_result_sensi <- t(fin_result_sensi)
fin_result_sensi <- fin_result_sensi[, c("50%", "2.5%", "97.5%", "50%1", "2.5%1", "97.5%1", "50%2", "2.5%2", "97.5%2")]
# write.csv(fin_result_sensi, "fin_result_sensi_new.csv", row.names = F)
# write.csv(fin_result, "fin_result.csv", row.names = F)
# write.csv(fin_result, "fin_result_noclinical_new.csv", row.names = F)


# we approximate our final result  to calculate the ratios (HFR, CFR ,MFR, CHR, MHR) 
# and rate (symptomatic rate, medical consultation rate, hospitalization rate)
# if the number smaller than 1000 will round to the nearest ten (e.g.872 to 870)
# if the number bigger than 1000 will round to the nearest hundred (e.g.1154 to 1200)

# The rate was calculated by dividing the population size in Wuhan 
#(0-19yrs: 2002195; 20-39yrs: 3885005; 40-59yrs: 3363579; >=60yrs: 1446191; total:10696970)



