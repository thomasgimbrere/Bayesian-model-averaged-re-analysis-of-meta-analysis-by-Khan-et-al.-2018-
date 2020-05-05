# Loading Packages and setting working directory

require(metaBMA)
require(tidyverse)
require(LaplacesDemon)
require(MASS)

setwd("C:/Users/Thomas/Documents/Master BDS/Masterthese/")

# Loading data
data <- read.csv("PFS HR Khan et al..csv", header = TRUE, sep = ";")

#--------------------------------------------------PRIORS---------------------------------------------------

# PRIORS FOR EFFECT SIZE

# Zero-centered prior, adjusted for medium HR
zero_HR <- prior("norm", param = c(0,log(1.47)))


# Informed prior, based on Avelumab trial.

# Calculate the logHR and SE from Avelumab trial
logHR <- log(0.90)
SE <- (log(1.08)-log(0.75))/(1.96*2)

# Compute prior
informed <- prior("norm", param = c(logHR,SE))



# PRIORS FOR HETEROGENEITY

# van Erp / inverse gamma
inverse_gamma <- prior("invgamma", c(shape = 1, scale = 0.15))


# Informed heterogeneity prior, for cancer research, based on overview by Rhodes et al.(2015)

# Generate samples from log t-distr with given precision using LaplacesDemon
set.seed(123)
distr <- rstp(n = 1e4, mu = -0.88, tau = (2.55)^2, nu = 5)

# Transform the generated samples from log(tau) ^2 to tau
e_distr <- exp(distr)
sqrt_e_distr <- sqrt(e_distr)

# Plot histogram
hist <- hist(sqrt_e_distr, probability = TRUE,
             xlim = c(0, 2), ylim = c(0, 3),
             breaks = 30)

# Fit a t-distribution on the tau estimates.
fit <- fitdistr(sqrt_e_distr, "t")
param <- unname(c(fit$estimate[1:2], round(fit$estimate[3])))

# Create a prior from the estimated t-distribution
informed_tau <- prior("t", param = param, lower = 0)

 
#----------------------------------------ANALYSES PER PRIOR COMBINATION-------------------------------------------------


# ZERO CENTERED X INFORMED TAU

# Fixed-effects
mf_zeroHR_informed <- meta_fixed(log..HR., SE, study, data, 
                                   d = zero_HR)

# Random-effects
mr_zeroHR_informed <- meta_random(log..HR., SE, study, data, 
                                    d = zero_HR,
                                    tau = informed_tau)

# Averaged effect
mb_zeroHR_informed <- meta_bma(log..HR., SE, study, data, 
                               d = zero_HR,
                               tau =informed_tau)


# ZERO CENTERED X INVERSE GAMMA

# Fixed-effects
mf_zeroHR_inverse <- meta_fixed(log..HR., SE, study, data, 
                                   d = zero_HR)

# Random-effects
mr_zeroHR_inverse <- meta_random(log..HR., SE, study, data, 
                                    d = zero_HR,
                                    tau = inverse_gamma)

# Averaged effect
mb_zeroHR_inverse <- meta_bma(log..HR., SE, study, data, 
                               d = zero_HR,
                               tau =inverse_gamma)

# AVELUMAB X INFORMED TAU

# Fixed-effects
mf_informed_informed <- meta_fixed(log..HR., SE, study, data, 
                                     d = informed)

# Random-effects
mr_informed_informed <- meta_random(log..HR., SE, study, data, 
                                      d = informed,
                                      tau = informed_tau)

# Averaged effect
mb_informed_informed <- meta_bma(log..HR., SE, study, data, 
                                 d = informed,
                                 tau = informed_tau)

# AVELUMAB X INVERSE GAMMA

# Fixed-effects
mf_informed_inverse<- meta_fixed(log..HR., SE, study, data, 
                                     d = informed)

# Random-effects
mr_informed_inverse <- meta_random(log..HR., SE, study, data, 
                                   d = informed,
                                   tau = inverse_gamma)

# Averaged effect
mb_informed_inverse <- meta_bma(log..HR., SE, study, data, 
                                d = informed,
                                tau = inverse_gamma)



#---------------------------------------------PLOTS------------------------------------------------------


# PRIORS EFFECT SIZE
par(mfrow=c(1,2))
plot(zero_HR)
plot(informed)

# PRIORS HETEROGENEITY
par(mfrow=c(1,2))
plot(inverse_gamma, from = 0, to = 2)
plot(informed_tau, from = 0, to = 2)


#POSTERIORS

par(mfrow=c(2,2))
# Zero-centered x informed
plot_posterior(mb_zeroHR_informed, main = " Zero-centered Effect Size x Rhodes et al. (2015)", sub = "median = -0.160 [-0.440, 0.165]")

# Zero-centered x inverse gamma
plot_posterior(mb_zeroHR_inverse, main = "Zero-centered Effect Size x Inverse Gamma", sub = "median = -0.166 [-0.347, -0.013]")

# Avelumab informed x informed
plot_posterior(mb_informed_informed, main = "Avelumab Informed Effect Size x Rhodes et al. (2015)", sub = "median = -0.131 [-0.273, 0.028]")

# Avelumab informed x inverse gamma
plot_posterior(mb_informed_inverse, main = "Avelumab Informed Effect Size x Inverse Gamma", sub = "median = -0.144 [-0.267, -0.019]")


#FOREST PLOTS
dev.off()
# Zero-centered x informed
plot_forest(mb_zeroHR_informed)

# Zero-centered x inverse gamma
plot_forest(mb_zeroHR_inverse)

# Avelumab informed x informed
plot_forest(mb_informed_informed)

# Avelumab informed x inverse gamma
plot_forest(mb_informed_inverse)


#HR VALUES, Averaged model (same is done for fixed-effect model and random-effects model, for Table 1.)

# Zero-centered x informed
MEDIAN_D <- exp(-0.152)
LOW_D <- exp(-0.476)
UP_D <- exp(-.176)

c(MEDIAN_D, LOW_D, UP_D)

# Zero-centered x inverse gamma
HR_C <- exp(-0.169)
LOW_C <- exp(-0.354)
UP_C <- exp(-.001)

c(HR_C, LOW_C, UP_C)

# Avelumab informed x informed

MEDIAN_B <- exp(-0.132)
LOW_B <- exp(-0.271)
UP_B <- exp(0.015)

c(MEDIAN_B, LOW_B, UP_B)

# Avelumab informed x inverse gamma

HR_A <- exp(-0.145)
LOW_A <- exp(-0.270)
UP_A <- exp(-0.021)

c(HR_A, LOW_A, UP_A)



#AVERAGED BAYES FACTORS 

# Zero-centered x informed
BF_zeroHR_informed <- mb_zeroHR_informed$inclusion$incl.BF

# Zero-centered x inverse gamma
BF_zeroHR_inverse <- mb_zeroHR_inverse$inclusion$incl.BF

# Avelumab informed x informed
BF_informed_informed <- mb_informed_informed$inclusion$incl.BF

# Avelumab informed x inverse gamma
BF_informed_inverse <- mb_informed_inverse$inclusion$incl.BF

table <- matrix(NA,4,2)
bayesfactors <- c("Zero-centered x informed","Zero-centered x inverse gamma","Avelumab informed x informed", "Avelumab informed x inverse gamma")
values <- c(BF_zeroHR_informed,BF_zeroHR_inverse,BF_informed_informed,BF_informed_inverse)
table[,1] <- bayesfactors
table[,2] <- values
table


