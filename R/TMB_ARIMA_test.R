library(TMB)
library(data.table)
library(ggplot2)

# set.seed(10)
setwd('~/Documents/R Stuff/TMB ARIMA Test/')

source('~/Documents/GitHub/hiv-incidence-model/r/make_D_matrix.R')
source('~/Documents/GitHub/hiv-incidence-model/r/utils.R')

compile('TMB_ARIMA_Test.cpp', flags='-w')
dyn.load(dynlib('TMB_ARIMA_Test'))

p <- 3
pac <- c(0.5, runif(p-1, -1, 1))
stopifnot(p > 0)
phi <- work <- pac
if (p == 1) {
  phi <- pac
} else {
  for (j in 1L:(p - 1L)) {
    a <- phi[j + 1L]
    print(a)
    phi[1L:j] <- work[1L:j] <- work[1L:j] - a * phi[j:1L]
  }
}


N <- 200
d <- 1
mod <- list(
 order = c(p, d, 0),
 ar = phi
)
mu <- arima.sim(
  model = mod,
  n = N-d
)

y <- rnorm(N, mu, p * 10^(d-1))
in_dat <- list(
  y = y,
  D = make_D_matrix(N, d),
  constr_sd = 0.001,
  constr_type = 1
)
in_par <- with(in_dat, list(
  mu = rep(0, ncol(D)),
  log_sigma_mu = 0,
  logit_phi = rep(0, p),
  log_sigma_y = 0
))
in_ran <- c('mu')

obj <- MakeADFun(
  data = in_dat,
  parameters = in_par,
  random = in_ran,
  method='L-BFGS-B',
  DLL = 'TMB_ARIMA_Test'
)
# stop()
opt <- do.call(optim, obj)
print(sd_obj <- sdreport(obj, getJointPrecision = T))

par_samples <- MASS::mvrnorm(1000, organize_means(sd_obj), solve(sd_obj$jointPrecision))
mu_samples <- melt(data.table(i = 1:N, t(par_samples[, colnames(par_samples) == 'mu'])), id.vars='i')
ggplot() +
  geom_line(data=mu_samples, aes(x=i, y=value, group=variable), alpha=0.05) +
  geom_point(aes(x=1:N, y=y), color='red', shape=16) +
  geom_line(aes(x=1:N, y=obj$report()$mu), color='blue') +
  geom_line(aes(x=1:N, y=mu), color='red') +
  theme_bw() + coord_cartesian(expand = 0)

summary(sd_obj, select='report')
phi
