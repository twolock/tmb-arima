library(TMB)
library(data.table)
library(ggplot2)
library(tmbstan)

options(mc.cores = parallel::detectCores())
# set.seed(1)

source('R/make_D_matrix.R')
source('R/utils.R')

compile('cpp/TMB_ARIMA_test.cpp', flags='-w')
dyn.load(dynlib('cpp/TMB_ARIMA_test'))

p <- 1
pac <- runif(p, -1, 1)
# pac <- 0.8
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

N <- 100
N_IS <- 100
N_forecast <- 0
IS_i <- sort(sample(1:(N-N_forecast), N_IS))
d <- 1
mod <- list(
 order = c(p, d, 0),
 ar = phi
)
mu <- arima.sim(
  model = mod,
  n = N-d
)
y_sd <- 1
y <- rnorm(N, mu, y_sd)[IS_i]
IS_map <- matrix(0, N_IS, N)
IS_map[cbind(1:N_IS, IS_i)] <- 1

in_dat <- list(
  y = y,
  D = make_D_matrix(N, d),
  constr_sd = 1,
  constr_type = 2,
  IS_map = IS_map
)
in_par <- with(in_dat, list(
  mu = rep(0, ncol(D)),
  log_sigma_mu = 0,
  logit_phi = rep(0.01, p),
  log_sigma_y = log(y_sd)
))
in_ran <- c('mu')
in_map <- c()
in_map <- make_TMB_map(in_par, 'log_sigma_y')

obj <- MakeADFun(
  data = in_dat,
  parameters = in_par,
  random = in_ran,
  method='L-BFGS-B',
  DLL = 'TMB_ARIMA_test',
  map = in_map
)
stan_fit <- tmbstan(obj, chains=4,
                    iter=2000, warmup=1500,
                    control=list(
                      adapt_delta=0.9
                    ))

mu_samples <- melt(data.table(i = 1:N, t(extract(stan_fit, par='mu')$mu)), id.vars='i')
ggplot() +
  geom_line(data=mu_samples, aes(x=i, y=value, group=variable), alpha=0.05) +
  geom_point(aes(x=IS_i, y=y), color='red', shape=16) +
  # geom_line(aes(x=1:N, y=obj$report()$mu), color='blue') +
  geom_line(aes(x=1:N, y=mu), color='red') +
  theme_bw() + coord_cartesian(expand = 0)

par_samples <- do.call(cbind, extract(stan_fit, par=setdiff(names(in_par), names(in_map))))
phi_samples <- do.call(cbind, lapply(apply(par_samples, 1, obj$report), function(x){x$phi}))
phi

opt <- do.call(optim, obj)
print(sd_obj <- sdreport(obj, getJointPrecision = T))

par_samples <- MASS::mvrnorm(1000, organize_means(sd_obj), solve(sd_obj$jointPrecision))
mu_samples <- melt(data.table(i = 1:N, t(par_samples[, colnames(par_samples) == 'mu'])), id.vars='i')
ggplot() +
  geom_line(data=mu_samples, aes(x=i, y=value, group=variable), alpha=0.05) +
  geom_point(aes(x=IS_i, y=y), color='red', shape=16) +
  geom_line(aes(x=1:N, y=obj$report()$mu), color='blue') +
  geom_line(aes(x=1:N, y=mu), color='red') +
  theme_bw() + coord_cartesian(expand = 0)

t(apply(phi_samples, 1, median))
summary(sd_obj, select='report')
phi
