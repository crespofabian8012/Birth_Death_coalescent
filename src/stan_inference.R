library(tibble)
library(tidyr)
library(dplyr)

getCurrentFileLocation <-  function()
{
  require(tibble)
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value,
      into = c("key", "value"),
      sep = "=",
      fill = 'right'
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
path = getCurrentFileLocation()
source(paste0(path,"/CoalSimulationBirthDeath.R"))

list_packages <-
  c(
    "parallel",
    "doFuture",
    "bigstatsr",
    "doRNG",
    "microbenchmark",
    "matrixStats",
    "cmna",
    "NLRoot",
    "ggplot2",
    "future.apply",
    "ape",
    "tibble",
    "ggtree",
    "geiger",
    "dplyr",
    "pracma",
    "foreach",
    "doParallel",
    "stats",
    "doSNOW",
    "tcltk",
    "ggplot2",
    "tidyverse",
    "pryr",
    "lineup"
  )


install_required_packages(list_packages)
lapply(list_packages,
       library,
       character.only = TRUE)
#########################################################################
#### Simulate data from BD coalescent and deterministic 
#### M_K approximation model 
#########################################################################
#########################################################################
sim = 100
sample.size = 20
K = 2
path_to_save = getCurrentFileLocation()

GammaList = c(100, 10)
DeltaList = rep(1, length(GammaList))
#DeltaList = GammaList

Time.Origin.STD <-
  bigstatsr::FBM(length(DeltaList), sim , type = "double", init = 0)
coal.events.times.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))
number.ancestors.simA = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                   1))
number.ancestors.Transition = bigstatsr::FBM(length(DeltaList), sim)

listDataFramesA <-
  simulateA.parallel(
    DeltaList,
    GammaList,
    sim,
    sample.size,
    Time.Origin.STD,
    coal.events.times.simA,
    number.ancestors.simA,
    number.ancestors.Transition
  )
saveRDS(
  listDataFramesA,
  file = paste(
    path_to_save,
    "/",
    "listDataFramesA_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

saveRDS(
  Time.Origin.STD,
  file = paste(
    path_to_save,
    "/",
    "Time.Origin.STD_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

meansTorigin <-
  bigstatsr::big_apply(
    Time.Origin.STD,
    a.FUN = function(Time.Origin.STD, ind)
      rowMeans(Time.Origin.STD[ind,]),
    ind = rows_along(Time.Origin.STD),
    a.combine = 'c',
    block.size = 500
  )

saveRDS(
  meansTorigin,
  file = paste(
    path_to_save,
    "/",
    "meansToriginA_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)

listToriginCI <-
  bigstatsr::big_apply(
    Time.Origin.STD,
    a.FUN = function(Time.Origin.STD, ind) {
      quants <- c(0.025, 0.50, 0.975)
      matrixStats::rowQuantiles(Time.Origin.STD[ind,],  probs = quants)
    } ,
    ind = rows_along(Time.Origin.STD),
    a.combine = "rbind",
    block.size = 500
  )
saveRDS(
  listToriginCI,
  file = paste(
    path_to_save,
    "/",
    "listToriginCI_A_",
    sample.size,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)
coal.events.times.simB = bigstatsr::FBM(length(DeltaList), sim * (sample.size -
                                                                    1))                                                   
listDataFramesB <-
  simulateB_K.parallel(
    DeltaList,
    GammaList,
    sim,
    sample.size,
    Time.Origin.STD,
    coal.events.times.simB,
    K
  )
saveRDS(
  listDataFramesB,
  file = paste(
    path_to_save,
    "/",
    "listDataFramesB_",
    sample.size,
    "_K=",
    K,
    "_",
    sim,
    ".rds",
    sep = ""
  )
)
#########################################################################
#### Bayesian inference with stan
#########################################################################
#########################################################################
meanCoalescenTimes =listDataFramesB[[1]]$mean
medianCoalescenTimes =listDataFramesB[[1]]$median
meanTimeOrigin = meansTorigin[1]
medianTimeOrigin = listToriginCI[3]
trueDelta =100
library(extraDistr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(rstan)
library(bayesplot)
library(loo)
library("bayesplot")
library(ggmcmc)
library(ggplot2)
library(coda)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
path_to_stan_model_1population = paste0(path_to_save, "/stan_models/stan_model_one_population.stan")

#########################################################################
#### checking stan functions
#########################################################################
#########################################################################
expose_stan_functions(path_to_stan_model_1population)
conditionalDensityTOrigin_pdf <- function(y, delta) {
  exp(sapply(y, FUN = conditionalDensityTOrigin_lpdf, delta=delta ))
}
print(integrate(conditionalDensityTOrigin_pdf, lower = 0, upper = Inf, trueDelta))
yr <- conditionalDensityTOrigin_rng(delta=trueDelta)
all.equal(integrate(conditionalDensityTOrigin_pdf, lower = 0, upper = yr,  trueDelta)$value,
          conditionalDensityTOrigin_cdf(yr, delta=trueDelta))

all.equal(exp(conditionalDensityTOrigin_lcdf(yr, delta=trueDelta)),conditionalDensityTOrigin_cdf(yr,delta=trueDelta))
#########################################################################
#### checking ended
#########################################################################
#########################################################################
#########################################################################
#### Run stan model
#########################################################################
#########################################################################
K = 2
theta =1
scaled.meanCoalescenTimes = theta *meanCoalescenTimes
input_stan<-list(K=K, sample_size=sample.size, sorted_coalescent_times_scaled_by_theta= c(0.0,scaled.meanCoalescenTimes))
fit_stan <- stan(file=path_to_stan_model_1population, data=input_stan, refresh=1,
                chains=4, seed=803214053)
#########################################################################
#### Post processing the stan output
#########################################################################
#########################################################################
check_treedepth(fit_stan)
check_energy(fit_stan)
check_div(fit_stan)
summary(fit_stan)
rstan::plot(fit_stan)
bayesplot::mcmc_trace(as.matrix(fit_stan),pars=c("delta","torigin"),
                facet_args = list(nrow = 2))
bayesplot::mcmc_areas(as.matrix(fit_stan),pars=c("delta","torigin"),prob = 0.5) +
  ggtitle("Posterior distributions with medians and 50% intervals")

rstan::traceplot(fit_stan)
ainfo <- rstan::get_adaptation_info(fit_stan)
cat(ainfo[[1]])
seed <- rstan::get_seed(fit_stan)
sp <- rstan::get_sampler_params(fit_stan)
sp2 <- rstan::get_sampler_params(fit_stan, inc_warmup = FALSE)
head(sp[[1]])

lp <- rstan::log_prob(fit_stan, c(1, 2))
grad <- rstan::grad_log_prob(fit_stan, c(1, 2))
lp2 <- rstan::attr(grad, "log_prob") # should be the same as "lp"
lp_cp <- bayesplot::log_posterior(fit_stan)
head(lp_cp)
np_cp <- bayesplot::nuts_params(fit_stan)
head(np_cp)
posterior_cp <- as.array(fit_stan)
color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp)
bayesplot::mcmc_pairs(posterior_cp, np = np_cp, pars = c("delta","torigin"),
           off_diag_args = list(size = 0.75))
scatter_theta_cp <- bayesplot::mcmc_scatter(
  posterior_cp,
  pars = c("delta", "torigin"),
 # transform = list(delta = "log", torigin= "log"), # can abbrev. 'transformations'
  np = np_cp,
  size = 1
)
scatter_theta_cp
color_scheme_set("mix-brightblue-gray")
bayesplot::mcmc_trace(posterior_cp, pars = "delta", np = np_cp) +
  xlab("Post-warmup iteration")
color_scheme_set("mix-brightblue-gray")
bayesplot::mcmc_trace(posterior_cp, pars = "torigin", np = np_cp) +
  xlab("Post-warmup iteration")
color_scheme_set("red")
bayesplot::mcmc_nuts_divergence(np_cp, lp_cp)
color_scheme_set("red")
bayesplot::mcmc_nuts_energy(np_cp)
rhats <- bayesplot::rhat(fit_stan)
print(rhats)
color_scheme_set("brightblue") # see help("color_scheme_set")
bayesplot::mcmc_rhat(rhats)+  yaxis_text(hjust = 1)
ratios_cp <- bayesplot::neff_ratio(fit_stan)
print(ratios_cp)
bayesplot::mcmc_neff(ratios_cp, size = 2)
bayesplot::mcmc_acf(posterior_cp, pars = c("delta", "torigin"), lags = 10)


fit_params <- rstan::extract(fit_stan)

mcmcList<-rstan::As.mcmc.list(fit_stan, pars=c("delta", "torigin"))
deltaPar <- ggs(mcmcList)
ggmcmc(deltaPar, plot = c("density", "running", "caterpillar"))
ggs_running(deltaPar)
ggs_autocorrelation(deltaPar)
ggs_crosscorrelation(deltaPar, absolute_scale = FALSE)
ggs_Rhat(deltaPar)
ggs_geweke(deltaPar)
ggs_density(deltaPar)

