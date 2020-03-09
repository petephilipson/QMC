# Run applications to PBC data
library(joineRML)
data("pbc2")
# Log transform for bilirubin
pbc2$log.b <- log(pbc2$serBilir)
pbc2$trans.pro <- (0.1 * pbc2$prothrombin)^-4

# BIVARIATE models
# Vanilla fit
fit.joineRML_biv <- mjoint(
  formLongFixed = list("log.bil" = log.b ~ year + age + drug, "alb" = albumin ~ year + age + drug),
  formLongRandom = list("log.bil" = ~ year | id, "alb" = ~ year | id),
  formSurv = Surv(years, status2) ~ age + drug,
  data = pbc2,
  timeVar = "year", control = list(type = "montecarlo", nMCmax = 12000, burnin = 5,
                                   tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5))
# Antithetic fit
fit.joineRML2_biv <- mjoint(
  formLongFixed = list("log.bil" = log.b ~ year + age + drug, "alb" = albumin ~ year + age + drug),
  formLongRandom = list("log.bil" = ~ year | id, "alb" = ~ year | id),
  formSurv = Surv(years, status2) ~ age + drug,
  data = pbc2,
  timeVar = "year", control = list(type = "antithetic", nMCmax = 12000, burnin = 5,
                                   tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5))
# Quasi MC
fit.joineRML3_biv <- mjoint(
  formLongFixed = list("log.bil" = log.b ~ year + age + drug, "alb" = albumin ~ year + age + drug),
  formLongRandom = list("log.bil" = ~ year | id, "alb" = ~ year | id),
  formSurv = Surv(years, status2) ~ age + drug,
  data = pbc2,
  timeVar = "year", control = list(type = "sobol", nMCmax = 12000, burnin = 5,
                                   tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5))

# Extract time
as.numeric(fit.joineRML_biv$comp.time[2], units = "secs")
as.numeric(fit.joineRML2_biv$comp.time[2], units = "secs")
as.numeric(fit.joineRML3_biv$comp.time[2], units = "secs")

# TRIVARIATE models
# Use same set-up as in BMC paper
# Covariate are log(bilirubin), albumin and transformed prothrombin
# Vanilla fit
fit.joineRML_tri <- mjoint(
  formLongFixed = list("log.bil" = log.b ~ year + age + drug, "alb" = albumin ~ year + age + drug , 
                       "trans.pro" = trans.pro ~ year + age + drug),
  formLongRandom = list("log.bil" = ~ year | id, "alb" = ~ year | id, "trans.pro" = ~ year | id),
  formSurv = Surv(years, status2) ~ age + drug, data = pbc2, timeVar = "year",
  control = list(type = "montecarlo", nMCmax = 12000, burnin = 5,
                 tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5))
# Antithetic fit
fit.joineRML2_tri <- mjoint(
  formLongFixed = list("log.bil" = log.b ~ year + age + drug, "alb" = albumin ~ year + age + drug , 
                       "trans.pro" = trans.pro ~ year + age + drug),
  formLongRandom = list("log.bil" = ~ year | id, "alb" = ~ year | id, "trans.pro" = ~ year | id),
  formSurv = Surv(years, status2) ~ age + drug, data = pbc2, timeVar = "year",
  control = list(type = "antithetic", nMCmax = 12000, burnin = 5,
                 tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5))
# Quasi MC
fit.joineRML3_tri <- mjoint(
  formLongFixed = list("log.bil" = log.b ~ year + age + drug, "alb" = albumin ~ year + age + drug , 
                       "trans.pro" = trans.pro ~ year + age + drug),
  formLongRandom = list("log.bil" = ~ year | id, "alb" = ~ year | id, "trans.pro" = ~ year | id),
  formSurv = Surv(years, status2) ~ age + drug, data = pbc2, timeVar = "year",
  control = list(type = "sobol", nMCmax = 12000, burnin = 5,
                 tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5))

# Extract SEs from a generic fitted object
seFun <- function(u) {
  sqrt(diag((vcov(u))))
}



