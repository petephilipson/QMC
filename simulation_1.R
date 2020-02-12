# Simulation study for bivariate model
library(joineRML)
library(reshape2)
library(ggplot2)

# Specification of 'true' parameters

beta <- rbind(c(0, 1, 1, 1),
              c(0, -1, 0, 0.5))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5
all(eigen(D)$values > 0) # check covarinace matrix is PD

sigma2 <- c(0.25, 0.25)
gamma.x <- c(0, 1)
gamma.y <- c(-0.5, 1)

ltri <- lower.tri(D, diag = TRUE)
theta <- c(D[lower.tri(D, diag = TRUE)], beta[1, ], beta[2, ],
           sigma2, gamma.x, gamma.y)

#------------------------------------------------------------------------

# Fit a joint model to a single simulated dataset given the number of patients

simFit <- function(n, ...) {

  # Simulate a dataset
  sim <- simData(n = n, ntms = 6, beta = beta, D = D, sigma2 = sigma2,
                 gamma.x = gamma.x, gamma.y = gamma.y,
                 theta0 = -3.5, theta1 = 0.25, censlam = 0.05)

  #plot(survfit(Surv(survtime, cens) ~ 1, data = sim$survdat))

  # Fit the joint model (vanilla MC)
  fit_sim1 <- mjoint(
    formLongFixed = list("Y1" = Y.1 ~ time + ctsxl + binxl,
                         "Y2" = Y.2 ~ time + ctsxl + binxl),
    formLongRandom = list("Y1" = ~ time | id,
                          "Y2" = ~ time | id),
    formSurv = Surv(survtime, cens) ~ ctsx + binx,
    data = sim$longdat,
    survData = sim$survdat,
    timeVar = "time",
    control = list(type = "montecarlo", nMCmax = 12000, burnin = 5,
                   tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5),
    verbose = FALSE,
    ...)

  # Fit the joint model (antithetic MC)
  fit_sim2 <- mjoint(
    formLongFixed = list("Y1" = Y.1 ~ time + ctsxl + binxl,
                         "Y2" = Y.2 ~ time + ctsxl + binxl),
    formLongRandom = list("Y1" = ~ time | id,
                          "Y2" = ~ time | id),
    formSurv = Surv(survtime, cens) ~ ctsx + binx,
    data = sim$longdat,
    survData = sim$survdat,
    timeVar = "time",
    control = list(type = "antithetic", nMCmax = 12000, burnin = 5,
                   tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5),
    verbose = FALSE,
    ...)

  # Fit the joint model (QMC with Sobol)
  fit_sim3 <- mjoint(
    formLongFixed = list("Y1" = Y.1 ~ time + ctsxl + binxl,
                         "Y2" = Y.2 ~ time + ctsxl + binxl),
    formLongRandom = list("Y1" = ~ time | id,
                          "Y2" = ~ time | id),
    formSurv = Surv(survtime, cens) ~ ctsx + binx,
    data = sim$longdat,
    survData = sim$survdat,
    timeVar = "time",
    control = list(type = "sobol", nMCmax = 12000, burnin = 5,
                   tol0 = 1e-03, tol2 = 5e-03, nMCscale = 5),
    verbose = FALSE,
    ...)

  return(list("vanilla" = fit_sim1,
              "antithetic" = fit_sim2,
              "qmc" = fit_sim3))

}

#------------------------------------------------------------------------

# Run the simulation (serial mode)

N <- 100 # number of simulations
n <- 250 # number of subjects

ests <- array(dim = c(length(theta), N, 3))
ses <- array(dim = c(length(theta), N, 3))
props <- array(dim = c(7, N, 3))

set.seed(12345)

for (i in 1:N) {

  # Run simulation
  print(paste("Iteration:", i))
  x <- simFit(n = n)

  # Extract estimates
  parFun <- function(u) {
    with(u$coefficients, c(D[lower.tri(D, diag = TRUE)], beta, sigma2, gamma))
  }
  ests[, i, 1] <- parFun(x[[1]])
  ests[, i, 2] <- parFun(x[[2]])
  ests[, i, 3] <- parFun(x[[3]])

  # Extract SEs
  seFun <- function(u) {
    sqrt(diag((vcov(u))))
  }
  ses[, i, 1] <- seFun(x[[1]])
  ses[, i, 2] <- seFun(x[[2]])
  ses[, i, 3] <- seFun(x[[3]])

  # Extract properties
  propFun <- function(u) {
    c(as.numeric(u$comp.time[2], units = "secs"),
      ncol(u$history),
      u$finalnMC,
      u$conv,
      u$sfit$nevent / u$sfit$n, # avg. event rate
      u$dims$nk[1] / u$dims$n, # avg. number of measurements
      mean((u$sfit$y[, 1] < 6) & (u$sfit$y[, 2] == 0)) # avg. censoring before end of follow-up
      )
  }
  props[, i, 1] <- propFun(x[[1]])
  props[, i, 2] <- propFun(x[[2]])
  props[, i, 3] <- propFun(x[[3]])

  dimnames(ses)[[1]] <- names(diag(vcov(x[[1]])))
  dimnames(ests)[[1]] <- dimnames(ses)[[1]]
  dimnames(props)[[1]] <- c("time", "iter", "nMC", "conv", "event_rate", "avg_meas", "cens")

  dimnames(ses)[[3]] <- c("vanilla", "antithetic", "qmc")
  dimnames(ests)[[3]] <- dimnames(ses)[[3]]
  dimnames(props)[[3]] <- dimnames(ses)[[3]]
  
  gc()

}

#------------------------------------------------------------------------

# Bias plot

ests_cen <- ests - theta
est_cen_melt <- melt(ests_cen)
colnames(est_cen_melt) <- c("var", "iter", "method", "diff")

ggplot(aes(y = diff, x = var, fill = method), data = est_cen_melt) +
  geom_boxplot() +
  labs(x = "", y = "bias")

# Time plot

props_melt <- melt(props)
colnames(props_melt) <- c("prop", "iter", "method", "value")

ggplot(aes(y = value, x = method), 
       data = subset(props_melt, prop == "time")) +
  geom_boxplot() +
  labs(x = "", y = "convergence time (seconds)") +
  scale_y_log10()

