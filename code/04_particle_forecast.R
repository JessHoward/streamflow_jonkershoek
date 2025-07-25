# # Assume: cal_ddat loaded and the data wrangling has been done as per streamflow analysis (ie. data has been removed from before 2013-03-16)

# Extract posterior samples
post <- as.matrix(jags.out)

# Extract means for coefficients
mu0       <- mean(post[, "mu0"])
beta_decay <- mean(post[, "beta_decay"])
beta_rain  <- mean(post[, "beta_rain"])
beta_sin   <- mean(post[, "beta_season_sin"])
beta_cos   <- mean(post[, "beta_season_cos"])

# Process/obs SDs
sigma_proc <- 1 / sqrt(mean(post[, "tau_add"]))
sigma_obs  <- 1 / sqrt(mean(post[, "tau_obs"]))

# Covariates for particle filter
rain <- log(cal_ddat$rainfall_dayback + 1)   # or imputed rainfall if you prefer
doy <- as.numeric(format(cal_ddat$Date, "%j")) / 365
season_sin <- sin(2 * pi * doy)
season_cos <- cos(2 * pi * doy)

# Ensure all vectors are same length as y
stopifnot(length(y) == length(rain), length(y) == length(season_sin))

particle_filter_forecast_cov <- function(y, rain, season_sin, season_cos,
                                         mu0, beta_decay, beta_rain, beta_sin, beta_cos,
                                         sigma_proc, sigma_obs, N = 1000) {
  T <- length(y)
  particles <- matrix(NA, nrow = T, ncol = N)
  weights <- matrix(NA, nrow = T, ncol = N)
  estimates <- rep(NA, T)
  
  # Initialize particles using first observation
  particles[1, ] <- rnorm(N, mean = y[1], sd = sigma_obs)
  weights[1, ] <- dnorm(y[1], mean = particles[1, ], sd = sigma_obs)
  weights[1, ] <- weights[1, ] / sum(weights[1, ])
  estimates[1] <- sum(particles[1, ] * weights[1, ])
  
  for (t in 2:T) {
    # Resample
    resample_indices <- sample(1:N, size = N, replace = TRUE, prob = weights[t - 1, ])
    particles[t - 1, ] <- particles[t - 1, resample_indices]
    
    # Covariates at time t (handle NAs)
    rain_t <- ifelse(is.na(rain[t]), 0, rain[t])
    sin_t  <- ifelse(is.na(season_sin[t]), 0, season_sin[t])
    cos_t  <- ifelse(is.na(season_cos[t]), 0, season_cos[t])
    
    # Propagate
    mu_t <- mu0 + beta_decay * (particles[t - 1, ] - mu0) +
      beta_rain * rain_t + beta_sin * sin_t + beta_cos * cos_t
    
    # Check for NA in mu_t
    mu_t[is.na(mu_t)] <- mu0  # fallback if needed
    
    particles[t, ] <- rnorm(N, mean = mu_t, sd = sigma_proc)
    
    # Weighting
    if (!is.na(y[t])) {
      weights[t, ] <- dnorm(y[t], mean = particles[t, ], sd = sigma_obs)
      # Handle all-zero weights
      if (sum(weights[t, ], na.rm = TRUE) == 0 || any(is.na(weights[t, ]))) {
        weights[t, ] <- rep(1 / N, N)
      } else {
        weights[t, ] <- weights[t, ] / sum(weights[t, ])
      }
    } else {
      weights[t, ] <- rep(1 / N, N)
    }
    
    
    estimates[t] <- sum(particles[t, ] * weights[t, ])
  }
  
  return(list(particles = particles, weights = weights, estimates = estimates))
}

y <- cal_ddat$`Streamflow Ave`
dates <- cal_ddat$Date

pf <- particle_filter_forecast_cov(
  y = y,
  rain = rain,
  season_sin = season_sin,
  season_cos = season_cos,
  mu0 = mu0,
  beta_decay = beta_decay,
  beta_rain = beta_rain,
  beta_sin = beta_sin,
  beta_cos = beta_cos,
  sigma_proc = sigma_proc,
  sigma_obs = sigma_obs,
  N = 1000
)

ci <- apply(pf$particles, 1, quantile, probs = c(0.025, 0.5, 0.975))

plot(dates, y, type = "l", col = "gray", ylab = "Streamflow", xlab = "Date")
polygon(c(dates, rev(dates)), 
        c(ci[1, ], rev(ci[3, ])),
        col = adjustcolor("skyblue", alpha.f = 0.5), border = NA)
lines(dates, ci[2, ], col = "blue", lwd = 2)
abline(v = as.Date("2024-01-01"), col = "red", lty = 2)
legend("topright", legend = c("Observed", "PF Median", "95% CI", "Forecast Start"), 
       col = c("gray", "blue", "skyblue", "red"), lty = c(1, 1, NA, 2), pch = c(NA, NA, 15, NA))


####

particle_filter_forecast <- function(y, N = 1000, sigma_proc = sigma_proc, 
                                     sigma_obs = sigma_obs) {
  T <- length(y)
  particles <- matrix(NA, nrow = T, ncol = N)
  weights <- matrix(NA, nrow = T, ncol = N)
  estimates <- rep(NA, T)
  
  # Initialize particles using first observation
  particles[1, ] <- rnorm(N, mean = y[1], sd = sigma_obs)
  weights[1, ] <- dnorm(y[1], mean = particles[1, ], sd = sigma_obs)
  weights[1, ] <- weights[1, ] / sum(weights[1, ])
  estimates[1] <- sum(particles[1, ] * weights[1, ])
  
  for (t in 2:T) {
    # Resample particles
    resample_indices <- sample(1:N, size = N, replace = TRUE, prob = weights[t - 1, ])
    particles[t - 1, ] <- particles[t - 1, resample_indices]
    
    # Propagate particles
    particles[t, ] <- rnorm(N, mean = particles[t - 1, ], sd = sigma_proc)
    
    # Weight if observation is available
    if (!is.na(y[t])) {
      weights[t, ] <- dnorm(y[t], mean = particles[t, ], sd = sigma_obs)
      weights[t, ] <- weights[t, ] / sum(weights[t, ])
    } else {
      # No observation: keep uniform weights
      weights[t, ] <- rep(1/N, N)
    }
    
    # Estimate
    estimates[t] <- sum(particles[t, ] * weights[t, ])
  }
  
  return(list(
    particles = particles,
    weights = weights,
    estimates = estimates
  ))
}

# Use your calibration dataset with NA from 2024-01-01
y <- cal_ddat$`Streamflow Ave`
dates <- cal_ddat$Date

# Run particle filter
pf <- particle_filter_forecast(y, N = 1000, sigma_proc = sigma_proc, sigma_obs = sigma_obs)

# Plot result
plot(dates, ddat$`Streamflow Ave`, type = "l", col = "black", ylab = "Streamflow", xlab = "Date")
lines(dates, pf$estimates, col = "blue")
abline(v = as.Date("2024-01-01"), col = "red", lty = 2)
legend("topright", legend = c("Observed", "Particle Filter Estimate", "Forecast Start"), 
       col = c("black", "blue", "red"), lty = c(1, 1, 2))

# Get quantiles from particle cloud
ci <- apply(pf$particles, 1, quantile, probs = c(0.025, 0.5, 0.975))

# Plot with CI ribbon
plot(dates, ddat$`Streamflow Ave`, type = "l", col = "gray", ylab = "Streamflow", xlab = "Date")
polygon(c(dates, rev(dates)), 
        c(ci[1,], rev(ci[3,])), 
        col = adjustcolor("skyblue", alpha.f = 0.5), border = NA)
lines(dates, ci[2,], col = "blue", lwd = 2)
abline(v = as.Date("2024-01-01"), col = "red", lty = 2)
legend("topright", legend = c("Observed", "Filtered/Forecast Median", "95% CI", "Forecast Start"), 
       col = c("gray", "blue", "skyblue", "red"), lty = c(1, 1, NA, 2), pch = c(NA, NA, 15, NA))

