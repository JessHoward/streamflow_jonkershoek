# # Assume: cal_ddat loaded and the data wrangling has been done as per streamflow analysis (ie. data has been removed from before 2013-03-16)

# Extract posterior samples
post <- as.matrix(jags.out)

# Convert precisions to standard deviations
sigma_proc <- 1 / sqrt(median(post[, "tau_add"]))
sigma_obs  <- 1 / sqrt(median(post[, "tau_obs"]))


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

