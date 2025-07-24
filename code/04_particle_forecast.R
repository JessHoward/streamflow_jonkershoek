# # Assume: cal_ddat is already loaded and contains 'Date' and 'Streamflow Ave'
# 
# # Define constants
# forecast_horizon <- 1  # how far ahead you want to forecast each time
# start_forecast_date <- as.Date("2024-01-01")  # date to start iterative forecasts
# 
# # Storage for forecasts
# forecast_dates <- cal_ddat$Date[cal_ddat$Date >= start_forecast_date]
# n_forecasts <- length(forecast_dates)
# forecast_mean <- rep(NA, n_forecasts)
# forecast_ci <- matrix(NA, nrow = n_forecasts, ncol = 2) # 95% CI
# 
# for (i in seq_along(forecast_dates)) {
#   print(c(i/length(forecast_dates))*100)
#   
#   # Get data up to the current forecast point
#   forecast_time <- forecast_dates[i]
#   dat_sub <- cal_ddat |> filter(Date <= forecast_time)
#   
#   y <- dat_sub$`Streamflow Ave`
#   
#   # Create NA for the point you're forecasting
#   y[length(y)] <- NA
#   
#   # Set up JAGS model input
#   data <- list(
#     y = y,
#     n = length(y),
#     x_ic = 0.5,
#     tau_ic = 100,
#     a_obs = 1,
#     r_obs = 1,
#     a_add = 1,
#     r_add = 1
#   )
#   
#   # Initial values
#   init <- list()
#   for(j in 1:nchain){
#     y.samp <- sample(y, length(y), replace = TRUE)
#     init[[j]] <- list(
#       tau_add = 1 / var(diff(na.omit(y.samp))),
#       tau_obs = 5 / var(na.omit(y.samp))
#     )
#   }
#   
#   # Fit model
#   j.model <- jags.model(file = textConnection(RandomWalk),
#                         data = data,
#                         inits = init,
#                         n.chains = nchain,
#                         quiet = TRUE)
#   
#   update(j.model, 1000) # burn-in
#   
#   # Sample latent state
#   jags.out <- coda.samples(model = j.model,
#                            variable.names = "x",
#                            n.iter = 1000)
#   
#   # Extract posterior for last time step
#   xmat <- as.matrix(jags.out)
#   x_last <- xmat[, grep(paste0("x\\[", length(y), "\\]"), colnames(xmat))]
#   
#   forecast_mean[i] <- mean(x_last)
#   forecast_ci[i, ] <- quantile(x_last, probs = c(0.025, 0.975))
# }
# 
# # Combine forecasts with dates
# forecast_results <- tibble(
#   Date = forecast_dates,
#   Forecast = forecast_mean,
#   Lower = forecast_ci[,1],
#   Upper = forecast_ci[,2]
# )
# 
# # Plot
# ggplot() +
#   geom_line(data = cal_ddat, aes(x = Date, y = `Streamflow Ave`), color = "black") +
#   geom_line(data = forecast_results, aes(x = Date, y = Forecast), color = "blue") +
#   geom_ribbon(data = forecast_results,
#               aes(x = Date, ymin = Lower, ymax = Upper),
#               fill = "blue", alpha = 0.3) +
#   labs(title = "Iterative Forecasting", y = "Streamflow")
# 
# 


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
