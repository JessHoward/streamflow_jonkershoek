library(plotrix)
library(knitr)

##### Data validation
# Posterior predictions of x (latent state)
x_cols <- grep("^x\\[", colnames(out))
x_post <- exp(out[, x_cols])  # back-transform from log-scale
x_med <- apply(x_post, 2, median)
x_mean <- apply(x_post, 2, mean)

# Only for observed data
resids <- y[included] - x_med[included]

# Root Mean Squared Error
rmse <- sqrt(mean(resids^2))

# R-squared (squared correlation between observed and predicted)
r2 <- cor(y[included], x_med[included])^2

validation_metrics <- data.frame(
  Metric = c("RMSE", "R-squared"),
  Value = c(round(rmse, 3), round(r2, 3))
)

kable(validation_metrics, caption = "Model validation metrics: streamflow with rainfall, seasonality, and decay")


### Identify held out data
heldout.idx <- which(heldout)

### Extract modelled predictions for held out data
x.heldout.cols <- x.cols[heldout.idx]
ci.heldout <- apply(exp(out[, x.heldout.cols]), 2, quantile, c(0.025, 0.5, 0.975))


#### Compare predictions versus observed
validation_df <- data.frame(
  Date = time[heldout.idx],
  Observed = y_full[heldout.idx],         # Use full streamflow log here
  Predicted_median = ci.heldout[2,],
  Predicted_lower = ci.heldout[1,],
  Predicted_upper = ci.heldout[3,]
)

kable(head(validation_df, 10), digits = 3, caption = "Comparison of Model Predictions to Held-Out Observations")

# Remove any rows with missing values (just in case)
valid_data <- na.omit(validation_df)

# Extract observed and predicted vectors
obs <- valid_data$Observed
pred <- valid_data$Predicted_median

# Calculate R-squared
rsq <- cor(obs, pred)^2
rmse <- sqrt(mean((obs - pred)^2))

validation_metrics <- data.frame(
  Metric = c("RMSE", "R-squared"),
  Value = c(round(rmse, 3), round(rsq, 3))
)

kable(validation_metrics, caption = "Validation Metrics on Held-Out Data")

## outputs on non-log data

obs_real <- exp(valid_data$Observed)
pred_real <- valid_data$Predicted_median  # Already in real scale

rsq_real <- cor(obs_real, pred_real)^2
rmse_real <- sqrt(mean((obs_real - pred_real)^2))

validation_metrics_real <- data.frame(
  Metric = c("RMSE (real scale)", "R-squared (real scale)"),
  Value = c(round(rmse_real, 3), round(rsq_real, 3))
)

kable(validation_metrics_real, caption = "Validation on Held-Out Data (Real Scale)")


###Plotting outputs
# Add a tiny jitter to zero or near-zero values to avoid log(0)
jitter_amount <- 1e-6

valid_data <- validation_df %>%
  mutate(Observed_real = exp(Observed),
         Observed_jittered = ifelse(Observed_real <= 0, 1e-6, Observed_real),
         Predicted_jittered = ifelse(Predicted_median <= 0, 1e-6, Predicted_median),
         Observed_log = log10(Observed_jittered),
         Predicted_log = log10(Predicted_jittered))

valid_data$Predicted_log <- log10(valid_data$Predicted_median + 1e-6)

#linear scale streamflow
ggplot(valid_data, aes(x = Observed_real, y = Predicted_median)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Observed Streamflow (m³/s)",
    y = "Predicted Streamflow (m³/s)",
    title = "Observed vs Predicted Streamflow (Linear Scale)"
  ) +
  theme_minimal()

#log_scale
ggplot(valid_data, aes(x = Observed_log, y = Predicted_log)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "log10(Observed Streamflow)",
    y = "log10(Predicted Streamflow)",
    title = "Observed vs Predicted Streamflow (Log10 Scale)"
  ) +
  theme_minimal()


### Add Taylor plot
Observed <- validation_df$Observed
Predicted <- validation_df$Predicted_median
# Step 1: Clean / filter
qaqc <- complete.cases(Observed, Predicted)  # OR: qaqc <- !is.na(Observed) & !is.na(Predicted)

O <- Observed[qaqc]
E <- Predicted[qaqc]

# Step 2: Simulate ensemble matrix if you don't already have one
set.seed(123)
n_draws <- 100
jitter_sd <- sd(E - O, na.rm = TRUE)

stream <- replicate(n_draws, E + rnorm(length(E), 0, jitter_sd))

# Step 3: Taylor diagram
taylor.diagram(ref = O, model = E, normalize = TRUE, ref.sd = TRUE, col = "blue", pch = 16)

# Step 4: Add ensemble members
for(i in 1:ncol(stream)){
  taylor.diagram(ref = O, model = stream[, i],
                 col = rgb(1, 0, 0, 0.3), pch = ".", add = TRUE, normalize = TRUE)
}
