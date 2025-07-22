### Streamflow data exploration 

# Load packages
library(tidyverse)

# Read daily metadata
dmdat <- read_csv("data/metadata_daily_2025-07-09.csv")

# Read hourly metadata
hmdat <- read_csv("data/metadata_hourly_2025-07-09.csv")

# Read daily data
ddat <- read_csv("data/data_daily_2025-07-09.csv")

# Read hourly data
hdat <- read_csv("data/data_hourly_2025-07-09.csv")

# Data cleaning 
summary(ddat)

# Convert unusual values NA
ddat <- ddat %>%
  mutate(
    across(
      c(`Air Temperature Max`, `Air Temperature Min`, `Relative Humidity Min`),
      ~ ifelse(. < -50, NA, .)
    ),
    `Relative Humidity Max` = ifelse(`Relative Humidity Max` < -50 | `Relative Humidity Max` > 1000, NA, `Relative Humidity Max`),
    across(
      c(`Soil Moisture 30cm Max`, `Soil Moisture 20cm Max`, `Soil Moisture 10cm Max`,
        `Soil Moisture 30cm Min`, `Soil Moisture 20cm Min`, `Soil Moisture 10cm Min`),
      ~ ifelse(. < 1 | Date < as.Date("2015-01-01"), NA, .)
    )
  )

hdat <- hdat %>%
  mutate(
    across(
      c(`Air Temperature`),
      ~ ifelse(. < -50, NA, .)
    ),
    `Relative Humidity` = ifelse(`Relative Humidity` < -50 | `Relative Humidity` > 1000, NA, `Relative Humidity`)
  )

# I think we should check the data again as there might still be some abnormal values I've missed

# Pivot longer and visualize the daily data
ddat |>
  pivot_longer(cols = -c(Date),
               names_to = "phenomenon",
               values_to = "value") |>
  ggplot(aes(x = Date, y = value)) +
  geom_line() +
  facet_wrap(~phenomenon, scales = "free") +
  labs(title = "Daily Data",
       x = "Date",
       y = "Value") +
  theme_minimal()

# Pivot longer and visualize the hourly data
hdat |>
  pivot_longer(cols = -c(Date),
               names_to = "phenomenon",
               values_to = "value") |>
  ggplot(aes(x = Date, y = value)) +
  geom_line() +
  facet_wrap(~phenomenon, scales = "free") +
  labs(title = "Hourly Data",
       x = "Date",
       y = "Value") +
  theme_minimal()

# Plot daily streamflow data with horizontal line showing flood threshold
ddat |> ggplot() +
  geom_line(aes(y = `Streamflow Ave`, x = as.Date(Date))) +
  geom_hline(aes(yintercept = 4.076)) +
  ggtitle("Langrivier daily streamflow") +
  xlab("Date") +
  ylab("Streamflow (Cubic metres per second)")

# Plot hourly streamflow data with horizontal line showing flood threshold
hdat |> ggplot() +
  geom_line(aes(y = `Streamflow`, x = as.Date(Date))) +
  geom_hline(aes(yintercept = 4.076)) +
  ggtitle("Langrivier hourly streamflow") +
  xlab("Date") +
  ylab("Streamflow (Cubic metres per second)")


# # Assume you have a vector of streamflows called 'streamflow'
# threshold <- 4.076 # Example threshold
# exceedance_probability <- 1 - pnorm(threshold, mean(streamflow$value, na.rm = T), sd(streamflow$value, na.rm = T))
# print(exceedance_probability)
