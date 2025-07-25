### Streamflow data cleaning and exploration 

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
      ~ ifelse(. < -50 | . == 0, NA, .)
    ),
    `Relative Humidity Max` = ifelse(`Relative Humidity Max` < -50 | `Relative Humidity Max` > 1000 | `Relative Humidity Max` == 0, NA, `Relative Humidity Max`),
    across(
      c(`Soil Moisture 30cm Max`, `Soil Moisture 20cm Max`, `Soil Moisture 10cm Max`,
        `Soil Moisture 30cm Min`, `Soil Moisture 20cm Min`, `Soil Moisture 10cm Min`),
      ~ ifelse(. < 1 | Date < as.Date("2015-01-01"), NA, .)
    ),
    `Streamflow Ave` = ifelse(`Streamflow Ave` < 0.01, NA, `Streamflow Ave`)
  )
## Potentially problematic values: max and min humidity values of 100

hdat <- hdat %>%
  mutate(
    across(
      c(`Air Temperature`),
      ~ ifelse(. < -50, NA, .)
    ),
    `Relative Humidity` = ifelse(`Relative Humidity` < -50 | `Relative Humidity` > 1000 | `Relative Humidity` == 0, NA, `Relative Humidity`)
  )

# Potentially problematic values: humidity values of 100

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

# save cleaned data
# write_csv(ddat, "data/data_daily_cleaned.csv")
# write_csv(hdat, "data/data_hourly_cleaned.csv")


