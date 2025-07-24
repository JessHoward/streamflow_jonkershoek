### Covariate identification

# Load packages
library(tidyverse)

# Read daily stramflow data
ddat <- read_csv("data/data_daily_cleaned.csv")

# Create new variable - Rainfall Total "set back" one day - assumes there will be a one day
# lag between rain falling and streamflow increase
ddat$rainfall_dayback = ddat$`Rainfall Total`[c(2:length(ddat$`Rainfall Total`),NA)]

library(GGally)
ggpairs(ddat %>% select(-Date))
