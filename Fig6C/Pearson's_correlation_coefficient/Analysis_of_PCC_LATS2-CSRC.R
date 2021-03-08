# Load libraries
library(tidyverse)

# Clear workspace
rm(list = ls())

# Load datafile
PCC_rawdata <- read.csv("./PCC.csv")

# Check loaded data
glimpse(PCC_rawdata)
summary(PCC_rawdata)

# Draw histogram
ggplot(PCC_rawdata, aes(x = PCC, xmin = 0, xmax = 1)) +
  geom_histogram(binwidth = 0.01, 
                 colour = "black", 
                 fill = "lightgray") +
  xlab("Pearson's correlation coefficient") +
  theme_bw()
