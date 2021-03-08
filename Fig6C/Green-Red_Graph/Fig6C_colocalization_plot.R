# Load libraries
library(tidyverse)

# Clear workspace
rm(list = ls())

# Load datasets
SRC_rawdata <- read.csv("./Values_HA-C-SRC.csv")
LATS2_rawdata <- read.csv("./Values_Myc-LATS2.csv")

# Process datasets
colnames(SRC_rawdata) <- c("Distance", "C-SRC")
colnames(LATS2_rawdata) <- c("Distance", "LATS2")
FinalData <- SRC_rawdata %>% 
  left_join(LATS2_rawdata) %>% 
  pivot_longer(-Distance, names_to = "Protein", values_to = "Intensity")

# Plot curves
ggplot(FinalData, aes(x = Distance, y = Intensity, fill = Protein, colour = Protein)) +
  geom_line() +
  scale_color_manual(values = c("red", "green3")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        legend.title = element_blank(), 
        strip.background = element_blank())
