# Load libraries
library(tidyverse)

# Clear workspace
rm(list = ls())

# Load datasets
PP1A_rawdata <- read.csv("./Values_FHF-PP1A.csv")
LATS2_rawdata <- read.csv("./Values_Myc-LATS2.csv")

# Process datasets
colnames(PP1A_rawdata) <- c("Distance", "PP1A")
colnames(LATS2_rawdata) <- c("Distance", "LATS2")
FinalData <- PP1A_rawdata %>% 
  left_join(LATS2_rawdata) %>% 
  pivot_longer(-Distance, names_to = "Protein", values_to = "Intensity")

# Re-align the data order which would be appeared on the final graph
FinalData$Protein <- factor(FinalData$Protein, levels = c("PP1A", "LATS2"))

# Plot curves
ggplot(FinalData, aes(x = Distance, y = Intensity, fill = Protein, colour = Protein)) +
  geom_line() +
  scale_color_manual(values = c("red", "green3")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        legend.title = element_blank(), 
        strip.background = element_blank())