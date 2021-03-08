# Load libraries
library(tidyverse)
library(ggsignif)

# Clear workspace
rm(list = ls())

# Load datasets
HS0hr_rawdata <- read.csv("./HS0hr.csv")
HS30min_rawdata <- read.csv("./HS30min.csv")
HS2hr_rawdata <- read.csv("./HS2hr.csv")
HS4hr_rawdata <- read.csv("./HS4hr.csv")

# Check loaded datasets
glimpse(HS0hr_rawdata)
glimpse(HS30min_rawdata)
glimpse(HS2hr_rawdata)
glimpse(HS4hr_rawdata)

# Extract exact data values from raw data
HS0hr_processed_1 <- HS0hr_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)
HS30min_processed_1 <- HS30min_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)
HS2hr_processed_1 <- HS2hr_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)
HS4hr_processed_1 <- HS4hr_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)

# Check histogram of each data
hist(HS0hr_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))
hist(HS30min_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))
hist(HS2hr_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))
hist(HS4hr_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))

# Perform Welch's t-test
t.test(HS0hr_processed_1$icn.factor, HS2hr_processed_1$icn.factor, var.equal = F)
t.test(HS2hr_processed_1$icn.factor, HS4hr_processed_1$icn.factor, var.equal = F)

# Prepare final data
final_data <- data.frame(HS0hr = HS0hr_processed_1$icn.factor, 
                        HS0.5hr = HS30min_processed_1$icn.factor, 
                        HS2hr = HS2hr_processed_1$icn.factor, 
                        HS4hr = HS4hr_processed_1$icn.factor) %>% 
  rowid_to_column("id") %>% 
  pivot_longer(-id, names_to = "HeatShock", values_to = "RelANucCytoRatio") %>% 
  select(-id) %>% 
  mutate(HeatShock = gsub(HeatShock, pattern = "HS", replacement = "42°C ", ignore.case = TRUE)) %>% 
  mutate(HeatShock = gsub(HeatShock, pattern = "hr", replacement = " hr", ignore.case = TRUE)) %>% 
  mutate(Type = if_else(condition = str_detect(HeatShock, pattern = "0 hr"), 
                         true = "Non-treat",
                         false = "Heat shock"))

# Re-align the data order which would be appeared on the final graph
final_data$HeatShock <- factor(final_data$HeatShock, levels = c("42°C 0 hr", "42°C 0.5 hr", "42°C 2 hr", "42°C 4 hr"))
final_data$Type <- factor(final_data$Type, levels = c("Non-treat", "Heat shock"))

# Prepare final plot
final_graph <- ggplot(final_data, aes(x = HeatShock, y = RelANucCytoRatio, ymin = 0)) +
  ylim(0, 1.3) +
  geom_violin(colour = "gray30") +
  geom_boxplot(width = 0.1, 
               fill = "black", 
               outlier.colour = "NA") +
  stat_summary(fun = median, 
               geom = "point", 
               fill = "white", 
               shape = 21, 
               size = 2.5) +
  geom_signif(inherit.aes = FALSE,
              stat = "identity",
              data = data.frame(x = c(0.8, 0.8),
                                xend = c(2.2, 3.2),
                                y = c(1.0, 1.2),
                                annotation = c("bar1", "bar2")),
              aes(x = x, xend = xend, y = y, yend = y, annotation = annotation),
              tip_length = 0) +
  xlab("") +
  ylab("RelA Nuc/Cyto Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        legend.title = element_blank(), 
        strip.background = element_blank())
pg <- ggplot_build(final_graph)
pv <- pg$data[[4]]$annotation
pv <- c("n.s.", "***")
pg$data[[4]]$annotation <- pv
final_graph <- ggplot_gtable(pg)
plot(final_graph)

# Definition of ANOVA-Tukey-Cramer test function
tukey <- function(x){
  vx <- c(x$RelANucCytoRatio)
  fx <- factor(x$HeatShock)
  anova_res <- aov(vx ~ fx)
  tukey_res <- TukeyHSD(anova_res)
  return(as.data.frame(tukey_res[[1]]))
}

# Perform Tukey-Cramer test
tukeytest_result_of_final_data <- tukey(final_data)
