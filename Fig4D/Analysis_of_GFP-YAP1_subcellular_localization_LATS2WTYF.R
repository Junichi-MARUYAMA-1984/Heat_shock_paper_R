# Load libraries
library(tidyverse)
library(ggsignif)

# Clear workspace
rm(list = ls())

# Load datasets
HS0hrWT_rawdata <- read.csv("./combined_csv_file_WT0hr.csv")
HS2hrWT_rawdata <- read.csv("./combined_csv_file_WT2hr.csv")
HS0hrYF_rawdata <- read.csv("./combined_csv_file_YF0hr.csv")
HS2hrYF_rawdata <- read.csv("./combined_csv_file_YF2hr.csv")

# Check loaded datasets
glimpse(HS0hrWT_rawdata)
glimpse(HS2hrWT_rawdata)
glimpse(HS0hrYF_rawdata)
glimpse(HS2hrYF_rawdata)

# Extract exact data values from raw data
HS0hrWT_processed_1 <- HS0hrWT_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)
HS2hrWT_processed_1 <- HS2hrWT_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)
HS0hrYF_processed_1 <- HS0hrYF_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)
HS2hrYF_processed_1 <- HS2hrYF_rawdata %>% 
  select(icn.factor, X..nuclei, X..cytoplasm, av..nuclei.intensity, av..cytoplasm.intensity)

# Slice the HS0hrWT/HS0hrYF/HS2hrYF data due to the element number of HS2hrWT (39)
HS0hrWT_processed_1 <- HS0hrWT_processed_1 %>% 
  slice(1:39)
HS0hrYF_processed_1 <- HS0hrYF_processed_1 %>% 
  slice(1:39)
HS2hrYF_processed_1 <- HS2hrYF_processed_1 %>% 
  slice(1:39)

# Check histogram of each data
hist(HS0hrWT_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))
hist(HS2hrWT_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))
hist(HS0hrYF_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))
hist(HS2hrYF_processed_1$icn.factor, breaks = seq(0, 2.0, 0.1))

# Perform Welch's t-test
t.test(HS0hrWT_processed_1$icn.factor, HS2hrWT_processed_1$icn.factor, var.equal = F)
t.test(HS0hrYF_processed_1$icn.factor, HS2hrYF_processed_1$icn.factor, var.equal = F)

# Prepare final data
FinalData <- data.frame(WT0hr = HS0hrWT_processed_1$icn.factor, 
                        WT2hr = HS2hrWT_processed_1$icn.factor, 
                        YF0hr = HS0hrYF_processed_1$icn.factor,
                        YF2hr = HS2hrYF_processed_1$icn.factor) %>% 
  rowid_to_column("id") %>% 
  pivot_longer(-id, names_to = "HeatShock", values_to = "YAP1NucCytoRatio") %>% 
  select(-id) %>% 
  mutate(HeatShock = gsub(HeatShock, pattern = "WT", replacement = "WT ", ignore.case = TRUE)) %>%
  mutate(HeatShock = gsub(HeatShock, pattern = "YF", replacement = "YF ", ignore.case = TRUE)) %>% 
  mutate(HeatShock = gsub(HeatShock, pattern = "hr", replacement = " hr", ignore.case = TRUE)) %>% 
  mutate(Type = "NA")
for(i in 1:nrow(FinalData)){
  if(str_detect(FinalData$HeatShock[i], pattern = "YF") == TRUE){
    FinalData$Type[i] = "YF mutant"
  }else {
    FinalData$Type[i] = "Wild type"
  }
}

# Re-align the final data for plotting
FinalData$HeatShock <- factor(FinalData$HeatShock, levels = c("WT 0 hr", "WT 2 hr", "YF 0 hr", "YF 2 hr"))

# Prepare final plot
finalgraph <- ggplot(FinalData, aes(x = HeatShock, y = YAP1NucCytoRatio, ymin = 0)) +
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
              data = data.frame(x = c(0.8, 2.8),
                                xend = c(2.2, 4.2),
                                y = c(1.2, 1.2),
                                annotation = c("bar1", "bar2")),
              aes(x = x, xend = xend, y = y, yend = y, annotation = annotation),
              tip_length = 0) +
  xlab("") +
  ylab("YAP1 Nuc/Cyto Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        legend.title = element_blank(), 
        strip.background = element_blank())
pg <- ggplot_build(finalgraph)
pv <- pg$data[[4]]$annotation
pv <- c("***", "***")
pg$data[[4]]$annotation <- pv
finalgraph <- ggplot_gtable(pg)
plot(finalgraph)

# Definition of ANOVA-Tukey-Cramer test function
tukey <- function(x){
  vx <- c(x$YAP1NucCytoRatio)
  fx <- factor(x$HeatShock)
  anova_res <- aov(vx ~ fx)
  tukey_res <- TukeyHSD(anova_res)
  return(as.data.frame(tukey_res[[1]]))
}

# Perform Tukey-Cramer test
tukeytest_result_of_FinalData <- tukey(FinalData)

