# Load libraries
library(tidyverse)
library(ggsignif)

# Clear working memory
rm(list = ls())

# Load data sets
HS0hr_1 <- read_tsv("./0hr_1.tsv")
HS0hr_2 <- read_tsv("./0hr_2.tsv")
HS37C30min_1 <- read_tsv("./37C30min_1.tsv")
HS37C30min_2 <- read_tsv("./37C30min_2.tsv")
HS37C2hr_1 <- read_tsv("./37C2hr_1.tsv")
HS37C2hr_2 <- read_tsv("./37C2hr_2.tsv")
HS37C4hr_1 <- read_tsv("./37C4hr_1.tsv")
HS37C4hr_2 <- read_tsv("./37C4hr_2.tsv")
HS42C30min_1 <- read_tsv("./42C30min_1.tsv")
HS42C30min_2 <- read_tsv("./42C30min_2.tsv")
HS42C2hr_1 <- read_tsv("./42C2hr_1.tsv")
HS42C2hr_2 <- read_tsv("./42C2hr_2.tsv")
HS42C4hr_1 <- read_tsv("./42C4hr_1.tsv")
HS42C4hr_2 <- read_tsv("./42C4hr_2.tsv")

# Definition of data-extracting function
cleanUpData <- function(x, y){
  cleaned_data <- x %>%
    mutate(Sample = y) %>%  
    select(Sample, CircRingAvgIntenRatioCh2)
  cleaned_data <- cleaned_data[1:200,]
  if(str_detect(y, pattern = "42") == TRUE){
    cleaned_data <- cleaned_data %>% 
      mutate(Type = "Heat shock")
  }else {
    cleaned_data <- cleaned_data %>% 
      mutate(Type = "Non-treat")
  }
  return(cleaned_data)
}

# Preparing data set containing all data
all_data <- cleanUpData(HS0hr_1, "0hr 1")
all_data <- rbind(all_data, cleanUpData(HS37C30min_1, "37C 0.5hr 1"))
all_data <- rbind(all_data, cleanUpData(HS37C2hr_1, "37C 2hr 1"))
all_data <- rbind(all_data, cleanUpData(HS37C4hr_1, "37C 4hr 1"))
all_data <- rbind(all_data, cleanUpData(HS42C30min_1, "42C 0.5hr 1"))
all_data <- rbind(all_data, cleanUpData(HS42C2hr_1, "42C 2hr 1"))
all_data <- rbind(all_data, cleanUpData(HS42C4hr_1, "42C 4hr 1"))
all_data2 <- cleanUpData(HS0hr_2, "0hr 2")
all_data2 <- rbind(all_data2, cleanUpData(HS37C30min_2, "37C 0.5hr 2"))
all_data2 <- rbind(all_data2, cleanUpData(HS37C2hr_2, "37C 2hr 2"))
all_data2 <- rbind(all_data2, cleanUpData(HS37C4hr_2, "37C 4hr 2"))
all_data2 <- rbind(all_data2, cleanUpData(HS42C30min_2, "42C 0.5hr 2"))
all_data2 <- rbind(all_data2, cleanUpData(HS42C2hr_2, "42C 2hr 2"))
all_data2 <- rbind(all_data2, cleanUpData(HS42C4hr_2, "42C 4hr 2"))
all_data2_2 <- cleanUpData(HS0hr_2, "0hr")
all_data2_2 <- rbind(all_data2_2, cleanUpData(HS42C30min_2, "42C 0.5hr"))
all_data2_2 <- rbind(all_data2_2, cleanUpData(HS42C2hr_2, "42C 2hr"))
all_data2_2 <- rbind(all_data2_2, cleanUpData(HS42C4hr_2, "42C 4hr"))

# Normalize data
all_data_tmp <- all_data %>% 
  filter(Sample == "0hr 1")
mean_of_alldata_0hr <- mean(all_data_tmp$CircRingAvgIntenRatioCh2)
all_data <- all_data %>% 
  mutate(CircRingAvgIntenRatioCh2_Normalized = CircRingAvgIntenRatioCh2 / mean_of_alldata_0hr) %>% 
  select(1, 3, 4)
all_data2_tmp <- all_data2 %>% 
  filter(Sample == "0hr 2")
mean_of_alldata2_0hr <- mean(all_data2_tmp$CircRingAvgIntenRatioCh2)
all_data2 <- all_data2 %>% 
  mutate(CircRingAvgIntenRatioCh2_Normalized = CircRingAvgIntenRatioCh2 / mean_of_alldata2_0hr) %>% 
  select(1, 3, 4)
all_data2_2_tmp <- all_data2_2 %>% 
  filter(Sample == "0hr")
mean_of_alldata2_2_0hr <- mean(all_data2_2_tmp$CircRingAvgIntenRatioCh2)
all_data2_2 <- all_data2_2 %>% 
  mutate(CircRingAvgIntenRatioCh2_Normalized = CircRingAvgIntenRatioCh2 / mean_of_alldata2_2_0hr) %>% 
  select(1, 3, 4)

# Modulate all_data2_2's x-label
all_data2_2 <- all_data2_2 %>% 
  mutate(Sample = gsub(Sample, pattern = "42C", replacement = "42°C", ignore.case = TRUE)) %>% 
  mutate(Sample = gsub(Sample, pattern = "hr", replacement = " hr", ignore.case = TRUE))

# Re-align the data order which would be appeared on the final graph
all_data2_2$Type <- factor(all_data2_2$Type, levels = c("Non-treat", "Heat shock"))

# Plot graphs (using all_data)
ggplot(all_data, aes(x = Sample, y = CircRingAvgIntenRatioCh2_Normalized, ymin = 0, colour = Type)) +
  ylim(0, 5) +
  geom_violin() +
  scale_color_manual(values = c("red", "skyblue")) +
  geom_boxplot(width = 0.1, fill = "black", outlier.colour = "NA") +
  stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  xlab("Samples") +
  ylab("Nuc/Cyto Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30))

# Plot graphs (using all_data2)
ggplot(all_data2, aes(x = Sample, y = CircRingAvgIntenRatioCh2_Normalized, ymin = 0, colour = Type)) +
  ylim(0, 5) +
  geom_violin() +
  scale_color_manual(values = c("red", "skyblue")) +
  geom_boxplot(width = 0.1, fill = "black", outlier.colour = "NA") +
  stat_summary(fun = mean, geom = "point", fill = "white", shape = 21, size = 2.5) +
  xlab("Samples") +
  ylab("Nuc/Cyto Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30))

# Plot graphs (using all_data2_2)
finalgraph <- ggplot(all_data2_2, aes(x = Sample, y = CircRingAvgIntenRatioCh2_Normalized, ymin = 0)) +
  ylim(0, 5) +
  geom_violin(colour = "gray30") +
  geom_boxplot(width = 0.1, fill = "black", outlier.colour = "NA") +
  stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_signif(inherit.aes = FALSE,
              stat = "identity",
              data = data.frame(x = c(0.8),
                                xend = c(3.2),
                                y = c(4.5),
                                annotation = c("bar1")),
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
pv <- c("***")
pg$data[[4]]$annotation <- pv
finalgraph <- ggplot_gtable(pg)
plot(finalgraph)

# Definition of ANOVA-Tukey-Cramer test function
tukey <- function(x){
  vx <- c(x$CircRingAvgIntenRatioCh2_Normalized)
  fx <- factor(x$Sample)
  anova_res <- aov(vx ~ fx)
  tukey_res <- TukeyHSD(anova_res)
  return(as.data.frame(tukey_res[[1]]))
}

# Perform Tukey-Cramer test
tukeytest_result_of_all_data <- tukey(all_data)
tukeytest_result_of_all_data2 <- tukey(all_data2)
tukeytest_result_of_all_data2_2 <- tukey(all_data2_2)

