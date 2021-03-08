# Load libraries
library(tidyverse)
library(ggsignif)

# Clear working memory
rm(list = ls())

# Load data sets
HS0hr <- read_tsv("./0hr.tsv")
DMSO2hr <- read_tsv("./DMSO2hr.tsv")
DMSO4hr <- read_tsv("./DMSO4hr.tsv")
Herb2hr <- read_tsv("./Herb2hr.tsv")
Herb4hr <- read_tsv("./Herb4hr.tsv")

# Definition of data-extracting function
cleanUpData <- function(x, y){
  cleaned_data <- x %>%
    mutate(Sample = y) %>%  
    select(Sample, CircRingAvgIntenRatioCh2)
  cleaned_data <- cleaned_data[1:200,]
  if(str_detect(y, pattern = "HerbimycinA") == TRUE){
    cleaned_data <- cleaned_data %>% 
      mutate(Type = "HerbimycinA")
  }else {
    cleaned_data <- cleaned_data %>% 
      mutate(Type = "DMSO")
  }
  return(cleaned_data)
}

# Preparing data set containing all data
all_data <- cleanUpData(HS0hr, "0 hr")
all_data <- rbind(all_data, cleanUpData(DMSO2hr, "DMSO 2 hr"))
all_data <- rbind(all_data, cleanUpData(DMSO4hr, "DMSO 4 hr"))
all_data <- rbind(all_data, cleanUpData(Herb2hr, "HerbimycinA 2 hr"))
all_data <- rbind(all_data, cleanUpData(Herb4hr, "HerbimycinA 4 hr"))

# Normalize data
all_data_tmp <- all_data %>% 
  filter(Sample == "0 hr")
mean_of_alldata_0hr <- mean(all_data_tmp$CircRingAvgIntenRatioCh2)
all_data <- all_data %>% 
  mutate(CircRingAvgIntenRatioCh2_Normalized = CircRingAvgIntenRatioCh2 / mean_of_alldata_0hr) %>% 
  select(1, 3, 4)

# Plot graphs (using all_data)
finalgraph <- ggplot(all_data, aes(x = Sample, y = CircRingAvgIntenRatioCh2_Normalized, ymin = 0)) +
  ylim(0, 3) +
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
              data = data.frame(x = c(2.8),
                                xend = c(5.2),
                                y = c(2.8),
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

