# Input data files should be named as "1_WT_0hr.tsv".
# "1_" is the order of each sample presented in final graph.
# The mean value of the 1st sample will be used as a standard value for normalization.
# The file name of each sample will be applied to a x-label of final graph.
# The underscores in each file name will be converted into spaces when they are presented in the final graph as a x-label.
# You should change several VARIABLES among this script according to each experimental design.
# In case of U2OS-GFP-YAP1 cells, please comment out the line in 'Definition of data-extracting function'.

# Load libraries
library(tidyverse)
library(ggsignif)

# Clear working memory
rm(list = ls())

# Load data sets
files <- list.files(pattern = "\\.tsv$")
for (f in 1:length(files)){
  i <- gsub("(.*).tsv","\\1",files[f])
  if(f == 1){
    variable_name <- gsub("^[0-9]_", "", i)
  }else {
    variable_name <- c(variable_name, gsub("^[0-9]_", "", i))
  }
  if(str_detect(variable_name[f], pattern = "^[0-9]") == TRUE){
    variable_name[f] <- paste("tempprefix", variable_name[f], sep = "")
  }
  load_data <- paste(variable_name[f], " <- read_tsv(\"", files[f], "\")", sep = "")
  eval(parse(text = load_data))
}

# Definition of data-extracting function
cleanUpData <- function(x, y){
  cleaned_data <- x %>%
    filter(CircAvgIntenCh2 >= 300) %>% # This step must be commented out in case of U2OS-GFP-YAP1 cells
    mutate(Sample = y) %>%  
    select(Sample, CircRingAvgIntenRatioCh2)
  cleaned_data <- cleaned_data[1:200,]
  if(str_detect(y, pattern = "YF") == TRUE){ # VARIABLE-1
    cleaned_data <- cleaned_data %>% 
      mutate(Type = "YF mutant") # VARIABLE-2
  }else {
    cleaned_data <- cleaned_data %>% 
      mutate(Type = "Wild type") # VARIABLE-3
  }
  return(cleaned_data)
}

# Prepare data set containing all data
for (f in 1:length(files)){
  xlabel_name <- gsub("_", " ", variable_name[f])
  if(str_detect(xlabel_name, pattern = "^tempprefix") == TRUE){
    xlabel_name <- gsub("^tempprefix", "", xlabel_name)
  }
  if (f == 1){
    standard_sample <- xlabel_name
    prepare_data_set <- paste("all_data <- cleanUpData(", variable_name[f], ", \"", xlabel_name, "\")", sep = "")
  }else {
    prepare_data_set <- paste("all_data <- rbind(all_data, cleanUpData(", variable_name[f], ", \"", xlabel_name, "\"))", sep = "")
  }
  eval(parse(text = prepare_data_set))
}

# Normalize data
all_data_tmp <- all_data %>% 
  filter(Sample == standard_sample)
mean_of_alldata_0hr <- mean(all_data_tmp$CircRingAvgIntenRatioCh2)
all_data <- all_data %>% 
  mutate(CircRingAvgIntenRatioCh2_Normalized = CircRingAvgIntenRatioCh2 / mean_of_alldata_0hr) %>% 
  select(1, 3, 4) %>% 
  filter(Sample != "NA")

# Modulate all_data's x-label
all_data <- all_data %>% 
  mutate(Sample = gsub(Sample, pattern = "hr", replacement = " hr", ignore.case = TRUE))

# Plot graphs (using all_data)
finalgraph <- ggplot(all_data, aes(x = Sample, y = CircRingAvgIntenRatioCh2_Normalized, ymin = 0)) +
  ylim(0, 4) +
  geom_violin(colour = "gray30") +
  geom_boxplot(width = 0.1, 
               fill = "black", 
               outlier.colour = "NA") +
  stat_summary(fun = median, geom = "point", 
               fill = "white", 
               shape = 21, 
               size = 2.5) +
  geom_signif(inherit.aes = FALSE,
              stat = "identity",
              data = data.frame(x = c(0.8, 2.8),
                                xend = c(2.2, 4.2),
                                y = c(3.8),
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
  vx <- c(x$CircRingAvgIntenRatioCh2_Normalized)
  fx <- factor(x$Sample)
  anova_res <- aov(vx ~ fx)
  tukey_res <- TukeyHSD(anova_res)
  return(as.data.frame(tukey_res[[1]]))
}

# Perform Tukey-Cramer test
tukeytest_result_of_all_data <- tukey(all_data)

