library(tidyverse)
library(introdataviz)
library("ggpubr")

# load data
ith_data_root <- "E:/LungIMCData/HumanWholeIMC/Aggregation/ITH"
lesion_ith_path <- file.path(ith_data_root, "lesion_ith_adc_10.csv")
lesion_ith <- read.csv(lesion_ith_path)

# update lesion ITH
lesion_ith <- lesion_ith %>% as_tibble() %>% 
    mutate(SmokeStatus = factor(SmokeStatus, levels = c("Never", "Heavy"))) %>%
    mutate(RawITH = 1.0 - RawITH) %>%
    mutate(RandomITH = 1.0 - RandomITH)

# ggscatter(lesion_ith, x = "RawITH", y = "RandomITH", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Raw ITH", ylab = "Random ITH")

## plot
# ggplot(data = lesion_ith, mapping = aes(x = RawITH, y = RandomITH, color=SmokeStatus)) +
ggplot(data = lesion_ith, mapping = aes(x = RawITH, y = RandomITH)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) + 
    xlim(0.1, 0.5) +
    ylim(0.1, 0.5) 