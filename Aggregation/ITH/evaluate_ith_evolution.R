library(tidyverse)
library(introdataviz)

# load data
ith_data_root <- "E:/LungIMCData/HumanWholeIMC/Aggregation/ITH"
lesion_ith_path <- file.path(ith_data_root, "lesion_ith_entire_fea.csv")
lesion_ith <- read.csv(lesion_ith_path)

# update lesion ITH
lesion_ith <- lesion_ith %>% as_tibble() %>% 
    mutate(LesionStage = factor(LesionStage, levels = c("AAH", "AIS", "MIA", "ADC"))) %>%
    mutate(SmokeStatus = factor(SmokeStatus, levels = c("Never", "Heavy"))) %>%
    mutate(LesionITH = 1.0 - LesionITH) 

# compute Kruskal-Wallis test
kw_teset <- kruskal.test(LesionITH ~ LesionStage, data = lesion_ith)

# plot
ggplot(data = lesion_ith, mapping = aes(x = LesionStage, y = LesionITH, fill=SmokeStatus)) +
    geom_split_violin() +
    ylim(0.05, 0.45) +
    labs(title = paste("p-value:", kw_teset$p.value, "(Kruskal-Wallis test)"),
         x = "Pathological Stage",
         y = "Lesion Intratumoral Heterogeneity")

# # plot
# ggplot(data = lesion_ith, mapping = aes(x = LesionStage, y = LesionITH, fill=SmokeStatus)) +
#     geom_boxplot() +
#     ylim(0.05, 0.45)