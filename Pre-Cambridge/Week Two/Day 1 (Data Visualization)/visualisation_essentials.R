library(ggplot2)

# useful ggplot extensions
library(ggrepel)
library(ggpubr)
library(patchwork)

# Specialised visualisations
library(EnhancedVolcano)
library(pheatmap)
library(ComplexHeatmap)

my_theme <- function(){
      theme_bw() +
      theme(axis.title = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12))
}

cars_df <- mtcars
cars_df$mpg[5] <- 50
flower_df <- iris

ggplot(flower_df, aes(Sepal.Length, Sepal.Width, colour = Species)) +
  geom_jitter(size = 5, alpha = 0.5, height = 0.1, width = 0.1) +
  scale_color_brewer(type = "qual", palette = 2) +
  my_theme()
  


set.seed(123)
ggplot(cars_df, aes(as.factor(am), mpg)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height = 0, width = 0.3) +
  my_theme() +
  xlab("Transmission type (0: automatic, 1: manual)") +
  ylab("Miles / gallon")


