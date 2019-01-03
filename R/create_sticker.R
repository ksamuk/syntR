# create the syntR hex sticker

#devtools::install_github("GuangchuangYu/hexSticker")
library("hexSticker")
library("tidyverse")


x <- c(1:19) + rnorm(19, sd = 0.4)
y <- c(1:6, 13:7, 14:19) + rnorm(19, sd = 0.4)
clust_color <- c(rep("blue", 6), rep("red", 7), rep("green", 6))

df <- data.frame(x, y, clust_color)

sub_plot <- df %>%
  ggplot(aes(x = x, y = y, color = clust_color))+
  geom_point(size = 2.35, alpha = 1.0, shape = 16) +
  theme_classic()+
  ylab("") + xlab("") +
  scale_color_brewer(type = "qual", palette = "Set1")+
  #theme_void() +
  theme_transparent() +
  theme(axis.text = element_blank(),
        legend.position = "none")

sticker(sub_plot,
        package = "syntR", p_size = 10, p_x = 1, p_y = 1.5,
        s_x = 0.9, s_y = 0.7, s_width = 1.2, s_height = 1.2,
        spotlight = FALSE, l_x = 1.0, l_y = 1.0, l_width = 5, l_height = 5, l_alpha = 0.8,
        h_fill = "white", h_color = "black",
        #h_fill = "#3366CC", h_color = "#003399",
        p_color = "black",
        filename = "inst/figures/logo.png")
