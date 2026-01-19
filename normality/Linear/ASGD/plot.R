library(ggplot2)
library(ggmagnify)


setwd("/Users/weikuang/Desktop/UChicago/StoOpt/IMA/github/SketchedNT-Inf/normality/Linear")

df_diff = read.csv("ASGD/Solution/Figures/Diff_mat.csv")



ggplot(df_diff, aes(sample = diff_t1)) +
  stat_qq(color = "black", size = 1) +
  stat_qq_line(color = "red", linewidth = 0.8) +
  theme_bw() + 
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )
