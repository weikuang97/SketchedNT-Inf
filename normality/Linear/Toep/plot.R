library(ggplot2)
library(ggmagnify)
library(latex2exp)


setwd("/.../normality/Linear/Toep")

df_diff_asgd = read.csv("ASGD/Solution/Figures/Diff_mat.csv")
df_diff_sgd1 = read.csv("SGD1/Solution/Figures/Diff_mat.csv")
df_diff_sgd2 = read.csv("SGD2/Solution/Figures/Diff_mat.csv")
df_diff_nt1 = read.csv("SketchedNT1/Solution/Figures/tau0/Diff_mat.csv")
df_diff_nt2 = read.csv("SketchedNT2/Solution/Figures/tau0/Diff_mat.csv")

std_normal = rnorm(200)

df_all = data.frame(
  ASGD        = sort(df_diff_asgd$diff_vec),
  SGD1        = sort(df_diff_sgd1$diff_vec),
  SGD2        = sort(df_diff_sgd2$diff_vec),
  SketchedNT1 = sort(df_diff_nt1$diff_vec),
  SketchedNT2 = sort(df_diff_nt2$diff_vec),
  ASGD_std    = sort(df_diff_asgd$diff_std_vec),
  SketchedNT1_std = sort(df_diff_nt1$diff_std_vec),
  SketchedNT2_std = sort(df_diff_nt2$diff_std_vec),
  std_normal  = sort(std_normal)
)

ggplot(df_all, aes(x = SketchedNT2, y = SGD2)) +
  geom_point(size = 1, color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.8) +
#  geom_abline(intercept = 0, slope = 1/sqrt(2), color = "green", linewidth = 0.8) +
  coord_fixed(1, xlim = c(-5, 5), ylim = c(-5, 5)) +                 # 1:1
  theme_bw() +
  labs(
    title = "",
    x = TeX("Newton $\\sqrt{1/\\bar{\\alpha}_t}\\cdot 1^T (x_t - x^*)$ with $\\beta_t = 1/t$"),
    y = TeX("SGD $\\sqrt{1/{\\alpha}_t}\\cdot 1^T (x_t - x^*)$ with $\\alpha_t = 1/t$")
  )


