library(ggplot2)
library(ggmagnify)

setwd("/.../sketchingplot/Linear")


Buff = as.integer(1e5)



df_Cov_Equi_tau1_NT = read.csv("SketchedNT/Solution/Figures/q1/tau1/Cov_ave.csv")
df_Cov_Equi_tau2_NT = read.csv("SketchedNT/Solution/Figures/q1/tau2/Cov_ave.csv")
df_Cov_Equi_tau3_NT = read.csv("SketchedNT/Solution/Figures/q1/tau3/Cov_ave.csv")
df_Cov_Equi_tau4_NT = read.csv("SketchedNT/Solution/Figures/q1/tau4/Cov_ave.csv")
df_Cov_Equi_tau5_NT = read.csv("SketchedNT/Solution/Figures/q1/tau5/Cov_ave.csv")

df_Err_SC_Equi_tau = data.frame(
  tau1 = df_Cov_Equi_tau1_NT$ErrSC,
  tau2 = df_Cov_Equi_tau2_NT$ErrSC,
  tau3 = df_Cov_Equi_tau3_NT$ErrSC,
  tau4 = df_Cov_Equi_tau4_NT$ErrSC,
  tau5 = df_Cov_Equi_tau5_NT$ErrSC
)



Err_p1 = ggplot(df_Err_SC_Equi_tau, aes(x = log(Buff + seq_along(tau1)))) +
  geom_line(aes(y = log(tau1)), color = "purple", linewidth = 1.2) +
  geom_line(aes(y = log(tau2)), color = "darkorange", linewidth = 1.2) +
  geom_line(aes(y = log(tau3)), color = "green", linewidth = 1.2) + 
  geom_line(aes(y = log(tau4)), color = "red", linewidth = 1.2) + 
  geom_line(aes(y = log(tau5)), color = "blue", linewidth = 1.2) + 
  geom_vline(
    xintercept = 12.0,
    linetype = "dashed",
    color = "magenta",
    linewidth = 1.2
  ) +
  labs(
    x = "log(iteration t)",
    y = "log(relative error)"
  ) +
  coord_cartesian(ylim = c(-1.5, 0.1)) +
  theme(
    axis.text.x = element_text(family = "serif", size = 18, color = "black"),
    axis.text.y = element_text(family = "serif", size = 18, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )



ggsave("Err_p1.pdf", plot = Err_p1, width = 9, height = 6, dpi = 200)

