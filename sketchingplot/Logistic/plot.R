library(ggplot2)
library(ggmagnify)

setwd("/.../Regression/Plot/Logistic")

df_Cov_NT0 = read.csv("SketchedNT/Solution/Figures/tau0/Cov_ave.csv")
df_Cov_NT2 = read.csv("SketchedNT/Solution/Figures/tau2/Cov_ave.csv")


Buff = as.integer(1e5)



Err_p2 = ggplot(df_Cov_NT0, aes(x = log(Buff + seq_along(ErrPI)))) +
  geom_line(aes(y = log(ErrPI)), color = "darkorange", linewidth = 1.2) +
  geom_line(aes(y = log(ErrSC)), color = "green", linewidth = 1.2) + 
  labs(
    x = "log(iteration t)",
    y = "log(relative error)"
  ) +
  theme(
    axis.text.x = element_text(family = "serif", size = 23, color = "black"),
    axis.text.y = element_text(family = "serif", size = 23, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )



Err_p3 = ggplot(df_Cov_NT2, aes(x = log(Buff + seq_along(ErrPI)))) +
  geom_line(aes(y = log(ErrPI)), color = "darkorange", linewidth = 1.2) +
  geom_line(aes(y = log(ErrSC)), color = "green", linewidth = 1.2) + 
  labs(
    x = "log(iteration t)",
    y = "log(relative error)"
  ) +
  theme(
    axis.text.x = element_text(family = "serif", size = 23, color = "black"),
    axis.text.y = element_text(family = "serif", size = 23, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )

# zoom-in setup
from <- c(xmin = 12.05, xmax = 12.25, ymin = -0.9, ymax = -0.5)
to <- c(xmin = 11.85, xmax = 12.55, ymin = 0, ymax = 1.3)
Err_p3 = Err_p3 + geom_magnify(
  from = from, 
  to = to, 
  axes = "xy",
  inset.linetype = "dashed",   
  target.linetype = "dashed",   
  proj.linetype = "dotted",
  colour = "darkblue",
  linewidth = 0.8
  )



Cov_p2 = ggplot(df_Cov_NT0, aes(x = (Buff + seq_along(CovPI))/1e5))+
  geom_line(aes(y = CovPI), color = "darkorange", linewidth = 1.2) +
  geom_line(aes(y = CovSC), color = "green", linewidth = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 1.8) +  # Horizontal line at y = 0.95
  scale_y_continuous(
    limits = c(0, 1),            
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)  
  ) + 
  labs(
    x = expression("iteration t" ~ "(\u00d7" * 10^5 * ")"),
    y = "coverage rate"
  ) +
  theme(
    axis.text.x = element_text(family = "serif", size = 23, color = "black"),
    axis.text.y = element_text(family = "serif", size = 23, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )

Cov_p3 = ggplot(df_Cov_NT2, aes(x = (Buff + seq_along(CovPI))/1e5))+
  geom_line(aes(y = CovPI), color = "darkorange", linewidth = 1.2) +
  geom_line(aes(y = CovSC), color = "green", linewidth = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 1.8) +  # Horizontal line at y = 0.95
  scale_y_continuous(
    limits = c(0, 1),           
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)  
  ) + 
  labs(
    x = expression("iteration t" ~ "(\u00d7" * 10^5 * ")"),
    y = "coverage rate"
  ) +
  theme(
    axis.text.x = element_text(family = "serif", size = 23, color = "black"),
    axis.text.y = element_text(family = "serif", size = 23, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )




Orc_p2 = ggplot(df_Cov_NT0, aes(x = (Buff + seq_along(CovXistar))/1e5))+
  geom_line(aes(y = CovXistar), color = "violet", linewidth = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 1.8) +  # Horizontal line at y = 0.95
  scale_y_continuous(
    limits = c(0, 1),            
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)  
  ) + 
  labs(
    x = expression("iteration t" ~ "(\u00d7" * 10^5 * ")"),
    y = "coverage rate"
  ) +
  theme(
    axis.text.x = element_text(family = "serif", size = 23, color = "black"),
    axis.text.y = element_text(family = "serif", size = 23, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )

Orc_p3 = ggplot(df_Cov_NT2, aes(x = (Buff + seq_along(CovXistar))/1e5))+
  geom_line(aes(y = CovXistar), color = "violet", linewidth = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth = 1.8) +  # Horizontal line at y = 0.95
  scale_y_continuous(
    limits = c(0, 1),           
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)  
  ) + 
  labs(
    x = expression("iteration t" ~ "(\u00d7" * 10^5 * ")"),
    y = "coverage rate"
  ) +
  theme(
    axis.text.x = element_text(family = "serif", size = 23, color = "black"),
    axis.text.y = element_text(family = "serif", size = 23, color = "black"),
    axis.title.x = element_text(family = "serif", size = 40, color = "black"),
    axis.title.y = element_text(family = "serif", size = 40, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),  # White background inside the plot area
    plot.background = element_rect(fill = "white", color = NA),  # White background for the entire plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines (light gray)
    panel.grid.minor = element_line(color = "gray90"),           # Minor grid lines (lighter gray)
    axis.line = element_line(color = "black")                    # Black axis lines for both x and y axes
  )



ggsave("Err_p2.png", plot = Err_p2, width = 9, height = 6, dpi = 200)
ggsave("Err_p3.png", plot = Err_p3, width = 9, height = 6, dpi = 200)
ggsave("Cov_p2.png", plot = Cov_p2, width = 9, height = 6, dpi = 200)
ggsave("Cov_p3.png", plot = Cov_p3, width = 9, height = 6, dpi = 200)
ggsave("Orc_p2.png", plot = Orc_p2, width = 9, height = 6, dpi = 200)
ggsave("Orc_p3.png", plot = Orc_p3, width = 9, height = 6, dpi = 200)
