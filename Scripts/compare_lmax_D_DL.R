# Load required packages
library(ggplot2)
library(ggpmisc)

# Create the data frame
df <- data.frame(
  lmaxDL = c(502.85, 497.96, 497.44, 506.21, 517.92, 514.84, 504.14, 514.25, 498.93, 497.86, 503.31),
  lmaxD  = c(500.38, 493.70, 494.88, 500.88, 512.13, 511.66, 499.72, 512.57, 494.21, 495.29, 498.71)
)

# Fit the model
fit <- lm(lmaxD ~ lmaxDL, data = df)

# Create the plot
splot <- ggplot(df, aes(x = lmaxDL, y = lmaxD)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", linetype = "dashed") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE,
    label.x = "right", label.y = "top") +
  labs(
    title = NULL,
    x = expression(lambda[max]~"(Dark-Light, Gaussian fit)"),
    y = expression(lambda[max]~"(Dark, Stavenga A1 pigment fit)")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

# Plot to screen
plot(splot)

# Save to PDF
ggsave("lambda_max_scatter.pdf", splot, width = 6, height = 5)
