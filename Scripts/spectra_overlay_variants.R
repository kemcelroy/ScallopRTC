
library(ggplot2)

# Define which files to overlay and which to use as reference verticals
gaussian_ids <- c("Mut1_", "Mut2_", "Mut3_", "Mut4_") # Example: "Cfar-RTC_dG189S_", "Pyes-RTC_dS189G_", "Cfar-RTC_", "Pyes-RTC_"
ref_ids <- c("WT1_", "WT_") # Example: "Cfar-RTC_", "Pyes-RTC_"

# Store curves and vertical markers
curve_data <- list()
vlines <- data.frame(ID = character(), lambda = numeric(), stringsAsFactors = FALSE)

# Get Gaussian fits from references (for vertical lines)
for (ID in ref_ids) {
  result <- process_spectrum(ID, workpath = "./")
  fit <- fit_template_model(result$data, function(w, lmax) gaussian_template(w, lmax, width = 40), "Gaussian", result$peakWave)
  if (!is.null(fit)) {
    vlines <- rbind(vlines, data.frame(ID = ID, lambda = fit$lambda_max))
  }
}

# Get Gaussian curves for variants
for (ID in gaussian_ids) {
  result <- process_spectrum(ID, workpath = "./")
  fit <- fit_template_model(result$data, function(w, lmax) gaussian_template(w, lmax, width = 40), "Gaussian", result$peakWave)
  if (!is.null(fit)) {
    norm_pred <- fit$predicted / max(fit$predicted, na.rm = TRUE)
    curve_data[[ID]] <- data.frame(
      Wavelength = result$data$Wavelength,
      Absorbance = norm_pred,
      ID = ID
    )
  }
}

# Combine all curves into one data frame
all_curves <- do.call(rbind, curve_data)

# Plot it
p <- ggplot(all_curves, aes(x = Wavelength, y = Absorbance, color = ID)) +
  geom_line(size = 1) +
  geom_vline(data = vlines, aes(xintercept = lambda, linetype = ID), color = "black", show.legend = TRUE) +
  scale_linetype_manual(values = c("WT1_" = "dashed", "WT2_" = "dotted")) +
  coord_cartesian(xlim = c(415, 600)) +
  labs(title = "Normalized Gaussian Fits (450–500 nm)",
       x = "Wavelength (nm)", y = "Normalized Absorbance",
       color = "Variant", linetype = "Reference") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.ticks = element_line(color = "black"))

# Save to file
ggsave("overlay_variants_scaled.pdf", plot = p, width = 8, height = 6)
