library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)

# Modified version of plot_dark_light_with_inset() that returns a cowplot grob
build_dark_light_cowplot <- function(result, fits_diff) {
  colors <- c("Difference spectrum" = "black", "Gaussian model" = "gray")

  p_main <- ggplot(result$data, aes(x = Wavelength)) +
    geom_line(aes(y = Renorm_Absorbance, color = "Difference spectrum"), size = 0.8) +
    labs(title = paste(result$ID, "Difference Spectrum"),
         x = "Wavelength (nm)", y = "Normalized Absorbance") +
    coord_cartesian(xlim = c(435, NA)) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.75, 0.2),
      legend.title = element_blank(),
      #legend.background = element_rect(color = "black", fill = "white")
    )

  for (fit in fits_diff) {
    if (fit$name == "Gaussian") {
      fit_df <- data.frame(
        Wavelength = result$data$Wavelength,
        Pred = fit$predicted
      )
      p_main <- p_main +
        geom_line(data = fit_df, aes(x = Wavelength, y = Pred, color = "Gaussian model"), linewidth = 0.8)
    }
  }

  # Inset: dark/light spectra
  inset_data <- result$full_spectrum %>%
    dplyr::filter(Wavelength >= 400 & Wavelength <= 650) %>%
    dplyr::mutate(
      Absorbance_dark = Absorbance_dark - min(Absorbance_dark, na.rm = TRUE),
      Absorbance_blue = Absorbance_blue - min(Absorbance_blue, na.rm = TRUE)
    )
  inset_colors <- c("Dark" = "black", "Light" = "blue")

  p_inset <- ggplot(inset_data, aes(x = Wavelength)) +
    geom_line(aes(y = Absorbance_dark, color = "Dark"), linewidth = 0.7) +
    geom_line(aes(y = Absorbance_blue, color = "Light"), linewidth = 0.7) +
    scale_color_manual(values = inset_colors) +
    labs(x = "Wavelength (nm)", y = "Absorbance") +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.82, 0.45),
      legend.title = element_blank(),
      #legend.background = element_rect(color = "black", fill = "white")
    )

  # Vertical line
  gaussian_fit <- Filter(function(f) f$name == "Gaussian", fits_diff)
  if (length(gaussian_fit) == 1) {
    x <- gaussian_fit[[1]]$lambda_max
    vline_pos <- 0.978 * x + 7.3285
    
    # Create the fixed-position label grob
    label_grob <- textGrob(
      label = bquote(lambda[max] == .(round(vline_pos, 1)) ~ "nm"),
      x = unit(0.51, "npc"), y = unit(0.97, "npc"),
      just = c("left", "top"),
      gp = gpar(col = "red", fontsize = 6.5)
    )
    
    p_inset <- p_inset +
      geom_vline(xintercept = vline_pos, color = "red", linetype = "dashed") +
      annotation_custom(label_grob)
  }

  # Combine using cowplot
  combined <- ggdraw() +
    draw_plot(p_main) +
    draw_plot(p_inset, x = 0.48, y = 0.41, width = 0.50, height = 0.50)

  return(combined)
}

# Run batch and collect plots
plots <- list()
panel_ids <- c("SAMPLE1_", "SAMPLE2_", "SAMPLE3_") # Example: "Apur-RTC_dW188Y_", "Apur-RTC_dA193V_", "Apur-RTC_dW188Y_A193V_", "Cfar-RTC_dG189S_", "Pyes-RTC_dS189G_"

for (ID in panel_ids) {
  result <- process_spectrum(ID = ID, workpath = "./")
  fits_diff <- list(
    Gacic = fit_template_model(result$data, gacic_template_complex, "Gacic", result$peakWave),
    Govardovski = fit_template_model(result$data, govardovski_template_complex, "Govardovski", result$peakWave),
    Stavenga = fit_template_model(result$data, stavenga_template_complex, "Stavenga", result$peakWave),
    Gaussian = fit_template_model(result$data, function(w, lmax) gaussian_template(w, lmax, width = 40), "Gaussian", result$peakWave)
  )
  fits_diff <- Filter(Negate(is.null), fits_diff)
  plots[[length(plots) + 1]] <- build_dark_light_cowplot(result, fits_diff)
}

pdf("multi_panel_mutants.pdf", width = 8, height = 11)
grid.arrange(grobs = plots, nrow = 3, ncol = 2)
dev.off()
