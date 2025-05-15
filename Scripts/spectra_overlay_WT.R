library(ggplot2)
library(patchwork)

# List of IDs to plot (must have *_dark.TXT and *_blue.TXT files)
overlay_ids <- c("Airr-RTC_", "Apur-RTC_", "Aple-RTC_", "Pmax-RTC_", "Cfar-RTC_", "Pyes-RTC_")

# Collect spectra
spectra_list <- list()

for (ID in overlay_ids) {
  result <- process_spectrum(ID, workpath = "./")
  if (!is.null(result)) {
    spectra_list[[ID]] <- data.frame(
      Wavelength = result$data$Wavelength,
      Absorbance = result$data$Renorm_Absorbance,
      ID = ID
    )
  }
}

# Combine into one data frame
all_spectra <- do.call(rbind, spectra_list)

# Map your species names
species_labels <- c(
  "Airr-RTC_" = expression(italic("Ar. irradians")),
  "Apur-RTC_" = expression(italic("Ar. purpuratus")),
  "Aple-RTC_" = expression(italic("Am. pleuronectes")),
  "Pmax-RTC_" = expression(italic("P. maximus")),
  "Cfar-RTC_" = expression(italic("C. farreri")),
  "Pyes-RTC_" = expression(italic("M. yessoensis"))
)

# Plot
p <- ggplot(all_spectra, aes(x = Wavelength, y = Absorbance,
                             linetype = ID, color = ID)) +
  geom_line(size = 1) +
  scale_linetype_manual(values = c(
    "Airr-RTC_" = "solid",
    "Apur-RTC_" = "dashed",
    "Aple-RTC_" = "solid",
    "Pmax-RTC_" = "dashed",
    "Cfar-RTC_" = "dotted",
    "Pyes-RTC_" = "twodash"
  ), labels = species_labels) +
  scale_color_manual(values = c(
    "Airr-RTC_" = "black",
    "Apur-RTC_" = "black",
    "Aple-RTC_" = "gray40",
    "Pmax-RTC_" = "gray40",
    "Cfar-RTC_" = "black",
    "Pyes-RTC_" = "black"
  ), labels = species_labels) +
  guides(color = guide_legend(title = "Species"),
         linetype = guide_legend(title = "Species")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(400, 600)) +
  labs(x = "Wavelength (nm)", y = "Normalized Dark-Light Absorbance") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.ticks = element_line(color = "black"),
        plot.title = element_blank(),
        legend.text = element_text(face = "italic"))

# --- Function to approximate visible spectrum color based on wavelength (nm) ---
wavelength_to_hex <- function(wavelength) {
  if (wavelength >= 380 & wavelength <= 440) {
    R <- -(wavelength - 440) / (440 - 380)
    G <- 0
    B <- 1
  } else if (wavelength > 440 & wavelength <= 490) {
    R <- 0
    G <- (wavelength - 440) / (490 - 440)
    B <- 1
  } else if (wavelength > 490 & wavelength <= 510) {
    R <- 0
    G <- 1
    B <- -(wavelength - 510) / (510 - 490)
  } else if (wavelength > 510 & wavelength <= 580) {
    R <- (wavelength - 510) / (580 - 510)
    G <- 1
    B <- 0
  } else if (wavelength > 580 & wavelength <= 645) {
    R <- 1
    G <- -(wavelength - 645) / (645 - 580)
    B <- 0
  } else if (wavelength > 645 & wavelength <= 780) {
    R <- 1
    G <- 0
    B <- 0
  } else {
    R <- G <- B <- 0
  }
  
  # Intensity correction
  if (wavelength >= 380 & wavelength <= 419) {
    factor <- 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
  } else if (wavelength >= 420 & wavelength <= 700) {
    factor <- 1
  } else if (wavelength > 700 & wavelength <= 780) {
    factor <- 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
  } else {
    factor <- 0
  }
  
  rgb <- rgb(R * factor, G * factor, B * factor)
  return(rgb)
}

# --- Build wavelength-color data frame ---
spectrum_df <- data.frame(
  Wavelength = seq(380, 780, by = 1)
)
spectrum_df$Color <- sapply(spectrum_df$Wavelength, wavelength_to_hex)

# --- Spectrum bar plot ---
spectrum_bar <- ggplot(spectrum_df, aes(x = Wavelength, y = 1, fill = Color)) +
  geom_tile() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(410, 600)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# --- Combine with your existing plot (named `p`) ---
combined_plot <- p / spectrum_bar + plot_layout(heights = c(10, 1))

# --- Save ---
ggsave("darklight_overlay_all_WT_with_physical_bar.pdf", plot = combined_plot, width = 8, height = 4)
