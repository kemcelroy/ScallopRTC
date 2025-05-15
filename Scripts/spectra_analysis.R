library(readr)
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(grid)

# --- Template Functions ---

# --- Gacic Template ---
gacic_template_complex <- function(x, lambda_max) {
  A <- 27.45313; B <- 34.335; C <- 0.3809; D <- 3.6254
  term <- (x + B * log(C) - lambda_max) / B
  raw <- A * (1 + exp(term))^(-D) * exp(term) * C^(-D)
  
  # Normalize to peak = 1
  peak_term <- log(C)
  peak_value <- A * (1 + exp(peak_term))^(-D) * exp(peak_term) * C^(-D)
  
  return(raw / peak_value)
}

# --- Govardovski Template ---
govardovski_template_complex <- function(x, lambda_max) {
  term1 <- exp(69.7 * ((0.8795 + 0.0459 * exp(-1 * ((lambda_max - 300)^2) / 11940)) - (lambda_max / x)))
  term2 <- exp(28 * (0.922 - (lambda_max / x)))
  term3 <- exp(-14.9 * (1.104 - (lambda_max / x)))
  1 / (term1 + term2 + term3 + 0.674)
}

# --- Stavenga Template ---
stavenga_template_complex <- function(x, lambda_max) {
  width <- 1.0
  log_ratio1 <- log10(x / lambda_max)
  term1 <- exp(-width * 380 * log_ratio1^2 * (1 + 6.09 * log_ratio1 + 3 * (6.09^2) / 8 * log_ratio1^2))
  log_ratio2 <- log10(x / 350)
  term2 <- 0.29 * exp(-width * 247 * log_ratio2^2 * (1 + 3.59 * log_ratio2 + 3 * (3.59^2) / 8 * log_ratio2^2))
  term1 + term2
}

# --- Gaussian Template ---
gaussian_template <- function(wavelength, lambda_max, width = 20) {
  exp(-((wavelength - lambda_max)^2) / (2 * width^2))
}


# --- Header Detection ---
detect_header_end <- function(file) {
  lines <- readLines(file)
  start_line <- grep("^nm\\s+Abs", lines)
  return(start_line[1])
}

# --- Template Fitting ---
fit_template_model <- function(data, template_fn, template_name, start_lambda) {
  if (!"Wavelength" %in% names(data)) {
    stop("🚫 'Wavelength' column not found in data passed to fit_template_model()")
  }
  if (nrow(data) < 5) return(350)
  tryCatch({
    fit <- nlsLM(
      Renorm_Absorbance ~ a * template_fn(Wavelength, lambda_max),
      data = data,
      start = list(lambda_max = start_lambda, a = 1),
      lower = c(lambda_max = start_lambda - 25, a = 0.9),
      upper = c(lambda_max = start_lambda + 75, a = 1.1),
      control = nls.lm.control(maxiter = 1000)
    )
    predicted <- predict(fit)
    actual <- data$Renorm_Absorbance
    ss_res <- sum((actual - predicted)^2)
    ss_tot <- sum((actual - mean(actual))^2)
    r_squared <- 1 - ss_res / ss_tot
    return(list(
      name = template_name,
      fit = fit,
      predicted = predicted,
      lambda_max = coef(fit)[["lambda_max"]],
      r_squared = r_squared
    ))
  }, error = function(e) {
    message("⚠️ ", template_name, " fit failed: ", conditionMessage(e))
    return(NULL)
  })
}


# --- Plot Function for Dark-Light Spectrum Fit ---
plot_dark_light_fit <- function(filtered_data, fit_results, ID) {
  p <- ggplot(filtered_data, aes(x = Wavelength, y = Renorm_Absorbance)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    theme_minimal() +
    labs(title = paste(ID, "Dark - Blue Difference Fit"), x = "Wavelength (nm)", y = "Normalized Absorbance")
  
  colors <- c(Gacic = "darkorange", Govardovski = "darkgreen", Stavenga = "purple")
  for (fit in fit_results) {
    p <- p + geom_line(aes(y = fit$predicted), color = colors[fit$name], linewidth = 0.8)
  }
  print(p)
}

# --- Add Plot Function for Dark-Only with Inset ---
plot_dark_only_with_inset <- function(full_spectrum, dark_only_data, fit_results, peakWave, ID) {
  # Main plot: full spectrum (not normalized)
  main_plot <- ggplot(full_spectrum, aes(x = Wavelength)) +
    geom_line(aes(y = Absorbance_dark), color = "black", size = 0.8, alpha = 0.8) +
    geom_line(aes(y = Absorbance_blue), color = "blue", size = 0.8, alpha = 0.6) +
    labs(title = paste(ID, "Full Spectrum"), x = "Wavelength (nm)", y = "Absorbance") +
    theme_minimal()
  
  # Inset: dark-only fits (normalized)
  inset <- ggplot(dark_only_data, aes(x = Wavelength, y = Renorm_Absorbance)) +
    geom_line(color = "black", size = 0.8) +
    theme_minimal(base_size = 8) +
    labs(title = "Dark-Only Fit", x = NULL, y = NULL)
  
  colors <- c(
    Gacic = "darkorange",
    Govardovski = "darkgreen",
    Stavenga = "purple",
    Gaussian = "brown"
  )
  
  for (fit in fit_results) {
    inset <- inset + geom_line(aes(y = fit$predicted), color = colors[fit$name], linewidth = 0.8)
    inset <- inset + annotate("text",
                              x = fit$lambda_max + 5,
                              y = 0.8 - 0.1 * match(fit$name, names(fit_results)),
                              label = paste0(fit$name, ":\n", round(fit$lambda_max, 1), " nm"),
                              color = colors[fit$name], hjust = 0)
  }
  
  vp <- viewport(width = 0.45, height = 0.45, x = 0.75, y = 0.75)
  print(main_plot)
  print(inset, vp = vp)
}


# --- Ensure Proper Column Names in process_spectrum() ---
# --- Combined Plot: Main = Difference + Gaussian, Inset = Dark & Light 375–675nm ---
plot_dark_light_with_inset <- function(result, fits_diff) {
  # Main: difference spectrum with Gaussian fit
  colors <- c("Difference spectrum" = "black", "Gaussian model" = "gray")
  
  p_main <- ggplot(result$data, aes(x = Wavelength)) +
    geom_line(aes(y = Renorm_Absorbance, color = "Difference spectrum"), linewidth = 1.0) +
    labs(title = paste(result$ID, "Difference Spectrum"),
         x = "Wavelength (nm)", y = "Normalized Absorbance") +
    coord_cartesian(xlim = c(435, NA)) +
    scale_color_manual(values = colors, breaks = names(colors)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.085, 0.95),
      legend.title = element_blank()#,
      #legend.background = element_rect(color = NULL, fill = "white")
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
  
  # Inset: dark/light spectra from 375 to 675
  inset_colors <- c("Dark" = "black", "Blue" = "blue")
  inset_data <- result$full_spectrum %>%
    dplyr::filter(Wavelength >= 400 & Wavelength <= 650) %>%
    dplyr::mutate(
      Absorbance_dark = Absorbance_dark - min(Absorbance_dark, na.rm = TRUE),
      Absorbance_blue = Absorbance_blue - min(Absorbance_blue, na.rm = TRUE)
    )
  
  p_inset <- ggplot(inset_data, aes(x = Wavelength)) +
    geom_line(aes(y = Absorbance_dark, color = "Dark"), linewidth = 0.7) +
    geom_line(aes(y = Absorbance_blue, color = "Blue"), linewidth = 0.7) +
    labs(x = "Wavelength (nm)", y = "Absorbance") +
    scale_color_manual(values = inset_colors, breaks = names(inset_colors)) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.9, 0.9),
      legend.title = element_blank(),
      #legend.background = element_rect(color = "black", fill = "white")
    )
  
  # Vertical line based on Gaussian λmax
  gaussian_fit <- Filter(function(f) f$name == "Gaussian", fits_diff)
  if (length(gaussian_fit) == 1) {
    x <- gaussian_fit[[1]]$lambda_max
    vline_pos <- 0.978 * x + 7.3285
    
    p_inset <- p_inset +
      geom_vline(xintercept = vline_pos, color = "red", linetype = "dashed") +
      annotate("text", x = vline_pos + 5,
               y = max(inset_data$Absorbance_dark, na.rm = TRUE),
               label = paste0("λmax = ", round(vline_pos, 1), " nm"),
               color = "red", hjust = 0, size = 3)
  }
  
  vp <- viewport(width = 0.50, height = 0.50, x = 0.72, y = 0.68)
  print(p_main)
  print(p_inset, vp = vp)
}

process_spectrum <- function(ID, workpath,
                             symmetry_tolerance = 0.005,
                             hard_min = NULL, hard_max = 750) {
  dark_file <- file.path(workpath, paste0(ID, "dark_average.TXT"))
  blue_file <- file.path(workpath, paste0(ID, "blue_average.TXT"))
  
  darkskip <- detect_header_end(dark_file)
  blueskip <- detect_header_end(blue_file)
  
  darkdata <- suppressMessages(read_table(dark_file, skip = darkskip, col_names = FALSE))
  bluedata <- suppressMessages(read_table(blue_file, skip = blueskip, col_names = FALSE))
  
  darkdata <- darkdata[nrow(darkdata):1, ]; colnames(darkdata) <- c("Wavelength", "Absorbance")
  bluedata <- bluedata[nrow(bluedata):1, ]; colnames(bluedata) <- c("Wavelength", "Absorbance")
  
  df <- inner_join(darkdata, bluedata, by = "Wavelength", suffix = c("_dark", "_blue")) %>%
    mutate(
      darkstrip = Absorbance_dark - min(Absorbance_dark),
      bluestrip = Absorbance_blue - min(Absorbance_blue),
      Renorm_Absorbance = {
        spec <- darkstrip - bluestrip
        spec0 <- spec - min(spec)
        spec0 / max(spec0)
      }
    )
  
  visspec <- df$Renorm_Absorbance
  peak_index <- which.max(visspec)
  peak_value <- visspec[peak_index]
  
  left_index <- peak_index
  while (left_index > 1 && abs(visspec[left_index] - peak_value) > symmetry_tolerance) {
    left_index <- left_index - 1
  }
  right_index <- peak_index
  while (right_index < length(visspec) && abs(visspec[right_index] - peak_value) > symmetry_tolerance) {
    right_index <- right_index + 1
  }
  
  if ((right_index - left_index) < 5) {
    message("⚠️ Detected range too narrow — using 350–750 nm")
    wl_min <- 350
    wl_max <- 750
  } else {
    wl_min <- df$Wavelength[left_index]
    wl_max <- df$Wavelength[right_index]
  }
  
  # Auto-adjust hard_min if NULL using average of 700–750 nm tail
  if (is.null(hard_min)) {
    tail_avg <- {
      df_tail <- dplyr::filter(df, Wavelength >= 700 & Wavelength <= 750)
      mean(df_tail$Renorm_Absorbance, na.rm = TRUE)
    }
    
    left_match <- dplyr::filter(
      df,
      Wavelength < 500 & abs(Renorm_Absorbance - tail_avg) < 0.01
    )
    
    if (nrow(left_match) > 0) {
      hard_min <- max(left_match$Wavelength)
      message("🔧 Auto-adjusted hard_min based on tail average to ", hard_min)
    } else {
      hard_min <- 350
      message("🔧 No tail match found — fallback hard_min to 350")
    }
  }
  
  
  wl_min <- max(hard_min, wl_min)
  
  # Filter and manually construct the renormalized difference spectrum
  df_filtered <- dplyr::filter(df, Wavelength >= wl_min, Wavelength <= wl_max)
  
  wl_mask <- df$Wavelength >= wl_min & df$Wavelength <= wl_max
  spec <- df$darkstrip[wl_mask] - df$bluestrip[wl_mask]
  spec0 <- spec - min(spec)
  renorm <- spec0 / max(spec0)
  
  df_filtered$Renorm_Absorbance <- renorm
  visspec <- renorm
  index <- which.max(visspec)
  peakWave <- if (index > 2 & index < (length(visspec) - 2)) {
    mean(df_filtered$Wavelength[(index - 2):(index + 2)])
  } else {
    df_filtered$Wavelength[index]
  }
  
  
  
  dark_only <- darkdata %>%
    filter(Wavelength >= wl_min, Wavelength <= wl_max) %>%
    mutate(Absorbance = Absorbance - min(Absorbance),
           Renorm_Absorbance = Absorbance / max(Absorbance))
  
  dark_only_index <- which.max(dark_only$Renorm_Absorbance)
  peakWave_dark <- if (dark_only_index > 2 & dark_only_index < (nrow(dark_only) - 2)) {
    mean(dark_only$Wavelength[(dark_only_index - 2):(dark_only_index + 2)])
  } else {
    dark_only$Wavelength[dark_only_index]
  }
  
  full_spectrum <- inner_join(darkdata, bluedata, by = "Wavelength", suffix = c("_dark", "_blue"))
  
  filtered_data <- data.frame(Wavelength = df_filtered$Wavelength,
                              Renorm_Absorbance = visspec)
  dark_only_data <- data.frame(Wavelength = dark_only$Wavelength,
                               Renorm_Absorbance = dark_only$Renorm_Absorbance)
  
  return(list(data = filtered_data,
              dark_only_data = dark_only_data,
              full_spectrum = full_spectrum,
              peakWave = peakWave,
              peakWave_dark = peakWave_dark,
              ID = ID,
              wl_min = wl_min,
              wl_max = wl_max,
              hard_min = hard_min
  ))
}

# --- Summary Table ---
print_fit_summary <- function(fits, peakWave, label = "", hard_min = NULL) {
  if (length(fits) == 0) {
    cat("No fits available for", label, "
")
    return()
  }
  df <- data.frame(
    Model = sapply(fits, function(x) x$name),
    LambdaMax = sapply(fits, function(x) round(x$lambda_max, 2)),
    R_squared = sapply(fits, function(x) round(x$r_squared, 4))
  )
  cat("
Fit summary (", label, ")
", sep = "")
  print(df, row.names = FALSE)
  cat("
Peak wavelength: ", round(peakWave, 2), " nm
")
  if (!is.null(hard_min)) cat("Estimated hard_min: ", round(hard_min, 2), " nm
")
}



# --- Example Run Block ---
result <- process_spectrum(ID = "SAMPLENAME_", workpath = "./")

# --- Fit Templates to Difference Spectrum ---
fits_diff <- list(
  Gacic = fit_template_model(result$data, gacic_template_complex, "Gacic", result$peakWave),
  Govardovski = fit_template_model(result$data, govardovski_template_complex, "Govardovski", result$peakWave),
  Stavenga = fit_template_model(result$data, stavenga_template_complex, "Stavenga", result$peakWave),
  Gaussian = fit_template_model(result$data, function(w, lmax) gaussian_template(w, lmax, width = 40), "Gaussian", result$peakWave)
  )

fits_diff <- Filter(Negate(is.null), fits_diff)
fits_dark_only <- list(
  Gacic = fit_template_model(result$dark_only_data, gacic_template_complex, "Gacic", result$peakWave),
  Govardovski = fit_template_model(result$dark_only_data, govardovski_template_complex, "Govardovski", result$peakWave),
  Stavenga = fit_template_model(result$dark_only_data, stavenga_template_complex, "Stavenga", result$peakWave)
)

fits_dark_only <- Filter(Negate(is.null), fits_dark_only)

plot_dark_only_with_inset(result$full_spectrum, result$dark_only_data, fits_dark_only, result$peakWave, result$ID)

print_fit_summary(fits_dark_only, result$peakWave_dark, label = paste(result$ID, "(dark only)"), hard_min = result$hard_min)
print_fit_summary(fits_diff, result$peakWave, label = paste(result$ID, "(dark - blue diff)"), hard_min = result$hard_min)

plot_dark_light_fit(result$data, fits_diff, result$ID)
plot_dark_light_with_inset(result, fits_diff)
