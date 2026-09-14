plot_zscore_points <- function(normal_connections, tnbc_connections, 
                               gene_name = "Gene",
                               colors = c("Normal" = "#008080", "TNBC" = "darkred"),
                               point_alpha = 0.6,
                               point_size = 1.5) {
  
  # Add phenotype column to each dataset
  normal_connections$phenotype <- "Normal"
  tnbc_connections$phenotype <- "TNBC"
  
  # Combine the datasets
  combined_data <- rbind(normal_connections, tnbc_connections)
  
  # Perform Mann-Whitney U test
  test_result <- wilcox.test(normal_connections$zscore, tnbc_connections$zscore)
  p_value <- test_result$p.value
  
  # Format p-value for display
  if (p_value < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- paste0("p = ", round(p_value, 3))
  }
  
  # Create the plot
  p <- ggplot(combined_data, aes(x = phenotype, y = zscore, color = phenotype)) +
    geom_jitter(width = 0.2, alpha = point_alpha, size = point_size) +
    labs(title = paste0("Z-score Distribution - ", gene_name),
         subtitle = paste0("Mann-Whitney U test: ", p_label),
         x = "",
         y = "Z-score") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14)) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5)
  
  # Add significance bracket if significant
  if (p_value < 0.05) {
    y_max <- max(combined_data$zscore) * 1.1
    p <- p + 
      annotate("segment", x = 0.8, xend = 2.2, y = y_max, yend = y_max) +
      annotate("segment", x = 0.8, xend = 0.8, y = y_max * 0.98, yend = y_max) +
      annotate("segment", x = 2.2, xend = 2.2, y = y_max * 0.98, yend = y_max) +
      annotate("text", x = 1.5, y = y_max * 1.02, 
               label = ifelse(p_value < 0.001, "***", 
                              ifelse(p_value < 0.01, "**", "*")),
               size = 5)
  }
  
  # Print summary statistics for reference
  cat("\nSummary for", gene_name, ":\n")
  cat("Normal - n:", nrow(normal_connections), 
      ", median:", round(median(normal_connections$zscore), 4), "\n")
  cat("TNBC - n:", nrow(tnbc_connections), 
      ", median:", round(median(tnbc_connections$zscore), 4), "\n")
  cat("Mann-Whitney U test p-value:", p_value, "\n")
  
  return(p)
}
