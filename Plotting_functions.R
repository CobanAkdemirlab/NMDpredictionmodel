# Load required packages
library(ggplot2)
library(dplyr)
library(ggsignif)
library (mltools)
library(ggh4x)

plot_NMD_efficiency <- function(df, category_column, y_var = "ALLELE.RAT",
                                title = "NMD Efficiency Plot", 
                                colors = NULL,
                                signif_comparisons = NULL,
                                y_positions = NULL,
                                flip_axis = TRUE) {
  
  # Ensure category column is a factor with correct ordering
  df[[category_column]] <- factor(df[[category_column]], 
                                  levels = unique(df[[category_column]]))
  df <- df[order(df[[category_column]]), ]
  
  # Compute sample sizes for labeling
  sample_size_data <- df %>%
    group_by(.data[[category_column]]) %>%
    summarise(count = n(), .groups = "drop")
  
  # Default colors if not provided
  if (is.null(colors)) {
    colors <- c("upstream" = "red", "penultimate.last50bp" = "darkgrey", "last.exon" = "darkgrey")
  }
  
  # Create the plot
  bp <- ggplot(df, aes(x = .data[[category_column]], y = .data[[y_var]], fill = .data[[category_column]])) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    theme_minimal() +
    labs(
      title = title,
      x = "",
      y = "NMD Efficiency",
      fill = category_column
    ) +
    scale_fill_manual(values = colors) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      panel.background = element_blank()
    ) +
    ylim(0, 1.3) +  # Adjust for sample size labels
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5) + # Cutoff line
  
    # Add sample size labels
    geom_text(data = sample_size_data, aes(x = .data[[category_column]], y = 1.1, label = count),
              size = 4, position = position_dodge(width = 0.75), vjust = -0.5)
  
  # Flip coordinates if enabled
  if (flip_axis) {
    bp <- bp + coord_flip()
  }
  
  # Add significance comparisons if provided
  if (!is.null(signif_comparisons) && !is.null(y_positions)) {
    bp <- bp + geom_signif(
      comparisons = signif_comparisons,
      map_signif_level = TRUE,
      textsize = 4,
      y_position = y_positions
    )
  }
  
  return(bp)
}
