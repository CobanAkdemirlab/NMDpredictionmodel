
plot_NMD_efficiency_multi <- function(df, feature, feature_label, bin_breaks, bin_labels, 
                                      facet_var = NULL, reference_bin = NULL, color_map = NULL) {
  df[[feature]] <- as.numeric(df[[feature]])
  df$feature_binned <- cut(df[[feature]], breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE)
  df <- df[!is.na(df$feature_binned), ]
  df$feature_binned <- factor(df$feature_binned, levels = bin_labels)
  
  if (!is.null(facet_var)) df[[facet_var]] <- factor(df[[facet_var]])

  count_data <- df %>%
    group_by(across(all_of(c("feature_binned", facet_var)))) %>%
    summarise(count = n(), .groups = "drop")
  
  if (is.null(color_map)) {
    color_map <- rep("darkgrey", length(bin_labels))
    names(color_map) <- bin_labels
    color_map[length(bin_labels)] <- "red"
  }

  bp <- ggplot(df, aes(x = feature_binned, y = ALLELE.RAT, fill = feature_binned)) +
    geom_boxplot(width = 0.2, color = "black", alpha = 1) +
    scale_fill_manual(values = color_map) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16, face = "bold"),
      text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(paste("NMD Efficiency vs.", feature_label)) +
    xlab(feature_label) +
    ylab("NMD Efficiency") +
    ylim(0, 1.3) +  
    coord_flip() +
    geom_text(data = count_data, aes(x = feature_binned, y = 1.1, label = count),
              size = 4, position = position_dodge(width = 0.75), vjust = -0.5)

  if (!is.null(facet_var)) bp <- bp + facet_wrap(as.formula(paste("~", facet_var)))

  if (!is.null(reference_bin)) {
    labels <- levels(df$feature_binned)[levels(df$feature_binned) != reference_bin]
    y_positions <- seq(0.3, 0.9, length.out = length(labels))

    bp <- bp + geom_signif(
      comparisons = lapply(labels, function(label) c(reference_bin, label)),
      map_signif_level = TRUE,
      textsize = 3.5,
      y_position = y_positions
    )
  }

  return(bp)
}
