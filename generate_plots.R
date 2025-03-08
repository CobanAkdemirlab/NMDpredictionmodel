# List of features with their binning details
features_list <- list(
  list(feature = "canonical_rule", label = "Canonical Rule", breaks = c(0, 1), labels = c("Yes", "No"), reference = "Yes"),
  list(feature = "first_200bp", label = "First 200bp", breaks = c(0, 1), labels = c("No", "Yes"), reference = "No"),
  list(feature = "PTC_2_start", label = "PTC.2 Start", breaks = c(0, 100, 200, 300, 400), labels = c("0-100", "100-200", "200-300", "300-400"), reference = "0-100"),
  list(feature = "UTR3_length", label = "3'UTR Length", breaks = c(0, 10, 14), labels = c("0-10", "10-14"), reference = "0-10"),
  list(feature = "UTR5_uORF", label = "5'UTR uORF", breaks = c(0, 10, 14), labels = c("0-10", "10-14"), reference = "0-10"),
  list(feature = "coding_seq_length", label = "Coding Sequence Length", breaks = c(0, 100, 200, 300, 400, Inf), labels = c("0-100", "100-200", "200-300", "300-400", ">400"), reference = "0-100"),
  list(feature = "PTC_exon_length", label = "PTC Exon Length", breaks = c(2, 46, 67, 84, 97, 126, 158, 200, 316), labels = c("2", "46", "67", "84", "97", "126", "158", "200", "316"), reference = "2"),
  list(feature = "relative_PTC_location", label = "Relative PTC Location", breaks = c(0, 0.31, 0.6, 1), labels = c("0", "0.31", "0.6"), reference = "0"),
  list(feature = "distance_PTC_exon_end", label = "Distance from PTC to Exon End", breaks = c(0, 100, 200, 300, 400, 500, 600, 700), labels = c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", "600-700"), reference = "0-100")
)

# Facet variables to iterate over
facet_vars <- c("last.EJC", "pLI.cat", "Freq.cat")

# Loop through each feature and create the plots
for (feature_info in features_list) {
  for (facet_var in facet_vars) {
    # Skip facet_vars for canonical_rule (only one plot needed)
    if (feature_info$feature == "canonical_rule" && facet_var != "last.EJC") next
    
    # Generate plot
    bp <- plot_NMD_efficiency(
      df = df.sub,
      feature = feature_info$feature,
      feature_label = feature_info$label,
      bin_breaks = feature_info$breaks,
      bin_labels = feature_info$labels,
      facet_var = facet_var,
      reference_bin = feature_info$reference
    )
    
    # Save plot as image
    plot_filename <- paste0("plots/", feature_info$feature, "_", facet_var, ".png")
    ggsave(plot_filename, bp, width = 8, height = 6)
  }
}
