df.sub <- df[which(df$V2 == 'stopgain'), ]

# Ensure factor ordering for categories
df.sub$last.EJC <- factor(df.sub$last.EJC, levels = c("upstream", "penultimate.last50bp", "last.exon"))
df.sub <- df.sub[order(df.sub$last.EJC), ]

# Calculate sample size for each category
sample_size_data <- df.sub %>%
  group_by(last.EJC) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot with box and sample size annotations
bp <- df.sub %>% ggplot(aes(x = last.EJC, y = ALLELE.RAT, fill = last.EJC)) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +  # Boxplot only
  coord_flip() +  # Flip the axes for horizontal plot
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.background = element_blank()
  ) +
  ggtitle("Canonical Rule") +
  xlab("") +
  ylab("NMD Efficiency") +
  ylim(0, 1.1) +
  geom_signif(
    comparisons = list(c('penultimate.last50bp', 'last.exon'),
                       c('upstream', 'last.exon'),
                       c('penultimate.last50bp', 'upstream')),
    map_signif_level = TRUE,
    stat = "signif",
    position = "identity",
    test = "wilcox.test",
    textsize = 4,
    y_position = c(0.65, 0.8, 0.9)
  ) +
  scale_fill_manual(values = c("upstream" = "red",  # Color for 'upstream'
                               "penultimate.last50bp" = "darkgrey",  # Corrected color for 'penultimate.last50bp'
                               "last.exon" = "darkgrey")) +  # Color for 'last.exon'
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5) +  # Cutoff line at 0.5
  # Add sample size as text labels above each boxplot
  geom_text(data = sample_size_data, aes(x = last.EJC, y = 1.03, label = count),
            size = 4, position = position_dodge(width = 0.8), vjust = -0.5)

bp

