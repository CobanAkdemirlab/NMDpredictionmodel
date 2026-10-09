# ══════════════════════════════════════════════════════════════════════════════
# LMNA saturation genome editing (Cortázar et al. 2025) vs TrunKitten (8-feature)
#   Experimental NMD efficiency (norm_mean) by TrunKitten binary call:
#   NMD-triggering vs NMD-escape (escape_pred_at_youden, threshold_used in data)
# Input : SGE_LMNA_comparing_8trunkitten.Rdata  (object LMNA_validation_8f_84)
# Output: LMNA_TrunKitten_binary.pdf / .png
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggbeeswarm)
})
load("SGE_LMNA_comparing_8trunkitten.Rdata")
d <- as.data.frame(LMNA_validation_8f_84)

TAU <- unique(d$threshold_used)
d$call <- factor(ifelse(d$escape_pred_at_youden == 1, "NMD-escape", "NMD-triggering"),
                 levels = c("NMD-triggering", "NMD-escape"))
d$exon <- factor(d$mut.exon)

wt  <- wilcox.test(norm_mean ~ call, data = d)
nlab <- aggregate(norm_mean ~ call, d, length)
cat(sprintf("Threshold = %.3f\n", TAU)); print(table(d$call, d$exon))
print(aggregate(norm_mean ~ call, d, median))
cat(sprintf("Wilcoxon rank-sum P = %.2g\n", wt$p.value))

p_lab <- if (wt$p.value < 1e-4) {
  sprintf("italic(P) == %s %%*%% 10^%d",
          formatC(wt$p.value / 10^floor(log10(wt$p.value)), format = "f", digits = 1),
          floor(log10(wt$p.value)))
} else sprintf("italic(P) == %.3f", wt$p.value)

p <- ggplot(d, aes(call, norm_mean)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = "grey95",
               colour = "grey20", linewidth = 0.6) +
  geom_quasirandom(aes(fill = exon), shape = 21, size = 3, colour = "grey15",
                   stroke = 0.4, width = 0.18) +
  annotate("segment", x = 1, xend = 2, y = 1.10, yend = 1.10, linewidth = 0.5) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = 1.10, yend = 1.07, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.16, label = p_lab, parse = TRUE, size = 4.5) +
  geom_text(data = nlab, aes(x = call, y = -0.17, label = paste0("n = ", norm_mean)),
            size = 4.3) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.2, 1.2)) +
  labs(x = sprintf("TrunKitten prediction (threshold = %.2f)", TAU),
       y = "Experimental NMD efficiency", fill = "LMNA exon") +
  theme_classic(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "right")

ggsave("LMNA_TrunKitten_binary.pdf", p, width = 6.5, height = 4.6, device = cairo_pdf)
ggsave("LMNA_TrunKitten_binary.png", p, width = 6.5, height = 4.6, dpi = 300)
cat("Saved: LMNA_TrunKitten_binary.pdf / .png\n")
