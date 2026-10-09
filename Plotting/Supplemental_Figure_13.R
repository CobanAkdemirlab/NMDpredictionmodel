# ══════════════════════════════════════════════════════════════════════════════
# Supplemental Figure 13
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

####
# Supplemental Figure S12 — LMNA SGE (Cortázar et al. 2025) vs 8-feature TrunKitten
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggbeeswarm); library(patchwork)})
load("SGE_LMNA_comparing_8trunkitten.Rdata")
d <- as.data.frame(LMNA_validation_8f_84)
d$exon <- factor(d$mut.exon)
TAU <- unique(d$threshold_used)
d$call <- factor(ifelse(d$escape_pred_at_youden == 1, "NMD-escape", "NMD-triggering"),
                 levels = c("NMD-triggering", "NMD-escape"))
ct <- cor.test(d$predicted_NMD, d$norm_mean, method = "spearman", exact = FALSE)
wt <- wilcox.test(norm_mean ~ call, data = d)
plab <- function(p) { e <- floor(log10(p)); sprintf("italic(P) == %s %%*%% 10^%d", formatC(p / 10^e, format = "f", digits = 1), e) }
th <- theme_classic(base_size = 15) + theme(axis.text = element_text(colour = "black"),
                                            plot.tag = element_text(size = 22, face = "bold"))
pA <- ggplot(d, aes(predicted_NMD, norm_mean)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey75", linewidth = 1) +
  geom_point(aes(fill = exon), shape = 21, size = 3, colour = "grey15", stroke = 0.4) +
  annotate("text", x = 0.02, y = 1.12, hjust = 0, size = 4.5, parse = TRUE,
           label = sprintf("rho == '%.3f'*','~~%s", ct$estimate, plab(ct$p.value))) +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.75, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.2, 1.2)) +
  labs(x = "TrunKitten-predicted NMD probability", y = "Experimental NMD efficiency",
       fill = "LMNA exon", tag = "A") + th
nlab <- aggregate(norm_mean ~ call, d, length)
pB <- ggplot(d, aes(call, norm_mean)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = "grey95", colour = "grey20", linewidth = 0.6) +
  geom_quasirandom(aes(fill = exon), shape = 21, size = 3, colour = "grey15", stroke = 0.4, width = 0.18) +
  annotate("segment", x = 1, xend = 2, y = 1.08, yend = 1.08, linewidth = 0.5) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = 1.08, yend = 1.05, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.14, label = plab(wt$p.value), parse = TRUE, size = 4.5) +
  geom_text(data = nlab, aes(x = call, y = -0.17, label = paste0("n = ", norm_mean)), size = 4.3) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.2, 1.2)) +
  labs(x = sprintf("TrunKitten prediction (threshold = %.2f)", TAU), y = "Experimental NMD efficiency",
       fill = "LMNA exon", tag = "B") + th
fig <- (pA | pB) + plot_layout(guides = "collect", widths = c(1.15, 1))
ggsave("SupplementalFigure12_LMNA.pdf", fig, width = 12, height = 4.8, device = cairo_pdf)
ggsave("SupplementalFigure12_LMNA.png", fig, width = 12, height = 4.8, dpi = 300)
cat(sprintf("rho = %.3f, P = %.2g; Wilcoxon P = %.2g\n", ct$estimate, ct$p.value, wt$p.value))

