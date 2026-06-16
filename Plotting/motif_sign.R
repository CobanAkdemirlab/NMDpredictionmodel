valid_regions <- c("utr3_all","utr3_200","newutr_all","newutr_200",
                   "ptc_pm100","ejc_pm100","ptc_to_ejc")

motif_col_idx <- 158:1032

y_candidates  <- c("ALLELE.RAT","allele.rat","allele_rat",
                   "NMD_efficiency","nmd_efficiency","nmd.efficiency")

out_csv <- "TOPMed_motif_region_effect_with_CI_autoApril152026.csv"

recode_presence <- function(x){
  x <- as.character(x)
  x[grepl("^yes", x, ignore.case=TRUE)] <- "Yes"
  x[grepl("^no",  x, ignore.case=TRUE)] <- "No"
  x[x %in% c("1","TRUE","T","true","t")] <- "Yes"
  x[x %in% c("0","FALSE","F","false","f")] <- "No"
  x[!(x %in% c("Yes","No"))] <- NA
  factor(x, levels=c("No","Yes"))
}

get_region_from_col <- function(colname) sub("\\.[^.]+$", "", colname)
get_motif_from_col  <- function(colname) sub("^.*\\.",   "", colname)
compute_ci_from_raw <- function(df_raw){
  
  ## find outcome column
  y_col <- NA_character_
  for (nm in names(df_raw)) {
    if (tolower(nm) %in% tolower(y_candidates)) {
      y_col <- nm
      break
    }
  }
  
  if (is.na(y_col)) {
    stop("Could not find outcome column like ALLELE.RAT")
  }
  
  message("Using outcome column: ", y_col)
  
  ## motif columns
  motif_cols <- names(df_raw)[motif_col_idx]
  
  y <- df_raw[[y_col]]
  
  out_list <- list()
  k <- 1
  
  for (cn in motif_cols) {
    
    region <- get_region_from_col(cn)
    motif  <- get_motif_from_col(cn)
    
    if (!(region %in% valid_regions)) next
    
    grp <- recode_presence(df_raw[[cn]])
    
    keep <- !is.na(grp) & !is.na(y)
    grp  <- grp[keep]
    yy   <- y[keep]
    
    n1 <- sum(grp == "Yes")
    n0 <- sum(grp == "No")
    
    if (n1 < 2 || n0 < 2) next
    
    y1 <- yy[grp == "Yes"]
    y0 <- yy[grp == "No"]
    
    m1 <- mean(y1)
    m0 <- mean(y0)
    
    s1 <- sd(y1)
    s0 <- sd(y0)
    
    eff <- m1 - m0
    se  <- sqrt(s1^2/n1 + s0^2/n0)
    
    ciL <- eff - 1.96*se
    ciH <- eff + 1.96*se
    
    w.res <- wilcox.test(y1, y0)
    
    out_list[[k]] <- data.frame(
      motif=motif,
      region=region,
      feature_col=cn,
      n_present=n1,
      n_absent=n0,
      effect=eff,
      se_diff=se,
      ci_lo=ciL,
      ci_hi=ciH,
      p_wilcox=w.res$p.value
    )
    
    k <- k + 1
  }
  
  res <- do.call(rbind, out_list)
  
  res$q_value <- p.adjust(res$p_wilcox, method="BH")
  
  res <- res[order(res$q_value), ]
  
  return(res)
}
res <- compute_ci_from_raw(df)
sig <- subset(res, q_value < 0.05)
write.csv(res, out_csv, row.names = FALSE)
write.csv(sig, "TOPMed_significant_motifs_April2026.csv", row.names = FALSE)
