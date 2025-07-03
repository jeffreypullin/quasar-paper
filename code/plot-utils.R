
to_r_vec <- function(str) {
    str <- str_replace(str, "\\[", "\\(")
    str <- str_replace(str, "\\]", "\\)")
    str <- str_replace_all(str, "([^,\\(\\)]+)", "'\\1'")
    str <- paste0("c", str)
    out <- eval(parse(text = str))
    out <- str_trim(out)
    out
}

acat <- function(pvals) {
  n <- length(pvals)
  sum <- sum(qcauchy(pvals) / n) 
  pcauchy(sum) 
}

theme_jp <- function() {

  font <- "Helvetica"

  theme(
    plot.title = element_text(
      family =  "Helvetica", size = 24, color = "#222222"
    ),
    plot.subtitle = element_text(
      family = font, size = 20, margin = ggplot2::margin(0, 0, 10, 0)
    ),
    plot.caption = element_blank(),
    plot.title.position = "plot",
    legend.position = "top",
    legend.box.margin = margin(t = -5),
    legend.text.align = 0,
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(family = font, size = 14, color = "#222222"),
    axis.title = element_text(family = font, size = 14, color = "#222222"),
    axis.text = element_text(family = font, size = 14, color = "#222222"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major.y = element_line(color = "#cbcbcb"),
    panel.grid.minor.y = element_line(color = "#cdcdcd"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(family = font, size = 18)
  )
}

theme_jp_vgrid <- function() {
  theme_jp() %+replace%
    theme(
      panel.grid.major.x = element_line(color = "#cbcbcb"),
      panel.grid.minor.x = element_line(color = "#cdcdcd"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
    )
}

method_lookup <- c(
  "tensorqtl" = "tensorQTL (LM)",
  "tensorqtl_acat" = "tensorQTL (ACAT)",
  "tensorqtl_qvalue" = "tensorQTL (qvalue)",
  "tensorqtl_acat_qvalue" = "tensorQTL (ACAT, qvalue)",
  "apex" = "apex (LMM)",
  "jaxqtl" = "jaxQTL (NB-GLM)",
  "jaxqtl_cis" = "jaxQTL (cis)",
  "jaxqtl_cis_nominal" = "jaxQTL (cis nominal)",
  "tensorqtl_cis" = "tensorQTL (cis)",
  "tensorqtl_cis_nominal" = "tensorQTL (cis nominal)",
  "quasar-nb_glm" = "quasar (NB-GLM)",
  "quasar-nb_glm-apl" = "quasar (NB-GLM, APL)",
  "quasar-p_glm" = "quasar (P-GLM)",
  "quasar-lm" = "quasar (LM)",
  "quasar-p_glmm" = "quasar (P-GLMM)",
  "quasar-nb_glmm" = "quasar (NB-GLMM)",
  "quasar-lmm" = "quasar (LMM)"
)

model_lookup <- c(
  "nb_glm" = "NB GLM",
  "nb_glm-apl" = "NB GLM, APL",
  "p_glm" = "P GLM",
  "lm" = "LM",
  "lmm" = "LMM",
  "p_glmm" = "P GLMM",
  "nb_glmm" = "NB GLMM"
)

cell_type_cols <- c(
  "Plasma" = "#66CCEE",
  "B IN" = "#228833",
  "CD4 NC" = "#EE6677"
)

pvalue_method_lookup <- c(
  "pvalue_perm_q" = "Perm, q-value",
  "pvalue_acat_bh" = "ACAT, BH",
  "pvalue_acat_q" = "ACAT, q-value",
  "pvalue_perm_bh" = "Perm, BH"
)
