# Minimal path diagrams for piecewise SEMs — combined (a)/(b) panel with R² and indirect effects
library(DiagrammeR)
library(DiagrammeRsvg)  # for export
library(rsvg)           # for export

# Helper to build the combined Graphviz source
sem_grviz_pair <- function(
    # Panel (a): Richness (standardized betas)
  beta_div_syn_a, se_div_syn_a, sig_div_syn_a,
  beta_syn_stab_a, se_syn_stab_a, sig_syn_stab_a,
  beta_div_stab_a, se_div_stab_a, sig_div_stab_a,
  R2m_sync_a, R2c_sync_a,
  R2m_stab_a, R2c_stab_a,
  
  # Panel (b): Shannon (standardized betas)
  beta_div_syn_b, se_div_syn_b, sig_div_syn_b,
  beta_syn_stab_b, se_syn_stab_b, sig_syn_stab_b,
  beta_div_stab_b, se_div_stab_b, sig_div_stab_b,
  R2m_sync_b, R2c_sync_b,
  R2m_stab_b, R2c_stab_b,
  
  # Panel titles
  title_a = "(a) Richness (q = 0)", 
  title_b = "(b) Shannon (q = 1)"
){
  fmt   <- function(x, d=3) sprintf(paste0("%.", d, "f"), x)
  lab   <- function(b, se, sig) paste0(fmt(b), " \u00B1 ", fmt(se), if (isTRUE(sig)) "" else " (ns)")
  lt    <- function(sig) if (isTRUE(sig)) "solid" else "dashed"
  col   <- function(b) ifelse(b < 0, "red", "black")
  
  # Indirect effects (standardized)
  ind_a <- beta_div_syn_a * beta_syn_stab_a
  ind_b <- beta_div_syn_b * beta_syn_stab_b
  
  g <- paste0(
    'digraph SEM {
  graph [rankdir = LR, fontsize = 10]
  node  [shape = box, style = rounded, fontname = Helvetica, fontsize = 11]
  edge  [penwidth = 2]

  # ---------------- Panel (a) ----------------
  subgraph cluster_a {
    label="', title_a, '  |  Indirect (std) = ', fmt(ind_a), '";
    labelloc="t"; labeljust="l"; fontsize=12; fontname=Helvetica;
    style="rounded"; color="gray70";

    A1 [label="Richness (q = 0)"]
    B1 [label="Synchrony (log \u03C6)\nR²m = ', fmt(R2m_sync_a, 3), ', R²c = ', fmt(R2c_sync_a, 3), '"]
    C1 [label="Stability (log10 I_C)\nR²m = ', fmt(R2m_stab_a, 3), ', R²c = ', fmt(R2c_stab_a, 3), '"]

    A1 -> B1 [label="', lab(beta_div_syn_a, se_div_syn_a, sig_div_syn_a), '", color="', col(beta_div_syn_a), '", style="', lt(sig_div_syn_a), '"]
    B1 -> C1 [label="', lab(beta_syn_stab_a, se_syn_stab_a, sig_syn_stab_a), '", color="', col(beta_syn_stab_a), '", style="', lt(sig_syn_stab_a), '"]
    A1 -> C1 [label="', lab(beta_div_stab_a, se_div_stab_a, sig_div_stab_a), '", color="', col(beta_div_stab_a), '", style="', lt(sig_div_stab_a), '"]
  }

  # ---------------- Panel (b) ----------------
  subgraph cluster_b {
    label="', title_b, '  |  Indirect (std) = ', fmt(ind_b), '";
    labelloc="t"; labeljust="l"; fontsize=12; fontname=Helvetica;
    style="rounded"; color="gray70";

    A2 [label="Shannon (q = 1)"]
    B2 [label="Synchrony (log \u03C6)\nR²m = ', fmt(R2m_sync_b, 3), ', R²c = ', fmt(R2c_sync_b, 3), '"]
    C2 [label="Stability (log10 I_C)\nR²m = ', fmt(R2m_stab_b, 3), ', R²c = ', fmt(R2c_stab_b, 3), '"]

    A2 -> B2 [label="', lab(beta_div_syn_b, se_div_syn_b, sig_div_syn_b), '", color="', col(beta_div_syn_b), '", style="', lt(sig_div_syn_b), '"]
    B2 -> C2 [label="', lab(beta_syn_stab_b, se_syn_stab_b, sig_syn_stab_b), '", color="', col(beta_syn_stab_b), '", style="', lt(sig_syn_stab_b), '"]
    A2 -> C2 [label="', lab(beta_div_stab_b, se_div_stab_b, sig_div_stab_b), '", color="', col(beta_div_stab_b), '", style="', lt(sig_div_stab_b), '"]
  }

  # Notes
  note [label="Random intercepts: bay, season", shape=note, fontsize=10]
  note -> C1 [arrowhead=none, color="gray60"]
  note -> C2 [arrowhead=none, color="gray60"]
}'
  )
  
  DiagrammeR::grViz(g)
}

# ----- Fill with your standardized path coefficients and R² values -----
# From your outputs:
# Richness model (standardized):
#   richness -> log_phi = -0.1781 (p = 0.0083, sig)
#   log_phi  -> log_IC  = -0.9105 (p < 0.001, sig)
#   richness -> log_IC  =  0.0329 (p = 0.3068, ns)
#   R²: log_phi (m=0.0301, c=0.3213); log_IC (m=0.8018, c=0.8434)

# Shannon model (standardized):
#   shannon -> log_phi = -0.7010 (p < 0.001, sig)
#   log_phi -> log_IC  = -0.8510 (p < 0.001, sig)
#   shannon -> log_IC  =  0.1168 (p = 0.0014, sig)
#   R²: log_phi (m=0.4812, c=0.5482); log_IC (m=0.8033, c=0.8575)

# g_combined <- sem_grviz_pair(
#   # (a) Richness panel
#   beta_div_syn_a = -0.1781, sig_div_syn_a = TRUE,
#   beta_syn_stab_a= -0.9105, sig_syn_stab_a= TRUE,
#   beta_div_stab_a=  0.0329, sig_div_stab_a= FALSE,
#   R2m_sync_a = 0.03012467, R2c_sync_a = 0.3213154,
#   R2m_stab_a = 0.80180065, R2c_stab_a = 0.8434445,
#   # (b) Shannon panel
#   beta_div_syn_b = -0.7010, sig_div_syn_b = TRUE,
#   beta_syn_stab_b= -0.8510, sig_syn_stab_b = TRUE,
#   beta_div_stab_b=  0.1168, sig_div_stab_b = TRUE,
#   R2m_sync_b = 0.4811806, R2c_sync_b = 0.5481526,
#   R2m_stab_b = 0.8032674, R2c_stab_b = 0.8575214,
#   # Titles (kept default)
#   title_a = "(a) Richness (q = 0)",
#   title_b = "(b) Shannon (q = 1)"
# )

g_combined <- sem_grviz_pair(
  # (a) Richness panel
  beta_div_syn_a = -0.1781, se_div_syn_a = 0.067, sig_div_syn_a = TRUE,
  beta_syn_stab_a= -0.9105, se_syn_stab_a = 0.026, sig_syn_stab_a= TRUE,
  beta_div_stab_a=  0.0329, se_div_stab_a = 0.029, sig_div_stab_a= FALSE,
  R2m_sync_a = 0.03012467,  R2c_sync_a = 0.3213154,
  R2m_stab_a = 0.80180065,  R2c_stab_a = 0.8434445,
  
  # (b) Shannon panel
  beta_div_syn_b = -0.7010, se_div_syn_b = 0.045, sig_div_syn_b = TRUE,
  beta_syn_stab_b= -0.8510, se_syn_stab_b = 0.032, sig_syn_stab_b = TRUE,
  beta_div_stab_b=  0.1168, se_div_stab_b = 0.036, sig_div_stab_b = TRUE,
  R2m_sync_b = 0.4811806,   R2c_sync_b = 0.5481526,
  R2m_stab_b = 0.8032674,   R2c_stab_b = 0.8575214,
  
  title_a = "(a) Richness (q = 0)", title_b = "(b) Shannon (q = 1)"
)

# View in R
g_combined

svg_code <- export_svg(g_combined)
writeLines(svg_code, "sem_combined_ab.svg")

# Fallback: render PNG via {rsvg}, then convert to TIFF via {magick}
# install.packages(c("rsvg","magick"))  # if needed
library(rsvg); library(magick)

rsvg_png("Fig_4_sem_combined_ab.svg", "Fig_4_sem_combined_ab.png", width = 1800, height = 700)
image_write(image_read("Fig_4_sem_combined_ab.png"), path = "Fig_4_sem_combined_ab.tif",
            format = "tiff", compression = "lzw")

# Suggested caption (plain text):
# Minimal path diagrams for the piecewise SEMs. Values on arrows are standardized
# path coefficients (STD β); dashed lines indicate non-significant paths. Negative
# coefficients from diversity to synchrony indicate that greater diversity is
# associated with lower synchrony (i.e., more asynchrony). Indirect effect (std)
# is the product of diversity → synchrony and synchrony → stability. Models include
# random intercepts for bay and season; node labels show R² (marginal, conditional).
