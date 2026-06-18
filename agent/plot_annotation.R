#!/usr/bin/env Rscript
# Annotation map for mu1046085 / Zm00004b034131 in the gggenomes style of
# sawers-rellan-labs/LH244_CRISPR (blast_viz_oligos.Rmd):
#   geom_seq backbone, exons as goldenrod feats (y-0.15), amplicons as colored
#   feat bars stacked by yoff (y+yoff), primers as geom_gene arrows
#   (forward=#40A040 green, reverse=#9060C0 purple), bottom legend, dashed gene vlines.
suppressPackageStartupMessages({library(gggenomes); library(ggplot2); library(readr); library(dplyr)})

ROOT <- "/Users/fvrodriguez/Desktop/primermu"
WIN_START <- 150176549L
TSD_G <- 150183971L; TSD_LEN <- 9L
loc <- function(g) g - WIN_START + 1L                 # genomic -> window-local 1-based
LOCUS_LEN <- 9659L
seq_id <- "Zm00004b034131_locus"

gene  <- c(150178549L, 150184207L)                    # minus strand
exons <- list(c(150183805,150184207), c(150181691,150181980),
              c(150179879,150180169), c(150178549,150179052))

## ---- sequence backbone ----
seqs <- data.frame(seq_id = seq_id, length = LOCUS_LEN)

## ---- exons ----
exon_df <- do.call(rbind, lapply(exons, function(e)
  data.frame(seq_id, start = loc(e[1]), end = loc(e[2]), feat_type = "Exon")))

## ---- insertion (TSD) ----
ins_df <- data.frame(seq_id, start = loc(TSD_G), end = loc(TSD_G + TSD_LEN - 1),
                     feat_type = "Insertion", name = "mu1046085")

## ---- primers (arrows) ----
g <- read_tsv(file.path(ROOT, "results/gsp_table.tsv"), show_col_types = FALSE)
g$len <- as.integer(g$len)
prim <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
  r <- g[i, ]
  if (startsWith(r$side, "forward")) {        # gene-fwd: 5' high coord, points to lower (left)
    s5 <- r$`5p_chr9`; s3 <- r$`5p_chr9` - (r$len - 1); ft <- "forward_primer"
  } else {                                    # gene-rev: 5' low coord, points to higher (right)
    s5 <- r$`5p_chr9`; s3 <- r$`5p_chr9` + (r$len - 1); ft <- "reverse_primer"
  }
  data.frame(seq_id, start = loc(s5), end = loc(s3), feat_type = ft, name = r$name)
}))
pos <- function(nm) g$`5p_chr9`[g$name == nm]

## ---- amplicons (stacked bars above backbone) ----
mkamp <- function(name, gstart, gend, ft)
  data.frame(seq_id, start = loc(min(gstart,gend)), end = loc(max(gstart,gend)),
             feat_type = ft, name = name)
amp <- rbind(
  # 3 recommended WT control pairs (F+R), span R5'..F5'
  mkamp("WT F3+R2 377", pos("GSP-R2"), pos("GSP-F3"), "WT_amplicon"),
  mkamp("WT F2+R1 446", pos("GSP-R1"), pos("GSP-F2"), "WT_amplicon"),
  mkamp("WT F1+R3 410", pos("GSP-R3"), pos("GSP-F1"), "WT_amplicon"),
  # 5' junction amplicons (GSP-F + TIR6): genomic footprint = F5' -> insertion
  mkamp("F1+TIR6 ~260", pos("GSP-F1"), TSD_G, "junction_5"),
  mkamp("F2+TIR6 ~237", pos("GSP-F2"), TSD_G, "junction_5"),
  mkamp("F3+TIR6 ~199", pos("GSP-F3"), TSD_G, "junction_5"),
  # 3' junction amplicons (GSP-R + TIR6): genomic footprint = R5' -> insertion
  mkamp("R1+TIR6 ~282", pos("GSP-R1"), TSD_G + TSD_LEN - 1, "junction_3"),
  mkamp("R2+TIR6 ~251", pos("GSP-R2"), TSD_G + TSD_LEN - 1, "junction_3"),
  mkamp("R3+TIR6 ~223", pos("GSP-R3"), TSD_G + TSD_LEN - 1, "junction_3")
)
amp$yoff <- 0.18 + 0.12 * (seq_len(nrow(amp)) - 1)     # stack

## ---- build plot ----
p <- gggenomes(seqs = seqs, genes = prim,
               feats = list(amplicons = amp, exons = exon_df, insertion = ins_df)) +
  geom_seq() +
  geom_feat(data = feats(exons), aes(y = y - 0.15, yend = y - 0.15, color = feat_type),
            position = "identity", linewidth = 5, alpha = 0.85) +
  geom_feat(data = feats(insertion), aes(color = feat_type), linewidth = 7) +
  geom_feat_label(data = feats(insertion), aes(label = name), size = 2.6,
                  nudge_y = -0.12, color = "#D94040") +
  geom_feat(data = feats(amplicons), aes(y = y + yoff, yend = y + yoff, color = feat_type),
            position = "identity", linewidth = 2.6, alpha = 0.75) +
  geom_text(data = amp, aes(x = pmin(start, end), y = 1 + yoff + 0.045, label = name),
            hjust = 0, vjust = 0, size = 2.2, color = "grey20", inherit.aes = FALSE) +
  geom_gene(aes(fill = feat_type), color = NA, size = 11) +
  geom_gene_label(aes(label = name), size = 2.3, nudge_y = 0.05) +
  geom_vline(xintercept = c(loc(gene[1]), loc(gene[2])), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(name = NULL,
    values = c("Exon" = "#DAA520", "Insertion" = "#D94040",
               "WT_amplicon" = "lightskyblue", "junction_5" = "#40A040", "junction_3" = "#9060C0"),
    labels = c("Exon" = "Exon", "Insertion" = "Mu insertion (9-bp TSD)",
               "WT_amplicon" = "WT control (F+R)", "junction_5" = "5' junction (F+TIR6)",
               "junction_3" = "3' junction (R+TIR6)")) +
  scale_fill_manual(name = NULL,
    values = c("forward_primer" = "#40A040", "reverse_primer" = "#9060C0"),
    labels = c("forward_primer" = "Forward GSP", "reverse_primer" = "Reverse GSP")) +
  geom_vline(xintercept = loc(TSD_G), linetype = "dotted", color = "#D94040", linewidth = 0.3) +
  scale_x_continuous(labels = function(x) format(round(x + WIN_START - 1), big.mark = ",")) +
  coord_cartesian(xlim = c(loc(150183680), loc(150184230)), ylim = c(0.55, 2.45)) +
  labs(title = "Zm00004b034131 — UniformMu mu1046085: primers & amplicons",
       subtitle = "chr9 minus strand (zoom on insertion). Junction bars = genomic footprint to the TSD; product continues into Mu via TIR6.",
       x = "chr9 position (bp)") +
  theme(legend.position = "bottom", axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(ROOT, "results/annotation_map.png"), p, width = 12, height = 7.5, dpi = 150)
cat("wrote results/annotation_map.png\n")
