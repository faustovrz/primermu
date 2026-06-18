# primermu — GSP genotyping primers for UniformMu *mu1046085* in *Zm00004b034131* (W22)

Reproducible design of gene-specific primers (GSPs) to genotype the UniformMu insertion
**mu1046085** (stock UFMu-05701) in maize gene **Zm00004b034131**, in the W22 background
(Zm-W22-REFERENCE-NRGENE-2.0). GSPs are designed to be run in a multiplex with the Mu-TIR
primer **TIR6**, following the UniformMu PCR-confirmation protocol.

The whole pipeline lives in one Quarto notebook, **`mu1046085_primer_design.qmd`**, written in
R with external tools called via `system2()`.

## Quick start

```bash
quarto render mu1046085_primer_design.qmd
```

On first render it downloads the W22 reference from MaizeGDB into `data/` (≈604 MB genome +
gene GFF3), then runs the design. Subsequent renders reuse the cached reference and finish in
seconds. Renders to a self-contained `mu1046085_primer_design.html`.

## What it does

1. Verify tools (`samtools`, `primer3_core`, `e-PCR`, `curl`).
2. Fetch the W22 genome + gene GFF3 (guarded — skipped if present).
3. Parse the gene model; extract the locus (gene ± 2 kb) with `samtools faidx`.
4. Locate the 9-bp Target Site Duplication (TSD `GTTGCGCTG`) — unique in the window.
5. Design **3 forward + 3 reverse GSPs** with Primer3 (21–27 nt, GC 40–60 %, GC-clamp,
   ≥80 bp off the insertion, all WT F×R products < 500 bp).
6. Tabulate amplicons (WT control + 5′/3′ junction products with TIR6).
7. **Specificity:** forward `e-PCR` over the whole W22 genome — every pair must give a single product.
8. **Annotation map** (gggenomes, in the style of
   [sawers-rellan-labs/LH244_CRISPR](https://github.com/sawers-rellan-labs/LH244_CRISPR)).

## Repository layout

```
mu1046085_primer_design.qmd   # the pipeline (R + system calls) — source of truth
mu1046085_primer_design.html  # rendered report
results/                      # tracked deliverables
  gsp_table.tsv               #   6 GSPs (seq, Tm, GC, position, distance)
  amplicons.tsv               #   WT + junction amplicon sizes
  epcr_fwd.txt                #   e-PCR specificity hits
  primers.fa                  #   GSPs + TIR6
  annotation_map.png          #   gggenomes map
  REPORT.md                   #   standalone write-up
work/                         # only the small extracted W22 target sequences are tracked
  locus.fa                    #   gene ± 2 kb (chr9:150,176,549-150,186,207)
  gene.gff3                   #   Zm00004b034131 gene model
agent/                        # development scripts (superseded by the .qmd)
data/                         # reference downloads — NOT tracked (see .gitignore)
```

## What is NOT tracked

The reference **genome downloads** are git-ignored on purpose — they are large and fully
reproducible from the MaizeGDB URLs baked into the `.qmd`. Only the small **extracted W22 target
sequences** (`work/locus.fa`, `work/gene.gff3`) are committed.

## Requirements

- **R** (gggenomes, ggplot2, dplyr, readr, stringr, Biostrings)
- **Quarto**, **pandoc**
- CLI: **samtools**, **primer3_core**, NCBI **e-PCR**, **curl**, **bgzip**

## Key result

| GSP | Sequence (5′→3′) | Tm | 5′ @ chr9 |
|---|---|---|---|
| GSP-F1 | GACTATACCGACTACCAGCTGCAC | 62.6 | 150,184,198 |
| GSP-F2 | CGGGCGAGATAATATCGTGACATC | 61.3 | 150,184,175 |
| GSP-F3 | GGGTATCTCTTTCTCCGGGTAGTG | 61.8 | 150,184,137 |
| GSP-R1 | ACTAAGTTGTCTAGGAGCTCCACTG | 61.6 | 150,183,730 |
| GSP-R2 | GAGAGCAGTGCTAACTAAAAGGGC | 61.7 | 150,183,761 |
| GSP-R3 | GTAATGTCAGGCGTACCAGGAGG | 62.6 | 150,183,789 |

Multiplex with **TIR6** `AGAGAAGCCAACGCCAWCGCCTCYATTTCGTC` (Tm 71.7; TIR8 nested mix as fallback).
All pairs verified single-copy by e-PCR.

## References

- UniformMu PCR confirmation protocol (TIR6/TIR8, touchdown PCR): **UniformMu Resource (2011)**,
  McCarty lab — <https://download.maizegdb.org/Insertions/UniformMu/UniformMu_Resource_2011.pdf>
- Settles AM *et al.* (2004) *BMC Genomics* 5:84.
- W22 genome: Zm-W22-REFERENCE-NRGENE-2.0, MaizeGDB.
