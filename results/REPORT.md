# GSP genotyping primers for UniformMu mu1046085 in *Zm00004b034131* (W22)

**Date:** 2026-06-18 · **Reference:** Zm-W22-REFERENCE-NRGENE-2.0 · **Stock:** UFMu-05701

## 1. Insertion / target summary

| Item | Value |
|---|---|
| Gene | **Zm00004b034131**, chr9:150,178,549–150,184,207, **minus strand** |
| Insertion (catalog) | B73 v5 NAM chr9:151,737,374–151,737,382; `W22_id chr9_150183971` (High, 100%) |
| **9-bp Target Site Duplication (TSD)** | `GTTGCGCTG` |
| **Insertion site in W22** | chr9:**150,183,971–150,183,979** — **unique** exact match in gene±2 kb |
| Location in gene model | **Exon 1 / CDS** (exon1 = 150,183,805–150,184,207) — classic Mu 5′ bias, exonic anchor |

Reminder: "length 9 bp" in the track = the TSD, **not** the element size. Mu elements are 1.4–4.9 kb with ~215 bp TIRs; the TIR6 primer reads **outward** from both element ends → two distinct junction amplicons (5′ and 3′).

## 2. Designed gene-specific primers (3 forward + 3 reverse)

Forward = gene-forward (5′ side, higher coord); Reverse = gene-reverse (3′ side, lower coord). Primer3 length 21–27 nt, GC 40–60%, GC-clamp, no thermodynamic config (classic Tm model).

| Name | Sequence (5′→3′) | Len | Tm | GC% | 5′ @ chr9 | Region | dist. to insertion |
|---|---|---|---|---|---|---|---|
| **GSP-F1** | `GACTATACCGACTACCAGCTGCAC` | 24 | 62.6 | 54 | 150,184,198 | exon1 | 219 bp |
| **GSP-F2** | `CGGGCGAGATAATATCGTGACATC` | 24 | 61.3 | 50 | 150,184,175 | exon1 | 196 bp |
| **GSP-F3** | `GGGTATCTCTTTCTCCGGGTAGTG` | 24 | 61.8 | 54 | 150,184,137 | exon1 | 158 bp |
| **GSP-R1** | `ACTAAGTTGTCTAGGAGCTCCACTG` | 25 | 61.6 | 48 | 150,183,730 | intron1 | 241 bp |
| **GSP-R2** | `GAGAGCAGTGCTAACTAAAAGGGC` | 24 | 61.7 | 50 | 150,183,761 | intron1 | 210 bp |
| **GSP-R3** | `GTAATGTCAGGCGTACCAGGAGG` | 23 | 62.6 | 57 | 150,183,789 | intron1 | 182 bp |

All 6 in the recommended 60–65 °C window; pool Tm spread = 1.35 °C. Forwards sit in exon 1; reverses fall in intron 1 (the exonic stretch below the insertion is short/GC-rich). Designed on W22 directly, so the protocol's B73→W22 polymorphism caveat does not apply.

## 3. Fixed partner: Mu-TIR primer

**TIR6** `AGAGAAGCCAACGCCAWCGCCTCYATTTCGTC` (Tm 71.7, degenerate W=A/T, Y=C/T). The high Tm vs the ~62 °C GSPs is handled by the touchdown protocol (it pairs the same way in the published assay). **Fallback:** TIR8 nested mix (TIR8.1–8.4) for two-stage specificity if TIR6 gives no product.

## 4. Amplicon sizes

**WT control (GSP-F + GSP-R), exact, all <500 bp** — any forward × any reverse works:

|  | GSP-R1 | GSP-R2 | GSP-R3 |
|---|---|---|---|
| **GSP-F1** | 469 | 438 | 410 |
| **GSP-F2** | 446 | 415 | 387 |
| **GSP-F3** | 408 | 377 | 349 |

**Insertion junction bands (GSP + TIR6), approximate** — each = (dist. to insertion) + 9 (TSD) + 32 (TIR6) + a constant intra-TIR offset (same for all; resolvable on 1% agarose):
- 5′-junction: GSP-F1 ~260 · GSP-F2 ~237 · GSP-F3 ~199 bp
- 3′-junction: GSP-R1 ~282 · GSP-R2 ~251 · GSP-R3 ~223 bp

**Three Tm-matched pairs to order (ΔTm < 1 °C):**
| Pair | Forward | Reverse | WT product | ΔTm |
|---|---|---|---|---|
| A | GSP-F3 | GSP-R2 | 377 bp | 0.1 |
| B | GSP-F2 | GSP-R1 | 446 bp | 0.3 |
| C | GSP-F1 | GSP-R3 | 410 bp | 0.0 |

## 5. Specificity — forward e-PCR (genome-wide, W22)

`e-PCR -n 1 -g 1 -m 200 -w 7` over the whole W22 genome. **All 9 F×R combos → exactly one amplicon**, all at the chr9 target, 0 mismatches / 0 gaps. No paralogous or duplicate amplicons anywhere → no false-heterozygous risk from gene duplication.

## 6. Assay & interpretation

Multiplex per plant: **GSP-F + GSP-R + TIR6**.
- **WT band only** (377/410/446 bp) → homozygous wild-type
- **Junction band(s) only** (GSP+TIR6) → homozygous insertion
- **Both** → heterozygous

Workflow (McCarty/UniformMu): (1) confirm each GSP-F+GSP-R pair on **W22 control DNA**; (2) run **both** GSP-F+TIR6 and GSP-R+TIR6 on UFMu plant DNA (occasionally only one side amplifies); (3) **sequence** junction products to confirm the site (TSD `GTTGCGCTG` abutting the TIR).

**Touchdown PCR (Table 2):** 94 °C 1′; [94 °C 25″ / 62 °C 30″ / 72 °C 1′] ×8–10; [94 °C 25″ / 56 °C 30″ / 72 °C 1′] ×~27; 72 °C 5′; 4 °C. Add **5% DMSO** (fresh, −20 °C aliquots) or PCRx Enhancer for GC-rich maize template.

## 7. Files
- `results/primers.fa` — 6 GSPs + TIR6
- `results/gsp_table.tsv` — primer table
- `results/epcr_fwd.txt` — e-PCR specificity hits
- `work/locus.fa` — gene ±2 kb template · `work/p3_*` — Primer3 I/O · `agent/*` — scripts

## 8. References
- UniformMu PCR confirmation protocol (TIR6/TIR8, touchdown PCR, genotyping): **UniformMu Resource (2011)**, McCarty lab — https://download.maizegdb.org/Insertions/UniformMu/UniformMu_Resource_2011.pdf
- Settles AM et al. (2004) *BMC Genomics* 5:84.
- W22 genome: Zm-W22-REFERENCE-NRGENE-2.0, MaizeGDB.
