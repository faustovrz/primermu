# Designing Gene-Specific Primers Around a Mu Insertion (W22 maize)

A practical guide to interpreting a *Mutator* (Mu) insertion annotation and
designing diagnostic genotyping primers around it.

---

## 1. The "9 bp" in the annotation is **not** the transposon length

Mu elements are large — **1.4–4.9 kb** (e.g. Mu1 ≈ 1.4 kb; the autonomous
MuDR ≈ 4.9 kb). What is invariant across the entire Mu family is a pair of
**~215 bp Terminal Inverted Repeats (TIRs)** at each end.

When Mu inserts, it duplicates **9 bp of host genomic sequence** at the target
site — the **Target Site Duplication (TSD)**. This 9-bp duplication is the
hallmark signature of a Mu insertion and is what the genomic track is
annotating.

| Annotation says | Actually means |
|---|---|
| "insertion is 9 bp long" | The **9-bp host Target Site Duplication (TSD)**, the footprint of insertion — *not* the element size |
| "primer inside the transposon" | The **universal Mu-TIR primer**, which anneals within the ~215 bp TIR |

---

## 2. Can you deduce the insertion orientation from the annotation?

**No — and you don't need to.**

Because the two ends of Mu are **inverted** repeats, a primer sitting in the
TIR reads *outward into flanking genomic DNA at **both** ends of the element*,
no matter which way the element landed. Therefore:

1. A TIR primer cannot, by itself, report the element's orientation.
2. The 9-bp TSD pinpoints the **insertion site** but carries **no orientation
   information** — it is just duplicated host sequence. Any "strand" shown in
   the track usually reflects only which flank was sequenced when the insertion
   was recovered, not a constraint on your assay.
3. The assay is **diagnostic by primer pairing**: flank the site with two
   gene-specific primers and pair each with the TIR primer.

---

## 3. Structure of the insertion (text diagram)

Gene orientation = transcription left → right.

```
   GENE TRANSCRIPTION  ───────────────────────────────────────►

   ── UPSTREAM (5') flank ──┐                              ┌── DOWNSTREAM (3') flank ──
   genomic DNA              │                              │              genomic DNA
   ─────────────────[ TSD ]┌┴──────── Mu ELEMENT ─────────┴┐[ TSD ]────────────────
                     9 bp  │ ▐TIR_L▌══ internal seq ══▐TIR_R▌ │ 9 bp
                           └──────────────────────────────────┘
                                ▲                          ▲
                  ◄──Mu-TIR primer            Mu-TIR primer──►
                  (extends OUT, leftward)     (extends OUT, rightward)

   GSP-F ──►                                                        ◄── GSP-R
 (forward, in upstream flank,                        (reverse, in downstream flank,
  points toward insertion)                            points toward insertion)
```

**Reading the arrows**

- **TIR_L** and **TIR_R** are reverse-complements of each other (inverted), so
  the *same* TIR primer sequence anneals at both ends and **always extends
  outward** into the flank.
- **GSP-F** — gene-specific **forward** primer (sense of the gene), in the
  upstream flank, points rightward toward the insertion.
- **GSP-R** — gene-specific **reverse** primer, in the downstream flank, points
  leftward toward the insertion.

---

## 4. How the pairings work (the diagnostic logic)

| Primer pair | Detects | Product |
|---|---|---|
| **GSP-F + GSP-R** | Wild-type / no-insertion allele | Normal-size amplicon spanning the empty site. With an insertion present it fails or gives a much larger band. |
| **GSP-F + Mu-TIR** | Left junction of the insertion | Band only if insertion present |
| **GSP-R + Mu-TIR** | Right junction of the insertion | Band only if insertion present |

Because the TIR primer faces outward at **both** ends, **either** TIR pairing
can amplify across its junction — which is exactly why orientation is
unnecessary. Standard setup:

- **GSP-F + GSP-R** → wild-type allele
- **TIR + GSP-F** (and/or **TIR + GSP-R**) → insertion allele

**Genotype calls**

| Bands observed | Genotype |
|---|---|
| WT band only | Homozygous wild-type |
| TIR-junction band only | Homozygous insertion |
| Both bands | Heterozygous |

---

## 5. Practical primer-design tips

- Place each GSP **~100–300 bp** from the 9-bp TSD so the TIR-junction products
  are an easy size to score and resolve from the WT band.
- Match GSP **Tm to the Mu-TIR primer** (typically **~58–62 °C**).
- Design GSPs at **18–25 nt**.
- Center the WT amplicon on the TSD coordinates from the track.
- Run **both** TIR+GSP-F and TIR+GSP-R the first time — whichever gives a clean
  junction band becomes your routine genotyping pair.
- **Sequence the junction product** to confirm the exact insertion site: you
  will see the 9-bp TSD abutting the TIR sequence.

---

## 6. Glossary

| Term | Definition |
|---|---|
| **Mu / Mutator** | A family of high-copy DNA transposons in maize used for insertional mutagenesis |
| **TIR** | Terminal Inverted Repeat (~215 bp); conserved at both ends of all Mu elements; basis of the universal TIR primer |
| **TSD** | Target Site Duplication; 9 bp of host sequence duplicated upon Mu insertion |
| **Mu-TIR primer** | Universal primer annealing within the TIR; extends outward into flanking DNA at both ends |
| **GSP** | Gene-Specific Primer; anchored in unique flanking genomic sequence |
