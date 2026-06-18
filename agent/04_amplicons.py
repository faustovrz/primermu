#!/usr/bin/env python3
"""Report the diagnostic amplicon sizes for the current Primer3 design (p3_out.txt).
Three amplicon classes in the multiplex (GSP-F + GSP-R + TIR6):
  1. WT control      = GSP-F + GSP-R         (exact, on no-insertion allele)
  2. 5'-junction     = GSP-F + TIR6          (mutant allele, 5' side of gene)
  3. 3'-junction     = GSP-R + TIR6          (mutant allele, 3' side of gene)
Gene Zm00004b034131 = MINUS strand; insertion TSD GTTGCGCTG at chr9:150183971-150183979."""
OUT="/Users/fvrodriguez/Desktop/primermu/work/p3_out.txt"
WIN=150176549; TSD_G=150183971; TSD_LEN=9
TIR6_LEN=32   # primer length

d={}
for ln in open(OUT):
    ln=ln.strip()
    if "=" in ln and ln!="=": k,v=ln.split("=",1); d[k]=v
n=int(d.get("PRIMER_PAIR_NUM_RETURNED","0"))

print("Gene is MINUS strand: RIGHT primer (higher coord) = GSP-F (gene-forward),")
print("LEFT primer (lower coord) = GSP-R (gene-reverse).\n")
for i in range(n):
    ls=d[f"PRIMER_LEFT_{i}"].split(","); rs=d[f"PRIMER_RIGHT_{i}"].split(",")
    lg=WIN+int(ls[0]); rg=WIN+int(rs[0])          # 5' genomic coords
    wt=int(d[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])
    # genomic distance from each primer's 5' end to the near edge of the TSD/insertion
    dF=rg-(TSD_G+TSD_LEN-1)        # GSP-F (rg high) down to TSD 3' edge
    dR=(TSD_G)-lg                  # GSP-R (lg low) up to TSD 5' edge
    # junction product = genomic distance + 9bp TSD + element-side contribution
    #   element side = TIR6 primer (32) + unknown intra-TIR span between TIR6 3' end and junction (>=0)
    jF_min=dF+TSD_LEN+TIR6_LEN
    jR_min=dR+TSD_LEN+TIR6_LEN
    print(f"PAIR {i}:")
    print(f"  GSP-F 5'@chr9:{rg}  (={dF} bp from insertion)   GSP-R 5'@chr9:{lg}  (={dR} bp from insertion)")
    print(f"  1) WT control (F+R)        = {wt} bp   [exact]")
    print(f"  2) 5'-junction (F+TIR6)    = {jF_min} bp + intra-TIR offset   [>= {jF_min}]")
    print(f"  3) 3'-junction (R+TIR6)    = {jR_min} bp + intra-TIR offset   [>= {jR_min}]")
    print()
print("NOTE: WT sizes are exact. Junction sizes = (GSP->insertion distance) + 9 (TSD) + 32 (TIR6)")
print("plus the TIR span between the TIR6 3' end and the element/genomic junction (protocol notes")
print("TIR sequence exists downstream of the primer). That intra-TIR offset is constant for every")
print("GSP and needs the Mu reference TIR sequence to pin exactly.")
