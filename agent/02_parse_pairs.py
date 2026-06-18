#!/usr/bin/env python3
"""Parse primer3 output into GSP pairs with genomic coordinates and predicted
WT + TIR6 junction product sizes. Gene Zm00004b034131 is on the MINUS strand."""
import sys

OUT = "/Users/fvrodriguez/Desktop/primermu/work/p3_out.txt"
WIN_START = 150176549            # genomic coord of locus index 0 (+ strand)
TSD_G = 150183971                # genomic 1-based start of 9bp TSD
TSD_LEN = 9
TIR6 = "AGAGAAGCCAACGCCAWCGCCTCYATTTCGTC"  # 32 nt, reads OUTWARD from element
GENE_STRAND = "-"

# exon intervals (genomic, + coords) for annotation
EXONS = [(150183805,150184207),(150181691,150181980),(150179879,150180169),(150178549,150179052)]
GENE = (150178549,150184207)

def region(g):
    for i,(a,b) in enumerate(EXONS,1):
        if a<=g<=b: return f"exon{i}"
    if g<GENE[0]: return "downstream/3'-flank"   # minus strand: low coord = 3' side
    if g>GENE[1]: return "upstream/5'-flank(promoter)"
    return "intron"

d={}
for line in open(OUT):
    line=line.strip()
    if line and line!="=" and "=" in line:
        k,v=line.split("=",1); d[k]=v

n=int(d.get("PRIMER_PAIR_NUM_RETURNED","0"))
print(f"pairs returned: {n}\n")
tir_off = len(TIR6) + TSD_LEN   # bp the element side adds to a GSP+TIR6 junction product

rows=[]
for i in range(n):
    ls=d[f"PRIMER_LEFT_{i}"].split(",");  rs=d[f"PRIMER_RIGHT_{i}"].split(",")
    lstart=int(ls[0]); llen=int(ls[1]); rlast=int(rs[0]); rlen=int(rs[1])
    # primer3 indices are 0-based on + strand template
    lg = WIN_START + lstart                     # left primer 5' genomic (+ strand)
    rg = WIN_START + rlast                       # right primer 5' genomic (+ strand)
    lseq=d[f"PRIMER_LEFT_{i}_SEQUENCE"]; rseq=d[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    lTm=float(d[f"PRIMER_LEFT_{i}_TM"]); rTm=float(d[f"PRIMER_RIGHT_{i}_TM"])
    lGC=float(d[f"PRIMER_LEFT_{i}_GC_PERCENT"]); rGC=float(d[f"PRIMER_RIGHT_{i}_GC_PERCENT"])
    wt=int(d[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])
    # distance of each primer's near end to the TSD (insertion) site
    # left primer is on the HIGH-coord side? depends; compute outward distances
    dist_left  = abs((WIN_START+lstart) - (TSD_G))           # approx
    # junction product with TIR6 = (genomic dist from primer to insertion) + TSD + TIR6 len
    # left primer pairs with TIR6 reading toward it; right primer likewise
    # distance from left primer 5' to TSD start:
    dleft  = (TSD_G) - lg
    # distance from TSD end to right primer 5':
    dright = rg - (TSD_G+TSD_LEN-1)
    jl = dleft + TSD_LEN + len(TIR6)
    jr = dright + TSD_LEN + len(TIR6)
    rows.append(dict(i=i,wt=wt,
        lseq=lseq,lTm=lTm,lGC=lGC,llen=llen,lg=lg,lreg=region(lg),jl=jl,
        rseq=rseq,rTm=rTm,rGC=rGC,rlen=rlen,rg=rg,rreg=region(rg),jr=jr,
        dtm=abs(lTm-rTm)))

# On minus strand: the LEFT primer (lower template index, higher... ) — clarify F/R vs gene.
# Template is + strand. LEFT primer = + strand sense (extends toward higher coord).
# Gene is minus-strand => gene-FORWARD = decreasing coord = the RIGHT primer's direction.
# So relative to GENE orientation: RIGHT primer = gene-forward (GSP-F), LEFT primer = gene-reverse (GSP-R).
for r in rows:
    print(f"=== PAIR {r['i']}  WT(F+R) product = {r['wt']} bp   dTm={r['dtm']:.2f} ===")
    print(f"  + LEFT  primer (= gene-REVERSE / GSP-R) {r['lseq']}")
    print(f"      len {r['llen']}  Tm {r['lTm']:.1f}  GC {r['lGC']:.0f}%  5'@chr9:{r['lg']} [{r['lreg']}]  | +TIR6 junction ~{r['jl']} bp")
    print(f"  - RIGHT primer (= gene-FORWARD / GSP-F) {r['rseq']}")
    print(f"      len {r['rlen']}  Tm {r['rTm']:.1f}  GC {r['rGC']:.0f}%  5'@chr9:{r['rg']} [{r['rreg']}]  | +TIR6 junction ~{r['jr']} bp")
    print()
