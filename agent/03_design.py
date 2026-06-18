#!/usr/bin/env python3
"""Design 3 gene-FORWARD + 3 gene-REVERSE GSPs flanking the mu1046085 insertion.
Gene Zm00004b034131 is on the MINUS strand:
  - higher genomic coord = gene 5' side  -> GSP-F (primer3 RIGHT primer)
  - lower  genomic coord = gene 3' side  -> GSP-R (primer3 LEFT primer)
All forward x reverse combos kept <500 bp WT product. One pool run per side."""
import subprocess, os

ROOT="/Users/fvrodriguez/Desktop/primermu"
WIN_START=150176549
TSD_G=150183971; TSD_LEN=9
GENE=(150178549,150184207)
EXON1=(150183805,150184207)
TIR6="AGAGAAGCCAACGCCAWCGCCTCYATTTCGTC"

seq="".join(l.strip() for l in open(f"{ROOT}/work/locus.fa") if not l.startswith(">")).upper()
off=TSD_G-WIN_START
assert seq[off:off+TSD_LEN]=="GTTGCGCTG"

COMMON=[
 "PRIMER_TASK=generic","PRIMER_PICK_INTERNAL_OLIGO=0","PRIMER_FIRST_BASE_INDEX=0",
 "PRIMER_MIN_SIZE=21","PRIMER_OPT_SIZE=24","PRIMER_MAX_SIZE=27",
 "PRIMER_MIN_TM=60.0","PRIMER_OPT_TM=62.0","PRIMER_MAX_TM=65.0",
 "PRIMER_MIN_GC=40.0","PRIMER_OPT_GC_PERCENT=50.0","PRIMER_MAX_GC=60.0",
 "PRIMER_GC_CLAMP=1","PRIMER_MAX_POLY_X=4","PRIMER_MAX_NS_ACCEPTED=0",
 "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0","PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0",
 "PRIMER_EXPLAIN_FLAG=1","PRIMER_NUM_RETURN=40",
]

def run(side):
    MIN_OFF=80   # keep every GSP >=~80 bp from the TSD: clean, gel-resolvable junction bands
    if side=="F":   # gene-forward = RIGHT primer, ABOVE insertion, keep in exon1
        incl_start=(off+TSD_LEN-1)+MIN_OFF
        incl_len=(EXON1[1]-WIN_START)-incl_start+1     # up to exon1 5' boundary
        pick=["PRIMER_PICK_LEFT_PRIMER=0","PRIMER_PICK_RIGHT_PRIMER=1"]
    else:           # gene-reverse = LEFT primer, BELOW insertion (exon1 + intron1)
        incl_start=off-245
        incl_len=(off-MIN_OFF)-incl_start
        pick=["PRIMER_PICK_LEFT_PRIMER=1","PRIMER_PICK_RIGHT_PRIMER=0"]
    lines=[f"SEQUENCE_ID=mu1046085_{side}",f"SEQUENCE_TEMPLATE={seq}",
           f"SEQUENCE_INCLUDED_REGION={incl_start},{incl_len}"]+pick+COMMON+["="]
    inp=f"{ROOT}/work/p3_{side}.txt"; outp=f"{ROOT}/work/p3_{side}_out.txt"
    open(inp,"w").write("\n".join(lines)+"\n")
    subprocess.run(f"primer3_core {inp} > {outp}",shell=True,check=True)
    d={}
    for ln in open(outp):
        ln=ln.strip()
        if "=" in ln and ln!="=": k,v=ln.split("=",1); d[k]=v
    tag="RIGHT" if side=="F" else "LEFT"
    n=int(d.get(f"PRIMER_{tag}_NUM_RETURNED","0"))
    cands=[]
    for i in range(n):
        pos=d[f"PRIMER_{tag}_{i}"].split(","); p5=int(pos[0]); plen=int(pos[1])
        g5=WIN_START+p5                      # 5' genomic coord
        cands.append(dict(seq=d[f"PRIMER_{tag}_{i}_SEQUENCE"],tm=float(d[f"PRIMER_{tag}_{i}_TM"]),
            gc=float(d[f"PRIMER_{tag}_{i}_GC_PERCENT"]),length=plen,g5=g5,rank=i))
    return cands

def select(cands,n=3,min_gap=22,tm_band=5.0):
    """greedy: take by primer3 rank (best first), require 5' spacing so the 3 are
    genuinely different priming sites. Tm kept within the 60-65 window already;
    each recommended pair is matched to <1C below."""
    chosen=[]
    for c in sorted(cands,key=lambda x:x["rank"]):
        if any(abs(c["g5"]-o["g5"])<min_gap for o in chosen): continue
        trial=chosen+[c]
        if max(x["tm"] for x in trial)-min(x["tm"] for x in trial)>tm_band: continue
        chosen.append(c)
        if len(chosen)==n: break
    return chosen

def region(g):
    if EXON1[0]<=g<=EXON1[1]: return "exon1"
    if g>GENE[1]: return "promoter/5'"
    if g<GENE[0]: return "3'flank"
    return "intron1"

F=select(run("F")); R=select(run("R"))
print(f"selected {len(F)} forwards, {len(R)} reverses\n")

def jct(side,g5):
    if side=="F": d=g5-(TSD_G+TSD_LEN-1)
    else:         d=(TSD_G)-g5
    return d+TSD_LEN+len(TIR6)

names={}
print("GENE-FORWARD GSPs (5' side, pair with TIR6 -> 5'-junction amplicon):")
for k,c in enumerate(sorted(F,key=lambda x:-x["g5"]),1):
    nm=f"GSP-F{k}"; names[id(c)]=nm
    print(f"  {nm}  {c['seq']}  len{c['length']} Tm{c['tm']:.1f} GC{c['gc']:.0f}%  5'@chr9:{c['g5']} [{region(c['g5'])}]  +TIR6~{jct('F',c['g5'])}bp")
print("\nGENE-REVERSE GSPs (3' side, pair with TIR6 -> 3'-junction amplicon):")
for k,c in enumerate(sorted(R,key=lambda x:x["g5"]),1):
    nm=f"GSP-R{k}"; names[id(c)]=nm
    print(f"  {nm}  {c['seq']}  len{c['length']} Tm{c['tm']:.1f} GC{c['gc']:.0f}%  5'@chr9:{c['g5']} [{region(c['g5'])}]  +TIR6~{jct('R',c['g5'])}bp")

allp=[(names[id(c)],c['seq'],c['tm']) for c in F]+[(names[id(c)],c['seq'],c['tm']) for c in R]
print(f"\nTm spread across all selected: {min(p[2] for p in allp):.1f}-{max(p[2] for p in allp):.1f} (range {max(p[2] for p in allp)-min(p[2] for p in allp):.2f})")

# all 9 WT combos
print("\nWT (F x R) product sizes, bp:")
print("        "+"  ".join(f"GSP-R{j}" for j in range(1,len(R)+1)))
Fs=sorted(F,key=lambda x:-x["g5"]); Rs=sorted(R,key=lambda x:x["g5"])
worst=0
for i,f in enumerate(Fs,1):
    row=[]
    for r in Rs:
        sz=f["g5"]-r["g5"]+1; worst=max(worst,sz); row.append(f"{sz:5d}")
    print(f"  GSP-F{i} "+"  ".join(row))
print(f"largest WT product across combos: {worst} bp (target <500)")

# write primers fasta + re-PCR STS (all 9 combos)
os.makedirs(f"{ROOT}/results",exist_ok=True)
with open(f"{ROOT}/results/primers.fa","w") as fa:
    for c in Fs: fa.write(f">{names[id(c)]}\n{c['seq']}\n")
    for c in Rs: fa.write(f">{names[id(c)]}\n{c['seq']}\n")
    fa.write(f">TIR6\n{TIR6}\n")
with open(f"{ROOT}/work/epcr_sts.txt","w") as sts:
    for i,f in enumerate(Fs,1):
        for j,r in enumerate(Rs,1):
            # e-PCR STS: id  left  right  size-range   (tab separated)
            sts.write(f"F{i}xR{j}\t{f['seq']}\t{r['seq']}\t100-600\n")
with open(f"{ROOT}/results/gsp_table.tsv","w") as t:
    t.write("name\tside\tsequence\tlen\tTm\tGC%\t5p_chr9\tregion\tdist_to_insertion\tTIR6_junction_approx\n")
    for c in Fs:
        d=c['g5']-(TSD_G+TSD_LEN-1)
        t.write(f"{names[id(c)]}\tforward(gene-fwd)\t{c['seq']}\t{c['length']}\t{c['tm']:.1f}\t{c['gc']:.0f}\t{c['g5']}\t{region(c['g5'])}\t{d}\t{jct('F',c['g5'])}\n")
    for c in Rs:
        d=(TSD_G)-c['g5']
        t.write(f"{names[id(c)]}\treverse(gene-rev)\t{c['seq']}\t{c['length']}\t{c['tm']:.1f}\t{c['gc']:.0f}\t{c['g5']}\t{region(c['g5'])}\t{d}\t{jct('R',c['g5'])}\n")
print("\nwrote results/primers.fa, results/gsp_table.tsv and work/epcr_sts.txt (9 F x R combos)")
