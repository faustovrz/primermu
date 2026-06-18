#!/usr/bin/env python3
"""Build a consolidated amplicon table for the mu1046085 / Zm00004b034131 assay.
Three amplicon classes:
  WT      = GSP-F + GSP-R         on the no-insertion allele   (exact, from primer3/e-PCR)
  5'-jct  = GSP-F + TIR6          on the Mu allele             (approx, +const intra-TIR offset)
  3'-jct  = GSP-R + TIR6          on the Mu allele             (approx)
"""
ROOT="/Users/fvrodriguez/Desktop/primermu"

# load GSP table
g={}
for ln in open(f"{ROOT}/results/gsp_table.tsv"):
    f=ln.rstrip("\n").split("\t")
    if f[0]=="name": continue
    g[f[0]]=dict(side=f[1],seq=f[2],tm=float(f[4]),g5=int(f[6]),jct=int(f[9]))

# exact WT sizes from e-PCR output (chr9  F#xR#  strand from to size ...)
wt={}
for ln in open(f"{ROOT}/results/epcr_fwd.txt"):
    if ln.startswith("#") or not ln.strip(): continue
    c=ln.split("\t"); pair=c[1]; size=int(c[5].split("/")[0]); wt[pair]=size

Fs=["GSP-F1","GSP-F2","GSP-F3"]; Rs=["GSP-R1","GSP-R2","GSP-R3"]
rows=[]
# WT amplicons (9)
for i,F in enumerate(Fs,1):
    for j,R in enumerate(Rs,1):
        rows.append((f"{F} + {R}","WT control","no-insertion allele",str(wt[f"F{i}xR{j}"]),"exact"))
# junction amplicons
for F in Fs:
    rows.append((f"{F} + TIR6","5' junction","Mu insertion allele",f"~{g[F]['jct']}","approx"))
for R in Rs:
    rows.append((f"{R} + TIR6","3' junction","Mu insertion allele",f"~{g[R]['jct']}","approx"))

# write tsv
with open(f"{ROOT}/results/amplicons.tsv","w") as t:
    t.write("primer_pair\tamplicon_type\ttemplate\tsize_bp\tprecision\n")
    for r in rows: t.write("\t".join(r)+"\n")

# pretty print
w=[max(len(r[k]) for r in rows+[("primer_pair","amplicon_type","template","size_bp","precision")]) for k in range(5)]
hdr=("primer_pair","amplicon_type","template","size_bp","precision")
def line(r): print("  ".join(r[k].ljust(w[k]) for k in range(5)))
line(hdr); line(tuple("-"*w[k] for k in range(5)))
for r in rows: line(r)
print(f"\nwrote results/amplicons.tsv ({len(rows)} amplicons)")
