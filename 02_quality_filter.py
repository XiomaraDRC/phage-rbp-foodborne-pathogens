import os, glob
import pandas as pd
from Bio import SeqIO

RAW_DIR = "/mnt/proyectos/rbp_paper/data/raw"
OUT_DIR = "/mnt/proyectos/rbp_paper/data/filtered"
os.makedirs(OUT_DIR, exist_ok=True)
PATHOGENS = ["Salmonella", "Ecoli", "Listeria"]
MIN_CDS = 20
MIN_LENGTH = 15000
results = []
for pathogen in PATHOGENS:
    in_dir = os.path.join(RAW_DIR, pathogen)
    out_dir = os.path.join(OUT_DIR, pathogen)
    os.makedirs(out_dir, exist_ok=True)
    files = glob.glob(os.path.join(in_dir, "*.gb"))
    passed = 0
    failed = 0
    for f in files:
        try:
            rec = SeqIO.read(f, "genbank")
            cds_count = sum(1 for ft in rec.features if ft.type == "CDS")
            seq_len = len(rec.seq)
            if seq_len < MIN_LENGTH or cds_count < MIN_CDS:
                failed += 1
                continue
            out_path = os.path.join(out_dir, os.path.basename(f))
            SeqIO.write(rec, out_path, "genbank")
            results.append({"accession": rec.id, "pathogen": pathogen, "length_bp": seq_len, "cds_count": cds_count})
            passed += 1
        except:
            failed += 1
    print(pathogen, "pasaron:", passed, "descartados:", failed)
df = pd.DataFrame(results)
df.to_csv("/mnt/proyectos/rbp_paper/data/filtered/filtered_metadata.csv", index=False)
print("TOTAL:", len(df))
print(df.groupby("pathogen")["accession"].count())
