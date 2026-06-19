import os, glob
from Bio import SeqIO
FILTERED_DIR = "/mnt/proyectos/rbp_paper/data/filtered"
OUT_DIR = "/mnt/proyectos/rbp_paper/data/rbp_candidates"
os.makedirs(OUT_DIR, exist_ok=True)
PATHOGENS = ["Salmonella", "Ecoli", "Listeria"]
for pathogen in PATHOGENS:
    in_dir = os.path.join(FILTERED_DIR, pathogen)
    out_fasta = os.path.join(OUT_DIR, pathogen + "_proteins.faa")
    files = glob.glob(os.path.join(in_dir, "*.gb"))
    total = 0
    fh = open(out_fasta, "w")
    for f in files:
        try:
            rec = SeqIO.read(f, "genbank")
            gid = str(rec.id)
            for feat in rec.features:
                if feat.type == "CDS" and "translation" in feat.qualifiers:
                    seq = feat.qualifiers["translation"][0]
                    pid = str(feat.qualifiers.get("protein_id", [gid])[0])
                    prod = str(feat.qualifiers.get("product", ["unknown"])[0]).replace(" ","_")
                    fh.write(">" + pid + "|" + gid + "|" + prod + chr(10))
                    fh.write(seq + chr(10))
                    total += 1
        except:
            continue
    fh.close()
    print(pathogen, "proteinas:", total)
