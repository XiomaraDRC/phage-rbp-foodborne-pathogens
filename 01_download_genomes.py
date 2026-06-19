import os, time
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "xiomara.romero.calle.d@gmail.com"
OUTPUT_DIR = "/mnt/proyectos/rbp_paper/data/raw"
os.makedirs(OUTPUT_DIR, exist_ok=True)

QUERIES = {
    "Salmonella": "Salmonella phage complete genome",
    "Ecoli": "Escherichia phage complete genome",
    "Listeria": "Listeria phage complete genome",
}

MAX_GENOMES = {"Salmonella": 100, "Ecoli": 100, "Listeria": 60}
MIN_LENGTH = 10000
MAX_LENGTH = 500000

all_results = []
for pathogen, query in QUERIES.items():
    print(f"PATOGENO: {pathogen}")
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=MAX_GENOMES[pathogen]*3)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"  Encontrados: {len(ids)} IDs")
    out_dir = os.path.join(OUTPUT_DIR, pathogen)
    os.makedirs(out_dir, exist_ok=True)
    downloaded = 0
    for uid in ids:
        if downloaded >= MAX_GENOMES[pathogen]:
            break
        try:
            handle = Entrez.efetch(db="nucleotide", id=uid, rettype="gb", retmode="text")
            rec = SeqIO.read(handle, "genbank")
            handle.close()
            seq_len = len(rec.seq)
            if not (MIN_LENGTH <= seq_len <= MAX_LENGTH):
                time.sleep(0.4)
                continue
            out_path = os.path.join(out_dir, f"{rec.id}.gb")
            SeqIO.write(rec, out_path, "genbank")
            all_results.append({"accession": rec.id, "description": rec.description, "length_bp": seq_len, "pathogen": pathogen})
            downloaded += 1
            if downloaded % 10 == 0:
                print(f"    {downloaded} descargados...")
            time.sleep(0.4)
        except Exception as e:
            time.sleep(1)
    print(f"  Total: {downloaded}")

if all_results:
    df = pd.DataFrame(all_results)
    df.to_csv(os.path.join(OUTPUT_DIR, "genome_metadata.csv"), index=False)
    print(f"TOTAL FINAL: {len(df)} genomas")
    print(df.groupby("pathogen")["accession"].count())
