import requests
import os

# Directorios
input_fasta = "/mnt/proyectos/rbp_paper/results/RBPs_cdhit90.faa"
output_dir = "/mnt/proyectos/rbp_paper/results/esmfold_structures"
os.makedirs(output_dir, exist_ok=True)

MAX_LENGTH = 400

def parse_fasta(fasta_file):
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences

def predict_structure(sequence, protein_id):
    # Limpiar nombre para nombre de archivo
    safe_id = protein_id.replace("|", "_").replace("/", "_")
    output_path = os.path.join(output_dir, f"{safe_id}.pdb")

    if os.path.exists(output_path):
        print(f"Ya existe: {protein_id}")
        return

    try:
        response = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=sequence,
            timeout=120
        )
        if response.status_code == 200:
            with open(output_path, "w") as f:
                f.write(response.text)
            print(f"OK: {protein_id}")
        else:
            print(f"ERROR {response.status_code}: {protein_id}")
    except Exception as e:
        print(f"FALLO: {protein_id} — {e}")

# Ejecutar solo proteinas <= 400 aa
sequences = parse_fasta(input_fasta)
total = len(sequences)
procesadas = 0
saltadas = 0

for i, (pid, seq) in enumerate(sequences.items()):
    if len(seq) > MAX_LENGTH:
        saltadas += 1
        continue
    print(f"[{i+1}/{total}] Procesando: {pid} ({len(seq)} aa)")
    predict_structure(seq, pid)
    procesadas += 1

print(f"\nCompletado. Procesadas: {procesadas} | Saltadas por longitud: {saltadas}")
