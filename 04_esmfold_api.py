import requests
import os

# Directorios
input_fasta = "/mnt/proyectos/rbp_paper/results/RBPs_cdhit90.faa"
output_dir = "/mnt/proyectos/rbp_paper/results/esmfold_structures"
os.makedirs(output_dir, exist_ok=True)

# Leer secuencias
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

# Predecir estructura con ESMFold API
def predict_structure(sequence, protein_id):
    output_path = os.path.join(output_dir, f"{protein_id}.pdb")
    
    # Si ya existe, saltamos
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

# Ejecutar
sequences = parse_fasta(input_fasta)
total = len(sequences)

for i, (pid, seq) in enumerate(sequences.items()):
    print(f"[{i+1}/{total}] Procesando: {pid}")
    predict_structure(seq, pid)

print("Completado.")
