import requests
import os
import time

input_fasta = "/mnt/proyectos/rbp_paper/results/RBPs_cdhit90.faa"
output_dir = "/mnt/proyectos/rbp_paper/results/esmfold_structures"
os.makedirs(output_dir, exist_ok=True)

MAX_LEN = 400
DELAY = 3
MAX_RETRIES = 3

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

def predict_structure(sequence, protein_id, max_retries=MAX_RETRIES):
    output_path = os.path.join(output_dir, f"{protein_id}.pdb")
    if os.path.exists(output_path):
        return "skip"

    # Truncar al C-terminal si es muy larga
    original_len = len(sequence)
    if len(sequence) > MAX_LEN:
        sequence = sequence[-MAX_LEN:]
        truncated = True
    else:
        truncated = False

    for attempt in range(max_retries):
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
                tag = f"OK [truncado C-term {original_len}→{MAX_LEN}aa]" if truncated else "OK [completa]"
                return tag
            else:
                print(f"  Intento {attempt+1} ERROR {response.status_code}")
                time.sleep(10)
        except Exception as e:
            print(f"  Intento {attempt+1} FALLO: {e}")
            time.sleep(15)
    return "failed"

sequences = parse_fasta(input_fasta)
total = len(sequences)
pendientes = [(pid, seq) for pid, seq in sequences.items()
              if not os.path.exists(os.path.join(output_dir, f"{pid}.pdb"))]

print(f"Total: {total} | Procesadas: {total - len(pendientes)} | Pendientes: {len(pendientes)}")
print(f"Proteínas >400aa que serán truncadas al C-terminal: {sum(1 for _, s in pendientes if len(s) > MAX_LEN)}")
print("-" * 60)

failed = []
for i, (pid, seq) in enumerate(pendientes):
    print(f"[{i+1}/{len(pendientes)}] {pid} ({len(seq)}aa)", end=" → ")
    result = predict_structure(seq, pid)
    print(result)
    if result == "failed":
        failed.append(pid)
    time.sleep(DELAY)

print(f"\n✓ Completado. Exitosas: {len(pendientes)-len(failed)} | Fallidas: {len(failed)}")
if failed:
    with open("failed_proteins.txt", "w") as f:
        f.write("\n".join(failed))
    print("Fallidas guardadas en: failed_proteins.txt")
