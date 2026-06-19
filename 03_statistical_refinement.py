import pandas as pd
import os

CSV = "/mnt/proyectos/rbp_paper/data/filtered/filtered_metadata.csv"
if not os.path.exists(CSV):
    print("Error: No existe el archivo de metadata.")
    exit()

df = pd.read_csv(CSV)

# Calculo de Densidad de CDS (Proteinas por cada 1000 bp)
df["cds_density"] = (df["cds_count"] / df["length_bp"]) * 1000

print("\n--- REPORTE ESTRATÉGICO DE DENSIDAD GENÓMICA ---")
print(f"Total de fagos analizados: {len(df)}")
print("\nDensidad promedio (CDS/kbp) por patógeno:")
summary = df.groupby("pathogen")["cds_density"].mean().round(3)
print(summary)

# Guardar datos finales para el paper
df.to_csv("/mnt/proyectos/rbp_paper/data/filtered/refined_phage_data.csv", index=False)
print("\n[OK] Datos refinados guardados. Listos para tablas y figuras.")
