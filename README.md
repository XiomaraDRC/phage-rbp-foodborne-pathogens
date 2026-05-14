# Phage RBP Identification Pipeline

Computational pipeline for receptor-binding protein (RBP) identification 
from phage genomes targeting Salmonella, E. coli, and Listeria monocytogenes.

## Pipeline steps
1. Genome quality filtering (Biopython)
2. Protein extraction
3. RBP prediction (PhageRBPdetect/ESM-2)
4. Structural prediction (ESMFold)
5. Structural validation (Foldseek vs PDB100)
6. Molecular docking (PatchDock)

## Requirements
Python 3.10+, Biopython 1.87, CD-HIT

## Citation
Romero-Calle et al. (2026). In preparation.
