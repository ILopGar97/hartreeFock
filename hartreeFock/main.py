# -*- coding: utf-8 -*-
from file_manager import load_system, csv_export_atomic_system  # Asegurate que el archivo se llame hf_loader.py
BASE_PATH = "COEFS"  # Carpeta con los archivos .ani, .neu, .cat, etc.

def main():
    
    neutral_atoms = load_system(Q=0, base_path=BASE_PATH)
    cation_atoms = load_system(Q=1, base_path=BASE_PATH)
    anion_atoms = load_system(Q=-1, base_path=BASE_PATH)
    N2_isoelectronic = load_system(Q=2, base_path=BASE_PATH)

    csv_export_atomic_system(
        atomic_data=neutral_atoms,
        output_file="neutral_atoms.csv",
        delimiter=";"
    )
    csv_export_atomic_system(
        atomic_data=cation_atoms,
        output_file="cation_atoms.csv",
        delimiter=";"
    )
    csv_export_atomic_system(
        atomic_data=anion_atoms,
        output_file="anion_atoms.csv",
        delimiter=";"
    )
    csv_export_atomic_system(
        atomic_data=N2_isoelectronic,
        output_file="N2_isoelectronic.csv",
        delimiter=";"
    )
    


if __name__ == "__main__":
    main()
