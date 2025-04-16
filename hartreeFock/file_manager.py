import os
import math
from typing import List, Dict, Tuple
from data_structures import AtomicData, Orbital
import csv
# ------------------------

L_MAP = {'S': 0, 'P': 1, 'D': 2, 'F': 3}


# ------------------------
# LECTURA DE ARCHIVOS
# ------------------------

def read_atoms(path: str) -> List[Tuple[int, str, int, int, int]]:
    atoms = []
    with open(path) as f:
        for line in f:
            t = line.strip().split()
            if len(t) >= 5:
                Z = int(t[0])
                config = t[1]
                PZ = int(t[2]) - 1  # Fortran → Python index
                QZ = int(t[3]) - 1
                SZ = int(t[4]) - 1
                atoms.append((Z, config, PZ, QZ, SZ))
    return atoms


def read_control(path: str, recl: int = 12) -> List[Tuple[int, int, str, int]]:
    control = []
    with open(path, "r") as f:
        while True:
            registro = f.read(recl)
            if not registro:
                break
            tokens = registro.strip().split()
            if len(tokens) >= 4:
                kz = int(tokens[0])
                indicador = int(tokens[1])
                nl = tokens[2].strip()
                mz = int(tokens[3])
                control.append((kz, indicador, nl, mz))
    return control


def read_nalfa(path: str) -> List[Tuple[int, float]]:
    with open(path) as f:
        tokens = f.read().split()
        return [(int(float(tokens[i])), float(tokens[i+1])) for i in range(0, len(tokens), 2)]


def read_coefficients(path: str) -> List[float]:
    with open(path) as f:
        return [float(t) for t in f.read().split()]


def read_energy(path: str) -> Dict[int, float]:
    with open(path) as f:
        return {int(t[0]): float(t[1]) for line in f if (t := line.strip().split())}


def read_valence(path: str) -> Dict[int, Tuple[int, int]]:
    with open(path) as f:
        return {int(t[0]): (int(t[1]), int(t[2])) for line in f if (t := line.strip().split())}


def read_potioniz(path: str) -> Dict[int, float]:
    with open(path) as f:
        return {int(t[0]): float(t[1]) for line in f if (t := line.strip().split())}

def extension_from_Q(NQ: int) -> str:
    """
    Retorna la extensión del archivo (.ani, .neu, .cat, .2 ... .9) dada Z y Q.
    Usa la cantidad de electrones N = Z + Q como selector principal.
    """
   
    if NQ == -1:
        return "ani"
    elif NQ == 0:
        return "neu"
    elif NQ == 1:
        return "cat"
    elif 2 <= NQ <= 9:
        return str(NQ)  # series isoelectrónicas
    else:
        raise ValueError(f"[ERROR] N={NQ} fuera del rango válido para familias (2–9)")

# ------------------------
# PROCESAMIENTO POR SISTEMA
# ------------------------

def assemble_atomic_data(
    Q: int,
    atoms: List[Tuple[int, str, int, int, int]],  # ← ahora incluye PZ, QZ, SZ
    control: List[Tuple[int, int, str, int]],
    nalfa: List[Tuple[int, float]],
    coef: List[float],
    energy: Dict[int, float],
    potioniz: Dict[int, float],
    valence: Dict[int, Tuple[int, int]],
) -> Dict[Tuple[int, int], AtomicData]:
    data = {}

    for Z, config, PZ, QZ, SZ in atoms:
        atom = AtomicData(Z=Z, Q=Q, N= Z-Q, config=config)
        atom.energy = energy.get(Z, 0.0)
        if potioniz and valence is not None:
            atom.ionization = potioniz.get(Z, 0.0)
            atom.valence_nl = valence.get(Z, (0, 0))
        else:
            atom.ionization = None
            atom.valence_nl = None

        base_N, base_A = [], []
        idx = PZ
        nalfa_index = QZ
        coef_index = SZ
        
        while idx < len(control):
            kz, indicador, nl, mz = control[idx]
            idx += 1

            if not nl[-1].isalpha() or not nl[:-1].isdigit():
                print(f"[WARN] NL inválido: {nl} (Z={Z})")
                continue

            n = int(nl[:-1])
            l = L_MAP.get(nl[-1].upper(), -1)
            if l == -1:
                print(f"[WARN] tipo l desconocido para {nl}")
                continue

            C = []
            #if Z == 103:
                #print(f"[DEBUG] coeficiente {ci} para Z={Z}, orbital={nl}, control={control[idx]} nalfa_index={nalfa_index} length={len(nalfa)} indicador={indicador}, kz={kz}")
            # Leer nueva base si corresponde
            
            base_N, base_A = [], []
            for _ in range(kz):
                try:
                    ni, ai = nalfa[nalfa_index]
                    base_N.append(ni)
                    base_A.append(ai)
                    nalfa_index += 1
                    
                except IndexError:
                    print(f"[ERROR] nalfa insuficiente para Z={Z}, orbital={nl}")
                    break

            # SIEMPRE leer nuevos coeficientes
            for _ in range(kz):
                try:
                    ci = coef[coef_index] * math.sqrt(mz)
                    C.append(ci)
                    coef_index += 1
                    
                except IndexError:
                    print(f"[ERROR] coeficiente faltante para Z={Z}, orbital={nl}")
                    break

            orbital = Orbital(n=n, l=l, mz=mz, N=base_N.copy(), A=base_A.copy(), C=C)
            atom.add_orbital(orbital)
            if indicador == -1:
                nalfa_index -= kz
            if indicador == 0:
                break
            

        data[(Z, Q)] = atom

    return data


def load_system(Q: int, base_path: str) -> Dict[Tuple[int, int], AtomicData]:
    
    suffix = extension_from_Q(Q)  # Determina la extensión del archivo según Z y Q
    if suffix is None:
        raise ValueError(f"Sistema Q={Q} no soportado")

    def p(name): return os.path.join(base_path, f"{name}.{suffix}")

    atoms = read_atoms(p("atoms"))
    control = read_control(p("control"))
    nalfa = read_nalfa(p("nalfa"))
    coef = read_coefficients(p("c"))
    energy = read_energy(p("energy"))
    if suffix not in ["ani", "neu", "cat"]: 
        potioniz = None
        valence = None 
    else:
        potioniz = read_potioniz(p("potioniz"))
        valence = read_valence(p("valen"))
        

    return assemble_atomic_data(
        Q=Q, atoms=atoms,
        control=control, nalfa=nalfa,
        coef=coef, energy=energy,
        potioniz=potioniz, valence=valence
    )

def csv_export_atomic_system(
    atomic_data: Dict[Tuple[int, int], AtomicData],
    output_file: str,
    delimiter: str = ";",
) -> None:
    with open(output_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=delimiter)

        # Cabecera
        writer.writerow([
            "Z", "Q", "N", "Config", "Energy", "Ionization",
            "Valence_n", "Valence_l", "n", "l", "mz",
            "N_list", "A_list", "C_list"
        ])

        # Escribir orbitales por átomo
        for atom in atomic_data.values():
            for orb in atom.orbitals:
                writer.writerow([
                    atom.Z,
                    atom.Q,
                    atom.N,
                    atom.config,
                    str(atom.energy).replace(".", ","),
                    str(atom.ionization).replace(".", ","),
                    atom.valence_nl[0] if atom.ionization is not None else "None",
                    atom.valence_nl[1] if atom.ionization is not None else "None",
                    orb.n,
                    orb.l,
                    orb.mz,
                    "[" + ", ".join(map(str, orb.N)) + "]",
                    "[" + ", ".join(f"{a:.6f}" for a in orb.A) + "]",
                    "[" + ", ".join(f"{c:.12f}" for c in orb.C) + "]"
                ])

    print(f"✅ Exportado correctamente: {output_file} con {len(atomic_data)} átomos")

