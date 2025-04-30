import os
import math
from typing import List, Dict, Tuple, Union
from data_structures import AtomicData, Orbital
import csv
# ------------------------

L_MAP = {'S': 0, 'P': 1, 'D': 2, 'F': 3}

elementos = {
    1: ("H", "Hydrogen"),
    2: ("He", "Helium"),
    3: ("Li", "Lithium"),
    4: ("Be", "Beryllium"),
    5: ("B", "Boron"),
    6: ("C", "Carbon"),
    7: ("N", "Nitrogen"),
    8: ("O", "Oxygen"),
    9: ("F", "Fluorine"),
    10: ("Ne", "Neon"),
    11: ("Na", "Sodium"),
    12: ("Mg", "Magnesium"),
    13: ("Al", "Aluminum"),
    14: ("Si", "Silicon"),
    15: ("P", "Phosphorus"),
    16: ("S", "Sulfur"),
    17: ("Cl", "Chlorine"),
    18: ("Ar", "Argon"),
    19: ("K", "Potassium"),
    20: ("Ca", "Calcium"),
    21: ("Sc", "Scandium"),
    22: ("Ti", "Titanium"),
    23: ("V", "Vanadium"),
    24: ("Cr", "Chromium"),
    25: ("Mn", "Manganese"),
    26: ("Fe", "Iron"),
    27: ("Co", "Cobalt"),
    28: ("Ni", "Nickel"),
    29: ("Cu", "Copper"),
    30: ("Zn", "Zinc"),
    31: ("Ga", "Gallium"),
    32: ("Ge", "Germanium"),
    33: ("As", "Arsenic"),
    34: ("Se", "Selenium"),
    35: ("Br", "Bromine"),
    36: ("Kr", "Krypton"),
    37: ("Rb", "Rubidium"),
    38: ("Sr", "Strontium"),
    39: ("Y", "Yttrium"),
    40: ("Zr", "Zirconium"),
    41: ("Nb", "Niobium"),
    42: ("Mo", "Molybdenum"),
    43: ("Tc", "Technetium"),
    44: ("Ru", "Ruthenium"),
    45: ("Rh", "Rhodium"),
    46: ("Pd", "Palladium"),
    47: ("Ag", "Silver"),
    48: ("Cd", "Cadmium"),
    49: ("In", "Indium"),
    50: ("Sn", "Tin"),
    51: ("Sb", "Antimony"),
    52: ("Te", "Tellurium"),
    53: ("I", "Iodine"),
    54: ("Xe", "Xenon"),
    55: ("Cs", "Cesium"),
    56: ("Ba", "Barium"),
    57: ("La", "Lanthanum"),
    58: ("Ce", "Cerium"),
    59: ("Pr", "Praseodymium"),
    60: ("Nd", "Neodymium"),
    61: ("Pm", "Promethium"),
    62: ("Sm", "Samarium"),
    63: ("Eu", "Europium"),
    64: ("Gd", "Gadolinium"),
    65: ("Tb", "Terbium"),
    66: ("Dy", "Dysprosium"),
    67: ("Ho", "Holmium"),
    68: ("Er", "Erbium"),
    69: ("Tm", "Thulium"),
    70: ("Yb", "Ytterbium"),
    71: ("Lu", "Lutetium"),
    72: ("Hf", "Hafnium"),
    73: ("Ta", "Tantalum"),
    74: ("W", "Tungsten"),
    75: ("Re", "Rhenium"),
    76: ("Os", "Osmium"),
    77: ("Ir", "Iridium"),
    78: ("Pt", "Platinum"),
    79: ("Au", "Gold"),
    80: ("Hg", "Mercury"),
    81: ("Tl", "Thallium"),
    82: ("Pb", "Lead"),
    83: ("Bi", "Bismuth"),
    84: ("Po", "Polonium"),
    85: ("At", "Astatine"),
    86: ("Rn", "Radon"),
    87: ("Fr", "Francium"),
    88: ("Ra", "Radium"),
    89: ("Ac", "Actinium"),
    90: ("Th", "Thorium"),
    91: ("Pa", "Protactinium"),
    92: ("U", "Uranium"),
    93: ("Np", "Neptunium"),
    94: ("Pu", "Plutonium"),
    95: ("Am", "Americium"),
    96: ("Cm", "Curium"),
    97: ("Bk", "Berkelium"),
    98: ("Cf", "Californium"),
    99: ("Es", "Einsteinium"),
    100: ("Fm", "Fermium"),
    101: ("Md", "Mendelevium"),
    102: ("No", "Nobelium"),
    103: ("Lr", "Lawrencium"),
}


# Definimos los orbitales en orden de llenado
orbitales = [
    (1, 's'), (2, 's'),
    (2, 'p'), (3, 's'),
    (3, 'p'), (4, 's'),
    (3, 'd'), (4, 'p'),
    (5, 's'), (4, 'd'), 
    (5, 'p'), (6, 's'),
    (4, 'f'), (5, 'd'), (6, 'p'),
    (7, 's'), (5, 'f'), (6, 'd'), (7, 'p')
]

# Máxima capacidad de cada subnivel
capacidad = {
    's': 2,
    'p': 6,
    'd': 10,
    'f': 14
}

# Clasificación de elementos basada en Z
no_metales = [1, 6, 7, 8, 9, 15, 16, 17, 34, 35, 53]
metaloides = [5, 14, 32, 33, 51, 52]
metales_postransicionales = [13, 31, 49, 50, 81, 82, 83, 84, 85]
lantanidos = list(range(57, 72))
actinidos = list(range(89, 104))

# Datos básicos de grupo y periodo (solo una parte, podemos expandir)
datos_basicos = {
    1: (1, 1), 2: (18, 1),
    3: (1, 2), 4: (2, 2), 5: (13, 2), 6: (14, 2), 7: (15, 2), 8: (16, 2), 9: (17, 2), 10: (18, 2),
    11: (1, 3), 12: (2, 3), 13: (13, 3), 14: (14, 3), 15: (15, 3), 16: (16, 3), 17: (17, 3), 18: (18, 3),
    19: (1, 4), 20: (2, 4), 21: (3, 4), 22: (4, 4), 23: (5, 4), 24: (6, 4), 25: (7, 4), 26: (8, 4), 27: (9, 4), 28: (10, 4),
    29: (11, 4), 30: (12, 4), 31: (13, 4), 32: (14, 4), 33: (15, 4), 34: (16, 4), 35: (17, 4), 36: (18, 4),
    37: (1, 5), 38: (2, 5), 39: (3, 5), 40: (4, 5), 41: (5, 5), 42: (6, 5), 43: (7, 5), 44: (8, 5), 45: (9, 5), 46: (10, 5),
    47: (11, 5), 48: (12, 5), 49: (13, 5), 50: (14, 5), 51: (15, 5), 52: (16, 5), 53: (17, 5), 54: (18, 5),
    55: (1, 6), 56: (2, 6), 57: (3, 6), 58: (None, 6), 59: (None, 6), 60: (None, 6), 61: (None, 6), 62: (None, 6),
    63: (None, 6), 64: (None, 6), 65: (None, 6), 66: (None, 6), 67: (None, 6), 68: (None, 6), 69: (None, 6), 70: (None, 6),
    71: (None, 6), 72: (4, 6), 73: (5, 6), 74: (6, 6), 75: (7, 6), 76: (8, 6), 77: (9, 6), 78: (10, 6),
    79: (11, 6), 80: (12, 6), 81: (13, 6), 82: (14, 6), 83: (15, 6), 84: (16, 6), 85: (17, 6), 86: (18, 6),
    87: (1, 7), 88: (2, 7), 89: (3, 7), 90: (None, 7), 91: (None, 7), 92: (None, 7), 93: (None, 7), 94: (None, 7),
    95: (None, 7), 96: (None, 7), 97: (None, 7), 98: (None, 7), 99: (None, 7), 100: (None, 7), 101: (None, 7),
    102: (None, 7), 103: (None, 7),
}


def configuracion_electronica_completa(Z, Q):
    Z_original = Z
    N = Z - Q
    N_original = N
    configuracion = []
    electrones_por_orbital = []

    for n, subnivel in orbitales:
        max_electrones = capacidad[subnivel]
        if N >= max_electrones:
            configuracion.append(f"{n}{subnivel}^{max_electrones}")
            electrones_por_orbital.append((n, subnivel, max_electrones))
            N -= max_electrones
        else:
            if N > 0:
                configuracion.append(f"{n}{subnivel}^{N}")
                electrones_por_orbital.append((n, subnivel, N))
                N = 0
            break

    configuracion_str = ' '.join(configuracion)

    grupo, periodo = datos_basicos.get(Z_original, (None, None))

    if not electrones_por_orbital:
        # 🛑 Si no hay electrones, advertencia y valores por defecto
        print(f"[WARN] Ion without electrons detected (Z={Z_original}, Q={Q})")
        return {
            'electronic_configuration': '',
            'full_orbital': "No electrons",
            'hybridization': "None",
            'group': grupo,
            'period': periodo,
            'clasification': "Ion with no electrons"
        }

    # ✅ Aquí sí hay electrones disponibles
    ultima_n, ultima_subnivel, ultima_cantidad = electrones_por_orbital[-1]
    llena = ultima_cantidad == capacidad[ultima_subnivel]
    semillena = ultima_cantidad == capacidad[ultima_subnivel] // 2

    if N_original in [2, 10, 18, 36, 54, 86]:
        tipo_capa = "Noble gas"
    elif llena:
        tipo_capa = "Full orbital"
    elif semillena:
        tipo_capa = "Half orbital"
    else:
        tipo_capa = "No special"

    if ultima_subnivel == 's':
        hibridacion = 'Likely without hybridisation'
    elif ultima_subnivel == 'p':
        hibridacion = 'sp' if ultima_cantidad <= 2 else 'sp2' if ultima_cantidad <= 4 else 'sp3'
    elif ultima_subnivel == 'd':
        hibridacion = 'sp3d'
    elif ultima_subnivel == 'f':
        hibridacion = 'Likely without hybridisation'

    if Z_original in no_metales:
        clasificacion = "Non metal"
    elif Z_original in metaloides:
        clasificacion = "Metaloid"
    elif Z_original in metales_postransicionales:
        clasificacion = "Post-transitional metal"
    elif Z_original in lantanidos:
        clasificacion = "Lanthanide"
    elif Z_original in actinidos:
        clasificacion = "Actinide"
    elif grupo == 1 and Z_original != 1:
        clasificacion = "Alkali metal"
    elif grupo == 2:
        clasificacion = "Alkaline earth metal"
    elif Z_original in [2, 10, 18, 36, 54, 86]:
        clasificacion = "Noble gas"
    else:
        clasificacion = "Transition metal"

    return {
        'electronic_configuration': configuracion_str,
        'full_orbital': tipo_capa,
        'hybridization': hibridacion,
        'group': grupo,
        'period': periodo,
        'clasification': clasificacion
    }



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
        extra_data = configuracion_electronica_completa(Z, Q)
        atom = AtomicData(
            Z=Z,
            Q=Q,
            N=Z-Q,
            config=config,
            Name=elementos[Z][1],
            Symbol=elementos[Z][0],
            electronic_configuration=extra_data["electronic_configuration"],
            full_orbital=extra_data["full_orbital"],
            hybridization=extra_data["hybridization"],
            period=extra_data["period"],
            clasification=extra_data["clasification"]
        )
        atom.group = int(extra_data["group"]) if extra_data["group"] is not None else None
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
            "Element name", "Symbol", "clasificacion", "Possible hybridization", "Z", "Q", "N", "Complete electronic configuration", "Special configuration", "Config", "Energy", "Ionization",
            "Valence_n", "Valence_l", "n", "l", "mz",
            "N_list", "A_list", "C_list"
        ])

        # Escribir orbitales por átomo
        for atom in atomic_data:
            for orb in atom.orbitals:
                writer.writerow([
                    atom.Name,
                    atom.Symbol,
                    atom.clasification,
                    atom.hybridization,
                    atom.Z,
                    atom.Q,
                    atom.N,
                    atom.electronic_configuration,
                    atom.full_orbital,
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

    print(f"✅ Exported correctly: {output_file} contein {len(atomic_data)} atoms")



def select_atoms_by_Z(
    atomic_data: Dict[Tuple[int, int], 'AtomicData'],
    Z_values: Union[List[int], str],
) -> List['AtomicData']:
    """
        Filters the atoms for their atomic number Z.
        If Z_values is "All", return all the atomicdata.
    """
    if Z_values == "all":
        return list(atomic_data.values())
    else:
        return [v for v in atomic_data.values() if v.Z in Z_values]
