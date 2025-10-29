import os
import math
from typing import List, Dict, Tuple, Union
from data_structures import AtomicData, Orbital
import csv
import json
import re
# ------------------------

L_MAP = {'S': 0, 'P': 1, 'D': 2, 'F': 3}

def load_pubchem_json(path: str) -> Dict[int, dict]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
        pubchem_data = {}

        columns = data["Table"]["Columns"]["Column"]
        for row in data["Table"]["Row"]:
            cells = row["Cell"]
            try:
                Z = int(cells[0])
                pubchem_data[Z] = {}

                for i, column_name in enumerate(columns):
                    value = cells[i] if i < len(cells) and cells[i] != "" else None

                    # Convert numerical fields where needed
                    if column_name in {"AtomicMass", "ElectronAffinity", "MeltingPoint", "BoilingPoint", "Density", "IonizationEnergy", "Electronegativity", "AtomicRadius"}:
                        try:
                            value = float(value) if value is not None else None
                        except ValueError:
                            value = None

                    pubchem_data[Z][column_name] = value

            except (ValueError, IndexError):
                continue

        return pubchem_data

# Prueba de carga
datos_atomicos = load_pubchem_json("PubChemElements_all.json")

elementos = {
    Z: (data["Symbol"], data["Name"])
    for Z, data in datos_atomicos.items()
    if Z > 0 and data.get("Symbol") and data.get("Name")
}
radios_atomicos = {
    Z: data["AtomicRadius"]
    for Z, data in datos_atomicos.items()
    if Z > 0 and data.get("AtomicRadius")
}
electronegatividad = {
    Z: data["Electronegativity"]
    for Z, data in datos_atomicos.items()
    if Z > 0 and data.get("Electronegativity")
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
        # Si no hay electrones, advertencia y valores por defecto
        print(f"[WARN] Ion without electrons detected (Z={Z_original}, Q={Q})")
        return {
            'electronic_configuration': '',
            'full_orbital': "No electrons",
            'hybridization': "None",
            'group': grupo,
            'period': periodo,
            'clasification': "Ion with no electrons"
        }

    # Aquí sí hay electrones disponibles
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

    group_block = datos_atomicos.get(Z_original, {}).get("GroupBlock")

    if group_block:
        clasificacion = group_block  # Usamos directamente la clasificación de PubChem
    else:
        clasificacion = "Unknown"

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
        result = {}
        for line in f:
            t = line.strip().split()
            if len(t) >= 2:
                result[int(t[0])] = float(t[1])
        return result


def read_valence(path: str) -> Dict[int, Tuple[int, int]]:
    with open(path) as f:
        result = {}
        for line in f:
            t = line.strip().split()
            if len(t) >= 3:
                result[int(t[0])] = (int(t[1]), int(t[2]))
        return result


def read_potioniz(path: str) -> Dict[int, float]:
    with open(path) as f:
        result = {}
        for line in f:
            t = line.strip().split()
            if len(t) >= 2:
                result[int(t[0])] = float(t[1])
        return result


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
            period=extra_data["period"],
            clasification=extra_data["clasification"]
        )
        atom.group = int(extra_data["group"]) if extra_data["group"] is not None else None
        atom.energy = energy.get(Z, 0.0)
        atom.atomic_radius = radios_atomicos.get(Z, None)
        atom.electronegativity = electronegatividad.get(Z, None)

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

    print(f"Exported correctly: {output_file} contein {len(atomic_data)} atoms")



def select_atoms_by_Z(
    atomic_data: Dict[Tuple[int, int], 'AtomicData'],
    Z_values: Union[List[int], str],
) -> List['AtomicData']:
    """
        Filters the atoms for their atomic number Z.
        If Z_values is "All", return all the atomicdata.
    """
    if Z_values == "all":
        selected = list(atomic_data.values())
    else:
        selected = [v for v in atomic_data.values() if v.Z in Z_values]

    for atom in selected:
        ordenar_orbitales_clasico(atom)
    
    return selected

def ordenar_orbitales_clasico(atom: AtomicData) -> None:
    # Precalcular el orden clásico de los orbitales
    ORDEN_ORBITALES = {(n, L_MAP[l.upper()]): idx for idx, (n, l) in enumerate(orbitales)}
    atom.orbitals.sort(key=lambda orb: ORDEN_ORBITALES.get((orb.n, orb.l), float('inf')))


MAGNITUDE_KEYS = ['JRD', 'JTD', 'QSI_alpha','EntropicMoment', 'ShannonEntropy', 'ShannonLenght', 
                  'EntropicPower', 'RenyiEntropy']

def safe_float(value, default=0.0):
    try:
        val = float(value)
        if math.isnan(val):
            return default
        return val
    except (ValueError, TypeError):
        return default

def format_table(title, entries):
    if not entries:
        return f"{title}\n(no data)\n\n"

    raw_keys = set().union(*(entry.keys() for entry in entries))
    priority = ['Z', 'Q', 'q', 'alpha'] + MAGNITUDE_KEYS + ['integration_error', 'space']
    keys = [k for k in priority if k in raw_keys] + [k for k in sorted(raw_keys) if k not in priority]

    header = keys
    rows = []

    for item in entries:
        row = []
        for key in keys:
            val = item.get(key, "")
            if isinstance(val, (float, int)) or (isinstance(val, str) and val.replace('.', '', 1).replace('-', '', 1).isdigit()):
                row.append(f"{safe_float(val):.3e}")
            else:
                row.append(str(val))
        rows.append(row)

    col_widths = [max(len(str(cell)) for cell in col) for col in zip(*([header] + rows))]

    table_lines = [title.strip(), "-" * (sum(col_widths) + 3 * len(col_widths))]
    table_lines.append(" | ".join(h.ljust(w) for h, w in zip(header, col_widths)))
    table_lines.append("-" * (sum(col_widths) + 3 * len(col_widths)))

    for row in rows:
        table_lines.append(" | ".join(cell.ljust(w) for cell, w in zip(row, col_widths)))

    table_lines.append("\n")
    return "\n".join(table_lines)

def sanitize_filename(name):
    return re.sub(r'[^\w\-_\. ]', '_', name).replace(' ', '_')

def json_to_multiple_txt(input_path):
    with open(input_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if isinstance(data, dict):
        base_dir = "tables_dict_txt"
        os.makedirs(base_dir, exist_ok=True)
        for title, content in data.items():
            clean_name = sanitize_filename(title.strip()) + ".txt"
            path = os.path.join(base_dir, clean_name)

            if isinstance(content, list):
                table_text = format_table(title, content)
            elif isinstance(content, dict):
                table_text = format_table(title, [content])
            elif isinstance(content, (int, float)):
                table_text = f"{title.strip()}: {content:.6e}\n\n"
            else:
                table_text = f"{title.strip()}: {content}\n\n"

            with open(path, "w", encoding="utf-8") as out:
                out.write(table_text)

    elif isinstance(data, list):
        for mag_key in MAGNITUDE_KEYS:
            filtered_by_mag = [item for item in data if mag_key in item]
            if not filtered_by_mag:
                continue

            # Crear subdirectorio por magnitud
            dir_name = f"tables_{mag_key}_txt"
            os.makedirs(dir_name, exist_ok=True)

            # Detectar parámetro de agrupación
            key_param = None
            if any("q" in item for item in filtered_by_mag):
                key_param = "q"
            elif any("alpha" in item for item in filtered_by_mag):
                key_param = "alpha"

            if key_param:
                values = sorted(set(item[key_param] for item in filtered_by_mag if key_param in item))
                for val in values:
                    subset = [item for item in filtered_by_mag if item.get(key_param) == val]
                    title = f"{mag_key} - {key_param}={val}"
                    filename = sanitize_filename(f"{mag_key}_{key_param}_{val:.2f}.txt")
                    path = os.path.join(dir_name, filename)
                    table_text = format_table(title, subset)
                    with open(path, "w", encoding="utf-8") as out:
                        out.write(table_text)
            else:
                # No hay q o alpha
                title = f"{mag_key}"
                filename = sanitize_filename(f"{mag_key}.txt")
                table_text = format_table(title, filtered_by_mag)
                path = os.path.join(dir_name, filename)
                with open(path, "w", encoding="utf-8") as out:
                    out.write(table_text)

    else:
        raise ValueError("El JSON no tiene un formato compatible. Debe contener un dict o una lista en la raíz.")

