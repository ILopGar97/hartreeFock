# hf_parser/data_structures.py

from dataclasses import dataclass, field
from typing import List, Tuple, Optional

# ---------- ESTRUCTURAS DE DATOS ----------

@dataclass
class Orbital:
    n: int
    l: int
    mz: int
    N: List[int]
    A: List[float]
    C: List[float]

@dataclass
class AtomicData:
    Name: str
    Symbol: str
    Z: int
    Q: int
    N : int
    electronic_configuration: str
    full_orbital: str
    hybridization: str
    config: str
    period: int
    clasification: str
    orbitals: List[Orbital] = field(default_factory=list)
    energy: float = 0.0
    group: Optional[int] = None
    ionization: Optional[float] = None
    valence_nl: Optional[float] = None

    def add_orbital(self, orbital: Orbital):
        self.orbitals.append(orbital)