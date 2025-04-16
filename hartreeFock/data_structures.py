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
    Z: int
    Q: int
    N : int
    config: str
    orbitals: List[Orbital] = field(default_factory=list)
    energy: float = 0.0
    ionization: Optional[float] = None
    valence_nl: Optional[float] = None

    def add_orbital(self, orbital: Orbital):
        self.orbitals.append(orbital)