from data_structures import *
from typing import List, Dict, Tuple, Union, Optional
from wave_function_constructor import Rnl_monoparticular, Ylm
import numpy as np
from scipy import integrate
import math


def radial_normalization(
    atoms: Union[List['AtomicData'], 'AtomicData'],
    space: Optional[Union[List[str], str]] = None
) -> List[Dict[str, Any]]:
    """
    Calculate the normalization constant for the radial wave function.
    The normalization constant is defined as:
    """
    if isinstance(atoms,list):
        if len(atoms)==0:
            raise ValueError("List of atoms is empty.")
        if space is None: 
            space = ["position"]*len(atoms)
            raise Warning("Space not provided. Defaulting to 'position' for all atoms.")
        if len(atoms) != len(space):
            raise ValueError("Length of atoms and space lists must be the same.")
        
        results = []
        for atom,spc in zip(atoms, space):
            if not isinstance(atom, AtomicData):
                raise TypeError("Expected an instance of AtomicData.")
            if spc not in ["position", "momentum"]:
                raise ValueError("Invalid space. Use 'position' or 'momentum'.")
            for orbital in atom.orbitals:
                integral, error = integrate.quad(lambda r: Rnl_monoparticular(r, atom.orbital, spc)**2 * r**2, 0, np.inf)
                if integral <= 0:
                    raise ValueError("Integral is not positive.")
                results.append({
                    'Z': atom.Z, 
                    'Q': atom.Q, 
                    'n': orbital.n,
                    'l': orbital.l,
                    'integral': integral,
                    'error': error,
                    'space': spc
                    })
        return results
    
    elif isinstance(atoms, AtomicData):
        if space is None: 
            space = "position"
            raise Warning("Space not provided. Defaulting to 'position'.")
        if not isinstance(atoms, AtomicData):
            raise TypeError("Expected an instance of AtomicData.")
        if space not in ["position", "momentum"]:
            raise ValueError("Invalid space. Use 'position' or 'momentum'.")
        
        results = []
        for orbital in atoms.orbitals:
            integral, error = integrate.quad(lambda r: Rnl_monoparticular(r, atoms.orbital, space)**2 * r**2, 0, np.inf)
            if integral <= 0:
                raise ValueError("Integral is not positive.")
            results.append({
                'Z': atoms.Z, 
                'Q': atoms.Q, 
                'n': orbital.n,
                'l': orbital.l,
                'integral': integral,
                'error': error,
                'space': space
                })
        return results
    else:
        raise TypeError("Expected an instance of AtomicData or a list of AtomicData.")    
    
        
        
    