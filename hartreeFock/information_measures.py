from data_structures import *
from typing import Any, List, Dict, Tuple, Union, Optional
from wave_function_constructor import Rnl_monoparticular, rho_monoparticular, Ylm
import numpy as np
from scipy import integrate
import math
import warnings


def radial_normalization(
    atoms: Union[List['AtomicData'], 'AtomicData'],
    space: Optional[Union[List[str], str]] = None
) -> List[Dict[str, Any]]:
    """
    Calculate the normalization constant for the radial wave function.
    
    atoms: AtomicData or list of AtomicData, the atomic data object(s) containing the parameters
    space: str or list of str, either "position" or "momentum"

    RETURN: list of dictionaries, each containing the normalization constant and other information for each atom
    """
    if isinstance(atoms,list):
        """
            Check instantces of input
        """
        if len(atoms)==0:
            raise ValueError("List of atoms is empty.")
        if space is None: 
            space = ["position"]*len(atoms)
            warnings.warn("Space not provided. Defaulting to 'position' for all atoms.", UserWarning)
        if isinstance(space, str):
            space = [space]*len(atoms)
            warnings.warn("Space provided as string. Defaulting to same space for all atoms.", UserWarning)
        if len(atoms) != len(space):
            raise ValueError("Length of atoms and space lists must be the same.")
        """
        -----------------------------------------------------------------------------------
        """
        results = []
        for atom,spc in zip(atoms, space):
            if not isinstance(atom, AtomicData):
                raise TypeError("Expected an instance of AtomicData.")
            
            if spc not in ["position", "momentum"]:
                raise ValueError("Invalid space. Use 'position' or 'momentum'.")
            
            for orbital in atom.orbitals:
                integral, error = integrate.quad(lambda r: Rnl_monoparticular(r, orbital, spc)**2 * r**2, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
                
                results.append({
                    'Z': atom.Z, 
                    'Q': atom.Q, 
                    'n': orbital.n,
                    'l': orbital.l,
                    'integral': integral,
                    'integration_error': error,
                    'normalization_error': abs(1.0 - integral),
                    'space': spc,
                    'part': "radial"
                    })
        return results
    
    elif isinstance(atoms, AtomicData):
        """
        Check instantces of input
        """
        if space is None: 
            space = "position"
            warnings.warn("Space not provided. Defaulting to 'position'.", UserWarning)
        if isinstance(space, list):
            if len(space) != 1:
                raise ValueError("Length of space list must be 1 for a single atom.")
            else:
                space = space[0]
        if space not in ["position", "momentum"]:
            raise ValueError("Invalid space. Use 'position' or 'momentum'.")
        """
        -----------------------------------------------------------------------------------
        """
        results = []
        for orbital in atoms.orbitals:
            print(orbital)
            integral, error = integrate.quad(lambda r: Rnl_monoparticular(r, orbital, space)**2 * r**2, 0, np.inf)
            
            results.append({
                'Z': atoms.Z, 
                'Q': atoms.Q, 
                'n': orbital.n,
                'l': orbital.l,
                'integral': integral,
                'integration_error': error,
                'normalization_error': abs(1.0 - integral),
                'space': space,
                'part': "radial"
                })
        return results
    else:
        raise TypeError("Expected an instance of AtomicData or a list of AtomicData.")    

def angular_normalization(
    atoms: Union[List['AtomicData'], 'AtomicData']
) -> List[Dict[str, Any]]:
    """
    Calculate the normalization constant for the radial wave function.
    The normalization constant is defined as:
    """
    if isinstance(atoms,list):
        """
            Check instantces of input
        """
        if len(atoms)==0:
            raise ValueError("List of atoms is empty.")
        """
        -----------------------------------------------------------------------------------
        """
        results = []
        for atom in atoms:
            if not isinstance(atom, AtomicData):
                raise TypeError("Expected an instance of AtomicData.")
            
            for orbital in atom.orbitals:
                l = orbital.l
                for m in range(-l, l+1):
                    integral, error = integrate.quad(lambda th: 2*np.pi*Ylm(th, orbital, m)**2 * math.sin(th), 0, np.pi, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
                
                    results.append({
                        'Z': atom.Z, 
                        'Q': atom.Q, 
                        'n': orbital.n,
                        'l': orbital.l,
                        'm': m,
                        'integral': integral,
                        'integration_error': error,
                        'normalization_error': abs(1.0 - integral),
                        'part': "angular"
                        })
        return results
    
    elif isinstance(atoms, AtomicData):
        results = []
        for orbital in atoms.orbitals:
            l = orbital.l
            for m in range(-l, l+1):
                integral, error = integrate.quad(lambda th: 2*np.pi*Ylm(th, orbital, m)**2 * math.sin(th), 0, np.pi, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
            
                results.append({
                    'Z': atom.Z, 
                    'Q': atom.Q, 
                    'n': orbital.n,
                    'l': orbital.l,
                    'm': m,
                    'integral': integral,
                    'integration_error': error,
                    'normalization_error': abs(1.0 - integral),
                    'part': "angular"
                    })
        return results
    else:
        raise TypeError("Expected an instance of AtomicData or a list of AtomicData.")
    
def rho_normalization(
        atoms: Union[List['AtomicData'], 'AtomicData'],
        space: Optional[Union[List[str], str]] = None
) -> List[Dict[str, Any]]:
    """
    Calculate the normalization of the density of the atomic system."""

    if isinstance(atoms,list):
        """
            Check instantces of input
        """
        if len(atoms)==0:
            raise ValueError("List of atoms is empty.")
        if space is None: 
            space = ["position"]*len(atoms)
            warnings.warn("Space not provided. Defaulting to 'position' for all atoms.", UserWarning)
        if isinstance(space, str):
            space = [space]*len(atoms)
            warnings.warn("Space provided as string. Defaulting to same space for all atoms.", UserWarning)
        if len(atoms) != len(space):
            raise ValueError("Length of atoms and space lists must be the same.")
        """
        -----------------------------------------------------------------------------------
        """
        results = []
        for atom,spc in zip(atoms, space):
            if not isinstance(atom, AtomicData):
                raise TypeError("Expected an instance of AtomicData.")
            
            if spc not in ["position", "momentum"]:
                raise ValueError("Invalid space. Use 'position' or 'momentum'.")
            
            
            integral, error = integrate.quad(lambda r: rho_monoparticular(r, atom, spc)**2 * r**2, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
            
            results.append({
                'Z': atom.Z, 
                'Q': atom.Q, 
                'integral': integral,
                'integration_error': error,
                'normalization_error': abs(1.0 - integral),
                'space': spc,
                'part': "rho"
                })
        return results
    
    elif isinstance(atoms, AtomicData):
        """
        Check instantces of input
        """
        if space is None: 
            space = "position"
            warnings.warn("Space not provided. Defaulting to 'position'.", UserWarning)
        if isinstance(space, list):
            if len(space) != 1:
                raise ValueError("Length of space list must be 1 for a single atom.")
            else:
                space = space[0]
        if space not in ["position", "momentum"]:
            raise ValueError("Invalid space. Use 'position' or 'momentum'.")
        """
        -----------------------------------------------------------------------------------
        """
        results = []
        integral, error = integrate.quad(lambda r: rho_monoparticular(r, atoms, space)**2 * r**2, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
        results.append({
            'Z': atoms.Z, 
            'Q': atoms.Q, 
            'integral': integral,
            'integration_error': error,
            'normalization_error': abs(1.0 - integral),
            'space': space,
            'part': "rho"
            })
        return results
    else:
        raise TypeError("Expected an instance of AtomicData or a list of AtomicData.")
 
def ShannonEntropy(
    atoms: Optional[Union[List[AtomicData], AtomicData]] = None,
    orbitals: Optional[Union[List[Orbital], Orbital]] = None,
    space: Optional[Union[List[str], str]] = None,
    include_angular: Optional[bool] = False   
)-> List[Dict[str, Any]]: 
    """
    Calculate the Shannon entropy, Shannon Lenght and Entropic Power of the atomic system. 
    """
    normalized = True #Dado que las medidas de informacion requieren que la funcion de onda este normalizada, no se puede cambiar.

    if atoms is not None and orbitals is not None:
        raise ValueError("Cannot provide both atoms and orbitals. Choose one.")
        
    
    if atoms is not None:
        if isinstance(atoms, list):
            if len(atoms)==0:
                raise ValueError("List of atoms is empty.")
            if space is None: 
                space = ["position"]*len(atoms)
                warnings.warn("Space not provided. Defaulting to 'position' for all atoms.", UserWarning)
            if isinstance(space, str):
                space = [space]*len(atoms)
                warnings.warn("Space provided as string. Defaulting to same space for all atoms.", UserWarning)
            if len(atoms) != len(space):
                raise ValueError("Length of atoms and space lists must be the same.")

            """
            -----------------------------------------------------------------------------------
            """
            results = []
            for atom,spc in zip(atoms, space):
                if not isinstance(atom, AtomicData):
                    raise TypeError("Expected an instance of AtomicData.")
                
                if spc not in ["position", "momentum"]:
                    raise ValueError("Invalid space. Use 'position' or 'momentum'.") 

                if include_angular:
                    raise ValueError("include_angular=True not implemented yet.")
                else:   # Se añade el termino 4pi para la normalización del promedio esférico
                    integrando = lambda r: - 4*np.pi*r**2*rho_monoparticular(r, atom, spc, True, normalized) * np.log(rho_monoparticular(r, atom, spc, True, normalized)) if rho_monoparticular(r, atom, spc, True, normalized) > 0 else 0
                    integral, error = integrate.quad(integrando, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
                    results.append({
                        'Z':atom.Z, 
                        'Q':atom.Q,
                        'ShannonEntropy': integral,
                        'ShannonLenght': np.exp(integral), 
                        'EntropicPower':(1/(2*np.pi*np.e))*np.exp(2*integral/3),
                        'integration_error': error,
                        'space': spc
                    })
            return results
        if isinstance(atoms, AtomicData):
            """
            Check instantces of input
            """
            if space is None: 
                space = "position"
                warnings.warn("Space not provided. Defaulting to 'position'.", UserWarning)
            if isinstance(space, list):
                if len(space) != 1:
                    raise ValueError("Length of space list must be 1 for a single atom.")
                else:
                    space = space[0]
            if space not in ["position", "momentum"]:
                raise ValueError("Invalid space. Use 'position' or 'momentum'.")
            """
            -----------------------------------------------------------------------------------
            """
            results = []
            if include_angular:
                raise ValueError("include_angular=True not implemented yet.")
            else:   
                integrando = lambda r: - 4*np.pi*r**2*rho_monoparticular(r, atoms, space, True, normalized) * np.log(rho_monoparticular(r, atoms, space, True, normalized)) if rho_monoparticular(r, atoms, space, True, normalized) > 0 else 0
                integral, error = integrate.quad(integrando, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
                results.append({
                    'Z':atoms.Z,
                    'Q':atoms.Q,
                    'ShannonEntropy': integral,
                    'ShannonLenght': np.exp(integral), 
                    'EntropicPower':(1/(2*np.pi*np.e))*np.exp(2*integral/3),
                    'integration_error': error,
                    'space': space
                })        
            return results
        
    elif orbitals is not None:
        if isinstance(orbitals, list):
            if len(orbitals)==0:
                raise ValueError("List of orbitals is empty.")
            if space is None: 
                space = ["position"]*len(orbitals)
                warnings.warn("Space not provided. Defaulting to 'position' for all orbitals.", UserWarning)
            if isinstance(space, str):
                space = [space]*len(orbitals)
                warnings.warn("Space provided as string. Defaulting to same space for all orbitals.", UserWarning)
            if len(orbitals) != len(space):
                raise ValueError("Length of orbitals and space lists must be the same.")
            """
            -----------------------------------------------------------------------------------
            """
            results = []
            for orbital,spc in zip(orbitals, space):
                if not isinstance(orbital, Orbital):
                    raise TypeError("Expected an instance of Orbital.")
                
                if spc not in ["position", "momentum"]:
                    raise ValueError("Invalid space. Use 'position' or 'momentum'.") 

                if include_angular:
                    raise ValueError("include_angular=True not implemented yet.")
                else:   #aquí no se añade el término 4pi para la normalización del promedio esférico dado que cancela con la integracion angular, pero si se añade dentro del argumento del logaritmo
                    integrando = lambda r: - r**2*(Rnl_monoparticular(r, orbital, spc, True, normalized)**2) * np.log((1/(4*np.pi))*(Rnl_monoparticular(r, orbital, spc, True, normalized))**2) if Rnl_monoparticular(r, orbital, spc, True, normalized)**2 > 0 else 0
                    integral, error = integrate.quad(integrando, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
                    results.append({
                        'n':orbital.n,
                        'l':orbital.l,
                        'ShannonEntropy': integral,
                        'ShannonLenght': np.exp(integral), 
                        'EntropicPower':(1/(2*np.pi*np.e))*np.exp(2*integral/3),
                        'integration_error': error,
                        'space': spc
                    })
            return results
        if isinstance(orbitals, Orbital):
            """
            Check instantces of input
            """
            if space is None: 
                space = "position"
                warnings.warn("Space not provided. Defaulting to 'position'.", UserWarning)
            if isinstance(space, list):
                if len(space) != 1:
                    raise ValueError("Length of space list must be 1 for a single atom or orbital.")
                else:
                    space = space[0]
            if space not in ["position", "momentum"]:
                raise ValueError("Invalid space. Use 'position' or 'momentum'.")
            """
            -----------------------------------------------------------------------------------
            """
            results = []
            if include_angular:
                raise ValueError("include_angular=True not implemented yet.")
            else:   
                integrando = lambda r: - r**2(Rnl_monoparticular(r, orbitals, space, normalized)**2) * np.log((1/4*np.pi)*(Rnl_monoparticular(r, orbitals, space, True, normalized))**2) if Rnl_monoparticular(r, orbitals, space, True, normalized)**2 > 0 else 0
                integral, error = integrate.quad(integrando, 0, np.inf, epsabs=1.49e-14, epsrel=1.49e-14, limit=1000)
                results.append({
                    'n':orbitals.n,
                    'l':orbitals.l,
                    'ShannonEntropy': integral,
                    'ShannonLenght': np.exp(integral), 
                    'EntropicPower':(1/(2*np.pi*np.e))*np.exp(2*integral/3),
                    'integration_error': error,
                    'space': space
                })
            return results
    else:
        raise ValueError("Must provide either atoms or orbitals.")
    
