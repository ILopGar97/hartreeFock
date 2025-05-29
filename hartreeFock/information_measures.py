from data_structures import *
from typing import Any, List, Dict, Tuple, Union, Optional
from wave_function_constructor import Rnl_monoparticular, rho_monoparticular, Ylm, Gamma_biparticular
import numpy as np
from scipy import integrate
import math
import warnings

epsrel, epsrel_d= 1.49e-10 , 1.49e-6
epsabs, epsabs_d = 1.49e-10, 1.49e-6

def Rnl_normalization(
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
                integral, error = integrate.quad(lambda r: Rnl_monoparticular(r, orbital, spc)**2 * r**2, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
                
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
            ##print(orbital)
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
                    integral, error = integrate.quad(lambda th: 2*np.pi*Ylm(th, orbital, m)**2 * math.sin(th), 0, np.pi, epsabs=epsabs, epsrel=epsrel, limit=1000)
                
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
                integral, error = integrate.quad(lambda th: 2*np.pi*Ylm(th, orbital, m)**2 * math.sin(th), 0, np.pi, epsabs=epsabs, epsrel=epsrel, limit=1000)
            
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
            
            
            integral, error = integrate.quad(lambda r: rho_monoparticular(r, atom, spc)**2 * r**2, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
            
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
        integral, error = integrate.quad(lambda r: rho_monoparticular(r, atoms, space)**2 * r**2, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
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
    include_angular: Optional[bool] = False,   
    two_electron_density: Optional[bool] = False
)-> List[Dict[str, Any]]: 
    """
    Calculate the Shannon entropy, Shannon Lenght and Entropic Power of the atomic system. 
    """
    normalized = True #Dado que las medidas de informacion requieren que la funcion de onda este normalizada, no se puede cambiar.

    if atoms is not None and orbitals is not None:
        raise ValueError("Cannot provide both atoms and orbitals. Choose one.")
    if orbitals is not None and two_electron_density:
        raise ValueError("Cannot provide orbitals and set two_electron_density=True. Choose one.")
        
    
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
                
                if spc not in ["position", "momentum", "my_version"]:
                    raise ValueError("Invalid space. Use 'position' or 'momentum'.") 

                if include_angular:
                    raise ValueError("include_angular=True not implemented yet.")
                
                if atom.Z == 1 and two_electron_density:
                    warnings.warn("The calculation of the two-body probability density calculation for Hydrogen (or systems with only one electron) does not make sense.", UserWarning)
                    continue

                else:
                    if two_electron_density:
                        #print(f"Calculating Shannon Entropy for two-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                        integrando = lambda r1, r2: - 16*np.pi**2 *r1**2 *r2**2 * Gamma_biparticular (r1, r2, atom, spc, True, normalized) * np.log(Gamma_biparticular(r1, r2, atom, spc, True, normalized)) if Gamma_biparticular(r1, r2, atom, spc, True, normalized) > 0 else 0
                        integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}, {'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}])

                    else:
                        #print(f"Calculating Shannon Entropy for mono-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                        integrando = lambda r: - 4*np.pi*r**2*rho_monoparticular(r, atom, spc, True, normalized) * np.log(rho_monoparticular(r, atom, spc, True, normalized)) if rho_monoparticular(r, atom, spc, True, normalized) > 0 else 0
                        integral, error = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
                    results.append({
                        'Z':atom.Z, 
                        'Q':atom.Q,
                        'ShannonEntropy': integral,
                        'ShannonLenght': np.exp(integral), 
                        'EntropicPower':(1/(2*np.pi*np.e))*np.exp(2*integral/3),
                        'integration_error': error,
                        'TwoElectronDensity': two_electron_density,
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
                if two_electron_density:
                    integrando = lambda r1, r2: - 16*np.pi**2 *r1**2 *r2**2 * Gamma_biparticular (r1, r2, atoms, space, True, normalized) * np.log(Gamma_biparticular(r1, r2, atoms, space, True, normalized)) if Gamma_biparticular(r1, r2, atoms, space, True, normalized) > 0 else 0
                    integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}, {'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}])
                else:
                    integrando = lambda r: - 4*np.pi*r**2*rho_monoparticular(r, atoms, space, True, normalized) * np.log(rho_monoparticular(r, atoms, space, True, normalized)) if rho_monoparticular(r, atoms, space, True, normalized) > 0 else 0
                    integral, error = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)

                results.append({
                    'Z':atoms.Z,
                    'Q':atoms.Q,
                    'ShannonEntropy': integral,
                    'ShannonLenght': np.exp(integral), 
                    'EntropicPower':(1/(2*np.pi*np.e))*np.exp(2*integral/3),
                    'integration_error': error,
                    'TwoElectronDensity': two_electron_density,
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
                    integral, error = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
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
                integrando = lambda r: - r**2*(Rnl_monoparticular(r, orbitals, space, normalized)**2) * np.log((1/4*np.pi)*(Rnl_monoparticular(r, orbitals, space, True, normalized))**2) if Rnl_monoparticular(r, orbitals, space, True, normalized)**2 > 0 else 0
                integral, error = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
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

def EntropicMoment(
    atoms: Optional[Union[List[AtomicData], AtomicData]] = None,
    #orbitals: Optional[Union[List[Orbital], Orbital]] = None,
    qs: Optional[Union[List[float], float]] = None,
    space: Optional[Union[List[str], str]] = None,
    include_angular: Optional[bool] = False,
    mono_electron_density: Optional[bool] = True,   
    two_electron_density: Optional[bool] = False, 
    product_density: Optional[bool] = False,
    average_density: Optional[bool] = False
)-> List[Dict[str, Any]]:
    """
    Calculate the frecuency moments of the atomic system. 
    """
    normalized = True

    if atoms is not None:
        if isinstance(atoms, AtomicData):
            if qs <= 0:
                raise ValueError("qs must be positive.")
            atoms = [atoms]
    if qs is not None:
        if isinstance(qs, float):
            qs = [qs]
        if isinstance(qs, list):
            for q in qs:
                if q <= 0:
                    raise ValueError("qs must be positive.")
        else:
            raise TypeError("Expected a number or a list of numbers.")
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
        
        if atom.Z == 1 and two_electron_density:
            warnings.warn("The calculation of the two-body probability density calculation for Hydrogen (or systems with only one electron) does not make sense.", UserWarning)
            continue
        else:
            for  q in qs:
                if space == "momentum" and q < 3/8:
                    raise ValueError("q must be greater than 3/8 in momentum space.")
                if mono_electron_density or product_density:
                    #print(f"Calculating Entropic Moment (q = {q}) for two-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                    integrando = lambda r: 4*np.pi*r**2*rho_monoparticular(r, atom, spc, True, normalized) ** q
                    integral_mono, error_mono = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)
                if two_electron_density:
                    #print(f"Calculating Entropic Moment (q = {q}) for two-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                    integrando = lambda r1, r2: 16*np.pi**2 *r1**2 *r2**2 * Gamma_biparticular (r1, r2, atom, spc, True, normalized)**q
                    integral_two, error_two = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}, {'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}])
                if mono_electron_density:
                    integral_prod = integral_mono ** 2
                    error_prod =  error_mono ** 2
                if average_density:
                    integrando = lambda r1, r2: 16*np.pi**2 *r1**2 *r2**2 * ((Gamma_biparticular (r1, r2, atom, spc, True, normalized) + rho_monoparticular(r1, atom, spc, True, normalized) * rho_monoparticular(r2, atom, spc, True, normalized))/2)**q
                    integral_avg, error_avg = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}, {'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}])

                results.append({
                    'Z':atom.Z, 
                    'Q':atom.Q,
                    'q': q, 
                    'EntropicMoment_mono': integral_mono if mono_electron_density else None,
                    'integration_error_mono': error_mono if mono_electron_density else None,
                    'EntropicMoment_two': integral_two if two_electron_density else None,
                    'integration_error_two': error_two if two_electron_density else None,
                    'EntropicMoment_product': integral_prod if product_density else None,
                    'integration_error_product': error_prod if product_density else None,
                    'EntropicMoment_average': integral_avg if average_density else None,
                    'integration_error_average': error_avg if average_density else None,
                    'space': spc
                })
    return results
    
        
def  RenyiEntropy(
    atoms: Optional[Union[List[AtomicData], AtomicData]] = None,
    #orbitals: Optional[Union[List[Orbital], Orbital]] = None,
    alpha: Optional[Union[List[float], float]]= None,
    space: Optional[Union[List[str], str]] = None,
    include_angular: Optional[bool] = False,   
    two_electron_density: Optional[bool] = False
    
)-> List[Dict[str, Any]]:
    """
    Calculate the Rènyi entropy of the atomic system. 
    """
    normalized = True #Dado que las medidas de informacion requieren que la funcion de onda este normalizada, no se puede cambiar.
    if alpha is not None:
        if isinstance(alpha, float):
            alpha = [alpha]
        if isinstance(alpha, list):
            for a in alpha:
                if a <= 0:
                    raise ValueError("alpha must be positive.")
    else:
        raise TypeError("Expected a number or a list of numbers.")
     
    if atoms is not None:
        if isinstance(atoms, AtomicData):
            atoms = [atoms]
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
                if atom.Z == 1 and two_electron_density:
                    warnings.warn("The calculation of the two-body probability density calculation for Hydrogen (or systems with only one electron) does not make sense.", UserWarning)
                    continue
                else:
                    for  a in alpha:
                        if a == 1: #Calculamos la entropia de Shannon que es el limite de la entropia de Renyi cuando alpha tiende a 1
                            
                            if two_electron_density:
                                #print(f"Calculating Renyi entropy (alpha = {a}) for two-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                                integrando = lambda r1, r2: - 16*np.pi**2 *r1**2 *r2**2 * Gamma_biparticular (r1, r2, atom, spc, True, normalized) * np.log(Gamma_biparticular(r1, r2, atom, spc, True, normalized)) if Gamma_biparticular(r1, r2, atom, spc, True,normalized) > 0 else 0
                                integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}, {'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}])
                            else:
                                #print(f"Calculating Renyi entropy (alpha = {a}) for mono-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                                integrando = lambda r: - 4*np.pi*r**2*rho_monoparticular(r, atom, spc, True, normalized) * np.log(rho_monoparticular(r, atom, spc, True, normalized)) if rho_monoparticular(r, atom, spc, True, normalized) > 0 else 0
                                integral, error = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)

                            results.append({
                                'Z':atom.Z,
                                'Q':atom.Q,
                                'alpha': a,
                                'RenyiEntropy': integral,
                                'integration_error': error,
                                'TwoElectronDensity': two_electron_density,
                                'space': spc
                            })
                        else:
                            if two_electron_density:
                                #print(f"Calculating Renyi entropy (alpha = {a}) for two-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                                integrando = lambda r1, r2: 16*np.pi**2 *r1**2 *r2**2 * Gamma_biparticular (r1, r2, atom, spc, True, normalized)**a
                                integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}, {'epsabs': epsabs_d, 'epsrel': epsrel_d, 'limit': 1000}])

                            else:
                                #print(f"Calculating Renyi entropy (alpha = {a}) for mono-particular probability density of {atom.Name} ({atom.Symbol}: Z={atom.Z}, Q={atom.Q})")
                                integrando = lambda r: 4*np.pi*r**2*rho_monoparticular(r, atom, spc, True, normalized) ** a
                                integral, error = integrate.quad(integrando, 0, np.inf, epsabs=epsabs, epsrel=epsrel, limit=1000)

                            results.append({
                                'Z':atom.Z,
                                'Q':atom.Q,
                                'alpha': a,
                                'RenyiEntropy': (1/(a-1))*np.log(integral),
                                'integration_error': error,
                                'TwoElectronDensity': two_electron_density,
                                'space': spc
                            })
            return results

def JenTsallisDivergence(
    atoms: Optional[Union[List[AtomicData], AtomicData]] = None,
    #orbitals: Optional[Union[List[Orbital], Orbital]] = None,
    alpha: Optional[Union[List[float], float]]= None,
    space: Optional[Union[List[str], str]] = None,
    include_angular: Optional[bool] = False,   
    
)-> List[Dict[str, Any]]:
    """
    Calculate the Jensen-Tsallis divergence of the atomic system (comparation between product of monoparticular densities and biparticular densities)
    """

    normalized = True #Dado que las medidas de informacion requieren que la funcion de onda este normalizada, no se puede cambiar.
    if alpha is not None:
        if isinstance(alpha, float):
            alpha = [alpha]
        if isinstance(alpha, list):
            for a in alpha:
                if a <= 0:
                    raise ValueError("alpha must be positive.")
    else:
        raise TypeError("Expected a number or a list of numbers.")
     
    if atoms is not None:
        if isinstance(atoms, AtomicData):
            atoms = [atoms]
        if isinstance(atoms, list):
            if len(atoms)==0:
                raise ValueError("List of atoms is empty.")
            if space is None: 
                space = ["position"]*len(atoms)
                #warnings.warn("Space not provided. Defaulting to 'position' for all atoms.", UserWarning)
            if isinstance(space, str):
                space = [space]*len(atoms)
                #warnings.warn("Space provided as string. Defaulting to same space for all atoms.", UserWarning)
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
                if atom.Z == 1:
                    warnings.warn("The calculation of the two-body probability density calculation for Hydrogen (or systems with only one electron) does not make sense.", UserWarning)
                    continue
                else:
                    for  a in alpha:
                        if a == 1: #Calculamos JSD que es el limite de la entropia de JRD cuando alpha tiende a 1
                            
                            def gamma(r1, r2): 
                                return Gamma_biparticular(r1, r2, atomic_data=atom, space = spc, spherical_averaged=True, normalized=normalized)
                            def rho1rho2(r1,r2):
                                return rho_monoparticular(r1, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)*rho_monoparticular(r2, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)
                           
                            def integrando(r1, r2): 
                                sumando1 = ((gamma(r1, r2)+rho1rho2(r1,r2))/2)*np.log((gamma(r1, r2)+rho1rho2(r1, r2))/2) if (gamma(r1, r2)+rho1rho2(r1,r2))/2 > 0 else 0
                                sumando2 = gamma(r1, r2)*np.log(gamma(r1, r2)) if gamma(r1, r2) > 0 else 0
                                sumando3 = rho1rho2(r1, r2)*np.log(rho1rho2(r1, r2)) if rho1rho2(r1, r2) > 0 else 0
                                return -(16*np.pi**2 * r1**2 * r2**2)*(sumando1 - 0.5*(sumando2 + sumando3))
                            
                            integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.1, 'epsrel': epsrel_d*0.1, 'limit': 2000}])

                            results.append({
                                'Z':atom.Z,
                                'Q':atom.Q,
                                'alpha': a,
                                'JTD': integral,
                                'integration_error': error,
                                'space': spc
                            })
                        else:
                            def gamma(r1, r2): 
                                return Gamma_biparticular(r1, r2, atomic_data=atom, space = spc, spherical_averaged=True, normalized=normalized)
                            def rho1rho2(r1,r2):
                                return rho_monoparticular(r1, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)*rho_monoparticular(r2, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)
                           
                            def integrando(r1, r2): 
                                sumando1 = ((gamma(r1, r2) + rho1rho2(r1,r2))/2)**a
                                sumando2 = gamma(r1, r2)**a
                                sumando3 = rho1rho2(r1, r2)**a
                                return (16*np.pi**2 * r1**2 * r2**2)*(1/(1-a))*(sumando1 - 0.5*(sumando2 + sumando3))
                            
                            integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                            

                            results.append({
                                'Z':atom.Z,
                                'Q':atom.Q,
                                'alpha': a,
                                'JTD': integral,
                                'integration_error': error,
                                'space': spc
                            })
            return results
        

def JenRenyiDivergence(
    atoms: Optional[Union[List[AtomicData], AtomicData]] = None,
    #orbitals: Optional[Union[List[Orbital], Orbital]] = None,
    alpha: Optional[Union[List[float], float]]= None,
    space: Optional[Union[List[str], str]] = None,
    include_angular: Optional[bool] = False,   
    
)-> List[Dict[str, Any]]:
    """
    Calculate the Jensen-Renyi divergence of the atomic system (comparation between product of monoparticular densities and biparticular densities)
    """

    normalized = True #Dado que las medidas de informacion requieren que la funcion de onda este normalizada, no se puede cambiar.
    if alpha is not None:
        if isinstance(alpha, float):
            alpha = [alpha]
        if isinstance(alpha, list):
            for a in alpha:
                if a <= 0:
                    raise ValueError("alpha must be positive.")
    else:
        raise TypeError("Expected a number or a list of numbers.")
     
    if atoms is not None:
        if isinstance(atoms, AtomicData):
            atoms = [atoms]
        if isinstance(atoms, list):
            if len(atoms)==0:
                raise ValueError("List of atoms is empty.")
            if space is None: 
                space = ["position"]*len(atoms)
                #warnings.warn("Space not provided. Defaulting to 'position' for all atoms.", UserWarning)
            if isinstance(space, str):
                space = [space]*len(atoms)
                #warnings.warn("Space provided as string. Defaulting to same space for all atoms.", UserWarning)
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
                if atom.Z == 1:
                    warnings.warn("The calculation of the two-body probability density calculation for Hydrogen (or systems with only one electron) does not make sense.", UserWarning)
                    continue
                else:
                    for  a in alpha:
                        if a == 1: #Calculamos JSD que es el limite de la entropia de JRD cuando alpha tiende a 1
                            
                            def gamma(r1, r2): 
                                return Gamma_biparticular(r1, r2, atomic_data=atom, space = spc, spherical_averaged=True, normalized=normalized)
                            def rho1rho2(r1,r2):
                                return rho_monoparticular(r1, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)*rho_monoparticular(r2, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)
                           
                            def integrando(r1, r2): 
                                sumando1 = ((gamma(r1, r2)+rho1rho2(r1,r2))/2)*np.log((gamma(r1, r2)+rho1rho2(r1, r2))/2) if (gamma(r1, r2)+rho1rho2(r1,r2))/2 > 0 else 0
                                sumando2 = gamma(r1, r2)*np.log(gamma(r1, r2)) if gamma(r1, r2) > 0 else 0
                                sumando3 = rho1rho2(r1, r2)*np.log(rho1rho2(r1, r2)) if rho1rho2(r1, r2) > 0 else 0
                                return -(16*np.pi**2 * r1**2 * r2**2)*(sumando1 - 0.5*(sumando2 + sumando3))
                            
                            integral, error = integrate.nquad(integrando, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.1, 'epsrel': epsrel_d*0.1, 'limit': 2000}])

                            results.append({
                                'Z':atom.Z,
                                'Q':atom.Q,
                                'alpha': a,
                                'JRD': integral,
                                'integration_error': error,
                                'space': spc
                            })
                        else:
                            def gamma(r1, r2): 
                                return Gamma_biparticular(r1, r2, atomic_data=atom, space = spc, spherical_averaged=True, normalized=normalized)
                            def rho1rho2(r1,r2):
                                return rho_monoparticular(r1, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)*rho_monoparticular(r2, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)
                           
                            def integrales():
                                error = [] 
                                Norm = 16*np.pi**2
                                sumando1 = lambda r1, r2:  Norm * r1**2 * r2**2*((gamma(r1, r2) + rho1rho2(r1,r2))/2)**a
                                integral1, error1 = integrate.nquad(sumando1, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                                error.append(error1)

                                sumando2 = lambda r1, r2: Norm * r1**2 * r2**2*gamma(r1, r2)**a
                                integral2, error2 = integrate.nquad(sumando2, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                                error.append(error2)

                                sumando3 = lambda r1, r2: Norm * r1**2 * r2**2*rho1rho2(r1, r2)**a
                                integral3, error3 = integrate.nquad(sumando3, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                                error.append(error3)

                                return (1/(1-a))*(np.log(integral1) - 0.5*(np.log(integral2) + np.log(integral3))), max(error)
                            
                            integral, error = integrales()
                            

                            results.append({
                                'Z':atom.Z,
                                'Q':atom.Q,
                                'alpha': a,
                                'JRD': integral,
                                'integration_error': error,
                                'space': spc
                            })
            return results        
        

def QSI_alpha(
    atoms: Optional[Union[List[AtomicData], AtomicData]] = None,
    #orbitals: Optional[Union[List[Orbital], Orbital]] = None,
    alpha: Optional[Union[List[float], float]]= None,
    space: Optional[Union[List[str], str]] = None,
    include_angular: Optional[bool] = False,   
    
)-> List[Dict[str, Any]]:
    """
    Calculate the Jensen-Renyi divergence of the atomic system (comparation between product of monoparticular densities and biparticular densities)
    """

    normalized = True #Dado que las medidas de informacion requieren que la funcion de onda este normalizada, no se puede cambiar.
    if alpha is not None:
        if isinstance(alpha, float):
            alpha = [alpha]
        if isinstance(alpha, list):
            for a in alpha:
                if a <= 0:
                    raise ValueError("alpha must be positive.")
    else:
        raise TypeError("Expected a number or a list of numbers.")
     
    if atoms is not None:
        if isinstance(atoms, AtomicData):
            atoms = [atoms]
        if isinstance(atoms, list):
            if len(atoms)==0:
                raise ValueError("List of atoms is empty.")
            if space is None: 
                space = ["position"]*len(atoms)
                #warnings.warn("Space not provided. Defaulting to 'position' for all atoms.", UserWarning)
            if isinstance(space, str):
                space = [space]*len(atoms)
                #warnings.warn("Space provided as string. Defaulting to same space for all atoms.", UserWarning)
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
                if atom.N == 1:
                    warnings.warn("The calculation of the two-body probability density calculation for Hydrogen (or systems with only one electron) does not make sense.", UserWarning)
                    continue
                else:
                    for  a in alpha:
                        def gamma(r1, r2): 
                            return Gamma_biparticular(r1, r2, atomic_data=atom, space = spc, spherical_averaged=True, normalized=normalized)
                        def rho1rho2(r1,r2):
                            return rho_monoparticular(r1, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)*rho_monoparticular(r2, atomic_data=atom, space=spc, spherical_averaged=True, normalized=normalized)
                        
                        def integrales():
                            error = [] 
                            Norm = 16*np.pi**2
                            numerador1 = lambda r1, r2:  Norm * r1**2 * r2**2*rho1rho2(r1, r2)**(a/2) * gamma(r1, r2)**(a/2)
                            integral1, error1 = integrate.nquad(numerador1, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                            error.append(error1)

                            denominador2 = lambda r1, r2: Norm * r1**2 * r2**2*gamma(r1, r2)**a
                            integral2, error2 = integrate.nquad(denominador2, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                            error.append(error2)

                            sumando3 = lambda r1, r2: Norm * r1**2 * r2**2*rho1rho2(r1, r2)**a
                            integral3, error3 = integrate.nquad(sumando3, [[0, np.inf], [0, np.inf]], opts= [{'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}, {'epsabs': epsabs_d*0.10, 'epsrel': epsrel_d*0.10, 'limit': 2000}])
                            error.append(error3)

                            return integral1/np.sqrt(integral2 * integral3), max(error)
                        
                        integral, error = integrales()
                        

                        results.append({
                            'Z':atom.Z,
                            'Q':atom.Q,
                            'alpha': a,
                            'QSI_alpha': integral,
                            'integration_error': error,
                            'space': spc
                        })
            return results