import numpy as np
import math
from data_structures import *
from scipy.special import lpmv 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import warnings


def Rnl_monoparticular(r: float, Orbital: Orbital, space: str, normalized = True) -> float:
    """
    Calculate the radial wave function for a given orbital.
    The radial wave function is defined as:

    r: float, the distance from the nucleus
    Orbital: Orbital, the orbital object containing the parameters
    space: str, either "position" or "momentum"
    normalized: bool, true if the function is normalized to 1, false if is normlized to number of electrons in the atom

    RETURN: float, the value of the radial wave function at distance r
    """
    A_list = Orbital.A
    C_list = Orbital.C
    N_list = Orbital.N
    result = 0.0
    mz = Orbital.mz
    l = Orbital.l

    if space == "position":
        
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):
            K_i = np.power(2*A_i, N_i + 0.5)/np.sqrt(math.factorial(2*N_i))
            if r == 0 and N_i == 1:
                result += C_i * K_i * np.exp(-A_i * r)
            elif r == 0 and N_i != 1:
                result += 0.0
            else:
                result += C_i *K_i* (r ** (N_i - 1)) * np.exp(-A_i * r) 

        if normalized:
            return result / np.sqrt(mz)
        else:
            return result 
            
    elif space == "momentum":
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):
            x = r / A_i
            y = 1 + x**2
            prefactor = C_i * x**l * 4**N_i  / np.sqrt(math.factorial(2 * N_i) * A_i**3)
            sfk = math.factorial(N_i) / y**(N_i + 1)
            sumk = 0.0
            for k in range(0, int((N_i - l) // 2) + 1):
                sumk += sfk
                numerator = (N_i - l - 2 * k) * (N_i - l - 2 * k - 1) * y
                denominator = (N_i - k) * (k + 1) * 4
                sfk *= -numerator / denominator
            result += prefactor * sumk
        result = result*2/np.sqrt(np.pi) # Este factor no aparece estar en el código Fortran original dado que solo  se calcula RHO y no la función de onda
        return  result /(np.sqrt( mz)) if normalized else result 
    
       
    elif space == "my_version":
        
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):

            Norm = C_i*np.power(2*A_i, N_i + 0.5)/(np.sqrt(math.factorial(2*N_i))) 

            Const1 = math.factorial(N_i)*math.sqrt(2/math.pi)/np.power(A_i, N_i + 2)

            Const2 = 1.0/(1+(r/A_i)**2)**(N_i + 1)

            sum = 0
            for k in range(0, int(N_i//2 + 1)):
                sum += math.comb(N_i+1, 2*k + 1)*np.power(-1, k)*np.power(r/A_i, 2*k)
            result += Norm*Const1*Const2*sum 

        if normalized:
            return result / np.sqrt(mz)
        else:
            return result
    else:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")
    

   
       
    
        
    
    
def Ylm(theta: float, orbital: Orbital,  m: int) -> float: #se omite la fase compleja
    """
    Calculate the spherical harmonics Y(l, m) for given orbital.
    
    theta: float, polar angle in radians
    orbital: Orbital, the orbital object containing the parameters
    m: int, magnetic quantum number

    RETURN: float, the value of the spherical harmonics at angle theta
    """
    l = orbital.l
    if abs(m) > l:
        raise ValueError("m no puede ser mayor que |l|")
    
    # P_l^m(cos(theta)) -> theta es el ángulo polar (desde z+ hacia abajo)
    P_lm = lpmv(m, l, math.cos(theta))

    # Normalización
    A_ml = math.sqrt((2*l + 1) * math.factorial(l - m) / (4 * math.pi*math.factorial(l + m)))

    return A_ml * P_lm

    
def rho_monoparticular(r: float, atomic_data: AtomicData, space: str,  spherical_averaged = False, normalized = True) -> float:
    """
    Calculate the density of the atomic system. By default, it is normalized to 1 and not spherical averaged.
    
    r: float, the distance from the nucleus
    atomic_data: AtomicData, the atomic data object containing the orbitals
    space: str, either "position" or "momentum"
    spherical_averaged: bool, true if the density is spherical averaged, false otherwise
    normalized: bool, true if the density is normalized to 1, false if is normalized to number of electrons in the atom

    RETURN: float, the value of the density at distance r
    """
    if space not in ["position", "momentum", "my_version"]:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")
      
    density = 0.0
    N = 0.0
    for orbital in atomic_data.orbitals:
        N += orbital.mz 
        density += Rnl_monoparticular(r, orbital, space, False)**2

    if normalized and spherical_averaged:
        return density / (4* np.pi * N)
    elif normalized and not spherical_averaged:
        return density / N
    elif not normalized and spherical_averaged:
        return density / (4*np.pi)
    else:
        return density         

"""def Gamma_biparticular(r1: float, r2: float, atomic_data: AtomicData, space: str, spherical_averaged=False, normalized=True) -> float:
    if space not in ["position", "momentum"]:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")

    orbitals = atomic_data.orbitals
    N = sum(orbital.mz for orbital in orbitals)  # Total de electrones reales
    

    if N <= 1:
        return 0.0  # No tiene sentido calcular densidad biparticular

    # Generamos una lista con cada electrón individual (expandido según ocupación)
    expanded_orbitals = []
    for orbital in orbitals:
        expanded_orbitals.extend([orbital] * orbital.mz)

    density = 0.0

    # Doble suma sobre electrones individuales (i != j)
    for i in range(len(expanded_orbitals)):
        for j in range(len(expanded_orbitals)):
            if i != j:
                o_i = expanded_orbitals[i]
                o_j = expanded_orbitals[j]

                
                
                # phi_i(r1), phi_j(r2)
                phi_i_r1 = Rnl_monoparticular(r1, o_i, space, True)
                phi_j_r2 = Rnl_monoparticular(r2, o_j, space, True)

                direct = (phi_i_r1**2) *(phi_j_r2**2)
                #Si los electrones tienen el mismo spin, no se considera el intercambio
                if i%2 != j%2:
                    #print(i, j)
                    density += direct
                    continue
                # Intercambio: phi_j(r1), phi_i(r2)
                phi_j_r1 = Rnl_monoparticular(r1, o_j, space, True)
                phi_i_r2 = Rnl_monoparticular(r2, o_i, space, True)

                exchange = (phi_i_r1 * phi_j_r1) * (phi_j_r2 * phi_i_r2)

                density += direct - exchange
                #print(f"i: {i}, j: {j}, phi^2_{o_i.n}{o_i.l}(r1) phi^2_{o_j.n}{o_j.l}(r2) - phi_{o_i.n}{o_i.l}(r1) phi_{o_j.n}{o_j.l}(r1) phi_{o_j.n}{o_j.l}(r2) phi_{o_i.n}{o_i.l}(r2)")

    density *= 1 / (N - 1)

    if normalized and spherical_averaged:
        return density / ((16*np.pi**2)*N)
    elif normalized and not spherical_averaged:
        return density / N
    elif not normalized and spherical_averaged:
        return density / (16*np.pi**2)
    else:
        return density"""   


# Adaptamos Gamma_biparticular_updated con la lógica exacta del código Fortran

def Gamma_biparticular(
    r1: float,
    r2: float,
    atomic_data: AtomicData,
    space: str,
    spherical_averaged: bool = False,
    normalized: bool = True
) -> float:
    """
    Calcula Γ(r1, r2) usando la lógica del código Fortran: sumando sobre orbitales únicos,
    aplicando intercambio solo si l_i == l_j y i ≠ j, con normalización final 1/(16π⁴ N(N-1))
    """

    if space not in ["position", "momentum", "my_version"]:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")

    orbitals = atomic_data.orbitals
    N = sum(o.mz for o in orbitals)
    if N <= 1:
        return 0.0

    gamma = 0.0
    cont = len(orbitals)

    for i in range(cont):
        o_i = orbitals[i]
        noci, li = o_i.mz, o_i.l
        Ri_r1 = Rnl_monoparticular(r1, o_i, space, normalized=False) / math.sqrt(noci)
        Ri_r2 = Rnl_monoparticular(r2, o_i, space, normalized=False) / math.sqrt(noci)

        for j in range(i, cont):
            o_j = orbitals[j]
            nocj, lj = o_j.mz, o_j.l
            Rj_r1 = Rnl_monoparticular(r1, o_j, space, normalized=False) / math.sqrt(nocj)
            Rj_r2 = Rnl_monoparticular(r2, o_j, space, normalized=False) / math.sqrt(nocj)

            if i == j:
                mult = noci * (noci - 1) / 2
            else:
                mult = noci * nocj

            gamma += mult * ((Ri_r1 * Rj_r2)**2 + (Rj_r1 * Ri_r2)**2)

            if li == lj and i != j:
                h = min(noci, nocj)
                gamma -= 2 * h * Ri_r1 * Rj_r1 * Ri_r2 * Rj_r2

    gamma /=  (N - 1)

    if normalized and spherical_averaged:
        return gamma/ ((16*np.pi**2)*N)
    elif normalized and not spherical_averaged:
        return gamma / N
    elif not normalized and spherical_averaged:
        return gamma/ (16*np.pi**2)
    else:
        return gamma  


def plot_density(orbital: Orbital, include_angular: Optional[bool] = None, m: Optional[int] = None, space: Optional[str] = None) -> None:
    """
    Plotea la densidad del orbital. Si include_angular=False, se hace un gráfico 2D de la parte radial.
    Si include_angular=True, se hace un gráfico 3D de la densidad combinada radial y angular.
    """
    if include_angular is None:
        include_angular = False  # Default to False if not provided
    else:
        if include_angular: 
            if m is None:
                raise ValueError("m debe ser proporcionado si include_angular es True.")
            if abs(m) > orbital.l:
                raise ValueError("|m| no puede ser mayor que l")
    
    if space is None:
        space = "position"
        warnings.warn("Space not provided. Defaulting to 'position'.", UserWarning)
    if space not in ["position", "momentum"]:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")
    
    if not include_angular:
        # Plot 2D de la densidad radial
        # Si el orbital tiene n grande o es muy difuso
        r_max = 5 + 2 * orbital.l
        r_vals = np.linspace(0.01, r_max, 100)
        density = [Rnl_monoparticular(r, orbital, space)**2 for r in r_vals]

        plt.figure()
        plt.plot(r_vals, density)
        plt.xlabel('r')
        plt.ylabel('|Rnl(r)|²')
        plt.title('Densidad Radial del Orbital')
        plt.grid(True)
        plt.show()

    else:
        l = orbital.l
        print(f"m = {0}, l = {l}")
        # Plot 3D de la densidad radial-angular
        # Si el orbital tiene n grande o es muy difuso
        r_max = 5 + 2 * orbital.l
        r_vals = np.linspace(0.01, r_max, 100)
        theta_vals = np.linspace(0, np.pi, 50)
        phi_vals = np.linspace(0, 2 * np.pi,50)

        r, theta, phi = np.meshgrid(r_vals, theta_vals, phi_vals, indexing='ij')

        
        # Cálculo de densidad total
        radial_part = np.vectorize(lambda r_: (Rnl_monoparticular(r_, orbital, space))**2)(r)
        angular_part = np.vectorize(lambda th: (Ylm(th, orbital, 1))**2)(theta)
        density = radial_part * angular_part

        # Coordenadas cartesianas
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # Filtrar por densidad umbral para visualizar mejor
        threshold = np.max(density) * 0.01
        mask = density > threshold

        # Graficar con Plotly
        fig = go.Figure(data=go.Scatter3d(
            x=x[mask].flatten(),
            y=y[mask].flatten(),
            z=z[mask].flatten(),
            mode='markers',
            marker=dict(
                size=2,
                color=density[mask].flatten(),
                colorscale='Viridis',
                opacity=0.4,
                colorbar=dict(title='Densidad')
            )
        ))

        fig.update_layout(
            title='Densidad Orbital 3D Interactiva',
            scene=dict(
                xaxis_title='x',
                yaxis_title='y',
                zaxis_title='z'
            ),
            margin=dict(l=0, r=0, b=0, t=30)
        )

        fig.show()



