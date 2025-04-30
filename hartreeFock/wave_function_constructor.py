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

    if space == "position":
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):
            K_i = np.power(2*A_i, N_i + 0.5)/np.sqrt(math.factorial(2*N_i))
            result += C_i *K_i* (r ** (N_i - 1)) * np.exp(-A_i * r) 

        if normalized:
            return result / np.sqrt(mz)
        else:
            return result 
        
    elif space == "momentum":
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):

            Norm = C_i*np.power(2*A_i, N_i + 0.5)/(np.sqrt(math.factorial(2*N_i))) 

            Const1 = math.factorial(N_i)*math.sqrt(2/math.pi)/np.power(A_i, N_i + 2)

            Const2 = 1.0/(1+(r/A_i)**2)**(N_i + 1)

            sum = 0
            for k in range(0, math.floor(N_i/2) + 1):
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

    
def rho_monoparticular(r: float, atomic_data: AtomicData, space: str, spherical_averaged = False, normalized = True) -> float:
    """
    Calculate the density of the atomic system. By default, it is normalized to 1 and not spherical averaged.
    
    r: float, the distance from the nucleus
    atomic_data: AtomicData, the atomic data object containing the orbitals
    space: str, either "position" or "momentum"
    spherical_averaged: bool, true if the density is spherical averaged, false otherwise
    normalized: bool, true if the density is normalized to 1, false if is normalized to number of electrons in the atom

    RETURN: float, the value of the density at distance r
    """
    if space not in ["position", "momentum"]:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")
      
    density = 0.0
    N = 0.0
    for orbital in atomic_data.orbitals:
        N += orbital.mz 
        if space == "position":
            density += Rnl_monoparticular(r, orbital, space)**2
        elif space == "momentum":
            density += Rnl_monoparticular(r, orbital, space)**2

    if normalized and spherical_averaged:
        return density / (4*np.pi*N)
    elif normalized and not spherical_averaged:
        return density / N
    elif not normalized and spherical_averaged:
        return density / (4*np.pi)
    else:
        return density         

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



