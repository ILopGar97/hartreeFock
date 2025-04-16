import numpy as np
import math
from data_structures import *
from scipy.special import lpmv 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go


def Rnl_monoparticular(r: float, Orbital: Orbital, space: str) -> float:
    """
    Calculate the radial wave function for a given orbital.
    The radial wave function is defined as:
    """
    A_list = Orbital.A
    C_list = Orbital.C
    N_list = Orbital.N
    result = 0.0
    mz = Orbital.mz

    if space == "position":
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):
            K_i = np.power(2*A_i, N_i + 0.5)/np.sqrt(math.factorial(2*N_i))
            result += C_i *K_i* (r ** (N_i - 1)) * np.exp(-A_i * r)/np.sqrt(mz) #Normalizado a 1
        return result
    elif space == "momentum":
        for N_i, A_i, C_i in zip(N_list, A_list, C_list):
            N = C_i*np.power(2*A_i, N_i + 0.5)/(np.sqrt(math.factorial(2*N_i))*np.sqrt(mz)) #Normalizado a 1
            Const1 = math.factorial(N_i)*math.sqrt(2*math.pi)/np.power(A_i, N_i + 2)
            Const2 = 1.0/(1+(r/A_i)**2)**(N_i + 1)
            sum = 0
            for k in range(0, math.floor(N_i/2)):
                sum += math.comb(N_i+1, 2*k + 1)*math.pow(-1, k)*math.power(r/A_i, 2*k)
            result += N*Const1*Const2*sum 
        return result
    else:
        raise ValueError("Invalid space. Use 'position' or 'momentum'.")
    
def Ylm(orbital: Orbital, theta: float, m: int) -> float: #se omite la fase compleja
    """
    Calculate the spherical harmonics Y(l, m) for given orbital.
    The spherical harmonics are defined as:
    """
    l = orbital.l
    if abs(m) > l:
        raise ValueError("m no puede ser mayor que |l|")

    m_abs = abs(m)

    # P_l^m(cos(theta)) -> theta es el ángulo polar (desde z+ hacia abajo)
    P_lm = lpmv(m_abs, l, math.cos(theta))

    # Normalización
    A_ml = math.sqrt((2*l + 1) * math.factorial(l - m_abs) / (4 * math.pi*math.factorial(l + m_abs)))

    return A_ml * P_lm
    


def plot_density(orbital: Orbital, include_angular=False, space='position'):
    """
    Plotea la densidad del orbital. Si include_angular=False, se hace un gráfico 2D de la parte radial.
    Si include_angular=True, se hace un gráfico 3D de la densidad combinada radial y angular.
    """
    print(f"Plotting density for orbital dergerg")
    if not include_angular:
        # Plot 2D de la densidad radial
        r_vals = np.linspace(0.01, 10, 500)
        density = [(r**2)*Rnl_monoparticular(r, orbital, space)**2 for r in r_vals]

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
        r_vals = np.linspace(0.01, 5, 50)
        theta_vals = np.linspace(0, np.pi, 50)
        phi_vals = np.linspace(0, 2 * np.pi,50)

        r, theta, phi = np.meshgrid(r_vals, theta_vals, phi_vals, indexing='ij')

        
        # Cálculo de densidad total
        radial_part = np.vectorize(lambda r_: (r**2)*(Rnl_monoparticular(r_, orbital, space))**2)(r)
        angular_part = np.vectorize(lambda th: math.sin(th)*(Ylm(orbital, th, 0))**2)(theta)
        density = radial_part * angular_part

        # Coordenadas cartesianas
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # Filtrar por densidad umbral para visualizar mejor
        threshold = np.max(density) * 0
        mask = density >= threshold

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



