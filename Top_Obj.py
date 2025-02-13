import numpy as np

from fea_input import fea_input


def top_obj(x, nelx, nely, nelz, p, phi,Hs,H):
    """
    Objective function for topology optimization.

    Parameters:
    x (ndarray): Design variables, including density and angles.

    Returns:
    tuple: Objective function value and its sensitivities.
    """
    global ieqn, penal, D
    
    penal = p
    # Define design variables
    rho = x[:len(x) // 2]
    theta = x[len(x) // 2:]


    # Reshaping rho and theta to 3D arrays (nely, nelx, nelz)
    rho = np.reshape(rho, (nely, nelx, nelz))
    theta = np.reshape(theta, (nely, nelx, nelz))

    # Read input file
    D, ieqn,utot, ptot, F, dcdrho, dcdtheta,X_X,ielem= fea_input(nelx, nely, nelz, rho, theta, phi)
    D['p'] = penal                            

    # Filter sensitivities
    dcdrho = check(rho, dcdrho, Hs, H, 1,nelx, nely, nelz)

    # Design sensitivity
    dF = np.concatenate((dcdrho.flatten(), dcdtheta.flatten()))

    return F, dF

# MESH-INDEPENDENCY FILTER for density
def check(rho, dcdrho, Hs, H, flag,nelx, nely, nelz):
    """
    Applies a sensitivity filter to ensure mesh-independence.

    Parameters:
    rho (ndarray): Current design variable density values.
    dcdrho (ndarray): Sensitivity values for density.
    Hs (ndarray): Sum of filter coefficients.
    H (csr_matrix): Sparse matrix for filtering.
    flag (int): Flag to determine the type of filter.

    Returns:
    ndarray: Filtered sensitivity values.
    """
    global H1, H_01

    if flag == 1:
        # Original sensitivity filter
        dcn = H.dot(rho.flatten() * dcdrho.flatten()) / Hs / rho.flatten()
    elif flag == 2:
        # Modified sensitivity filter
        dcn = H.dot(rho.flatten() * dcdrho.flatten()) / (H.dot(rho.flatten()))
    elif flag == 3:
        # Mean sensitivity filter
        dcn = H_01.dot(dcdrho.flatten()) / H1

    return dcn.reshape((nely, nelx, nelz))
