import numpy as np
from scipy.sparse import coo_matrix

def top_nonlcon(x):
    """
    Nonlinear constraint function for optimization.

    Parameters:
    x (np.ndarray): Design variable vector, where the first half represents rho values.

    Returns:
    c (np.ndarray): Inequality constraint values.
    ceq (np.ndarray): Equality constraint values.
    Gc (scipy.sparse.coo_matrix): Gradient of inequality constraints (if required).
    Geq (scipy.sparse.coo_matrix): Gradient of equality constraints (if required).
    """
    global passive_ele
    
    # Extract rho values
    rhop = x[:len(x) // 2]

    # Inequality constraints (none in this case)
    c = np.array([])

    # Equality constraints
    ceq = rhop[passive_ele] - 0.6

    # Gradient calculations if requested
    if 'Gc' in locals() or 'Geq' in locals():
        Gc = np.array([])  # No inequality constraints, so gradient is empty.
        
        # Gradient of equality constraints
        Geq = coo_matrix((np.ones(len(passive_ele)), 
                          (passive_ele, np.arange(len(passive_ele)))), 
                          shape=(len(x), len(ceq)))

        return c, ceq, Gc, Geq

    return c, ceq
