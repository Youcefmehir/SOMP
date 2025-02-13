import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import cv2  # For video writing
from Top_Obj import top_obj
from plot_layer import plot_layer
from outfun import outfun_wrapper
from fea_input import fea_input
import time
from plot_outplot import plot_outplot  # Assurez-vous que cette fonction est définie


optimValues = {
    'iteration': 0,
    'fval': None,
    'feval': 0
}

video = cv2.VideoWriter('need_change.avi', cv2.VideoWriter_fourcc(*'XVID'), 5, (640, 480))
# Dictionnaire global pour stocker les informations de l'optimisation
history = {}

# Variable globale pour stocker le temps initial (tic)
tic = None  # Initialisation à None pour vérifier si tic est déjà défini

def outfun(x,res, nelx, nely, nelz, phi):
    """
    Fonction de sortie pour le processus d'optimisation.
    
    Paramètres:
    x (np.ndarray): Valeurs actuelles des variables de conception.
    optimValues (dict): Dictionnaire contenant les informations d'optimisation.
    state (str): État actuel de l'optimisation ('init', 'iter', 'done').

    Retourne:
    bool: Indique si l'optimisation doit être arrêtée (False pour continuer).
    """
    global history, tic

    optimValues = {
        'iteration': res.nit,  # Remplacez par les vraies valeurs pendant l'optimisation
        'fval': res.fun,
        'vol': np.sum(x[:len(x)//2]) / (nelx * nely * nelz),
    }
    state=res.status
    
    # L'état peut être 'init', 'iter' ou 'done'
    state = 'iter'  # Par exemple, ici l'état serait 'iter' à chaque itération

    stop = False  # L'optimisation ne s'arrête pas par défaut

    if state == 'init':
        # Initialisation du chronomètre si ce n'est pas déjà fait
        tic = time.time()

    elif state == 'iter':
        if tic is None:
            tic = time.time()  # Si tic n'est pas initialisé, l'initialiser ici.
        
        # Calculer le temps écoulé
        t = time.time() - tic
        
        
        
        # Afficher les graphiques des variables de conception
        plot_outplot(x, t, optimValues, nelx, nely, nelz, phi,video, 1)

        # Calcul du volume des matériaux
        volume = np.sum(x[:len(x) // 2]) / (nelx * nely * nelz)  # Calcul du volume

        # Affichage des détails de l'itération
        print(f" Iter: {optimValues['iteration']:4d} "
              f"Obj: {optimValues['fval']:6.4f} "
              f"Vol: {optimValues['vol']:6.4f} "
              f"Time Elapsed: {t:6.4f} ")

    elif state == 'done':
        if tic is None:
            tic = time.time()  # Assurez-vous que tic est défini avant de l'utiliser

        # Enregistrer l'itération finale
        history['iter'].append(optimValues['iteration'])
        
        t = time.time() - tic
        print(f"Optimisation terminée. Temps total : {t:.2f} secondes.")
        print(f"Valeur finale de la fonction objectif : {optimValues['fval']:6.4f}")
        
    return stop  # Retourne False pour ne pas arrêter l'optimisation












def top_CFAO(nx, ny, nz, rho0, theta0, p, r):
    """
    Initialize the optimization process for topology optimization.

    Parameters:
    nx (int): Number of elements in x direction.
    ny (int): Number of elements in y direction.
    nz (int): Number of elements in z direction.
    rho0 (float): Initial density.
    theta0 (float): Initial angle.
    p (float): Penalization factor.
    r (float): Filter radius.
    """
    
    global penal, rmin, history, video, c0, passive_ele
    global H1, H_01
    
    
    


    history = {'iter': []}
    video = cv2.VideoWriter('need_change.avi', cv2.VideoWriter_fourcc(*'XVID'), 5, (640, 480))
    nelx, nely, nelz = nx, ny, nz
    penal = p
    rmin = r
    c0 = rho0

    
    # Precalculate Hs for sensitivity filter
    nele = nelx * nely * nelz
    iH = []
    jH = []
    sH = []
    H_1 = []
    
    k = 0
    for k1 in range(1, nelz + 1):  # Adjust for 1-based indexing (MATLAB)
        for i1 in range(1, nelx + 1):
            for j1 in range(1, nely + 1):
                e1 = (k1 - 1) * nelx * nely + (i1 - 1) * nely + j1
                for k2 in range(max(k1 - int(np.ceil(rmin)) + 1, 1), min(k1 + int(np.ceil(rmin)) - 1, nelz) + 1):
                    for i2 in range(max(i1 - int(np.ceil(rmin)) + 1, 1), min(i1 + int(np.ceil(rmin)) - 1, nelx) + 1):
                        for j2 in range(max(j1 - int(np.ceil(rmin)) + 1, 1), min(j1 + int(np.ceil(rmin)) - 1, nely) + 1):
                            e2 = (k2 - 1) * nelx * nely + (i2 - 1) * nely + j2
                            k += 1
                            iH.append(e1)
                            jH.append(e2)
                            sH.append(max(0, rmin - np.sqrt((i1 - i2) ** 2 + (j1 - j2) ** 2 + (k1 - k2) ** 2)))
                            H_1.append(1)
    
    # Convert lists to numpy arrays for efficient sparse matrix handling
    iH = np.array(iH, dtype=int)
    jH = np.array(jH, dtype=int)
    sH = np.array(sH, dtype=float)
    H_1 = np.array(H_1, dtype=float)
    
    # Build sparse matrix H
    H = coo_matrix((sH, (iH - 1, jH - 1)), shape=(nele, nele))  # Subtract 1 for 0-based indexing
    Hs = np.array(H.sum(axis=1)).flatten()
    
    # For mean sensitivity
    H_01 = coo_matrix((H_1, (iH - 1, jH - 1)), shape=(nele, nele))  # Subtract 1 for 0-based indexing
    H1 = np.array(H_01.sum(axis=1)).flatten()
    

    # Options for optimization
    # options = {
    #     'maxiter': 1000,
    #     'disp': True,
    #     'gtol': 1e-3,
    #     # 'callback': outfun
    # }
    
    options = {
    'disp': True,  # Affichage des informations pendant l'optimisation
    'maxiter': 100,  # Nombre maximum d'itérations
    'gtol': 1e-2,  # Tolérance sur la variation de la fonction objectif
    'xtol': 1e-2,  # Tolérance sur la variation des variables de conception
    }
    
    # Initialize densities and angles
    rho0 = rho0 * np.ones((nely, nelx, nelz))
    theta0 = theta0 * np.ones((nely, nelx, nelz))

    # Define phi (rotation about y axis)
    offset = 0
    a = -offset / 180 * np.pi
    b = offset / 180 * np.pi
    phi = a + (b - a) * np.random.rand(nele)
    phi = phi.reshape(nely, nelx, nelz)

    # Initial guess for optimization
    x0 = np.concatenate((rho0.flatten(), theta0.flatten()))
    lb = np.concatenate((1e-6 * np.ones_like(rho0.flatten()), -2 * np.pi * np.ones_like(theta0.flatten())))
    ub = np.concatenate((np.ones_like(rho0.flatten()), 2 * np.pi * np.ones_like(theta0.flatten())))

    # Equality constraint
    Aeq = np.concatenate((np.ones_like(rho0.flatten()), np.zeros_like(theta0.flatten())))
    beq = nelx * nely * nelz * c0

    # Optimize
    # res = minimize(lambda x: top_obj(x, nelx, nely, nelz, p, phi,Hs,H), x0, method='SLSQP', bounds=list(zip(lb, ub)), constraints={'type': 'eq', 'fun': lambda x: np.dot(Aeq, x) - beq}, options=options, callback=outfun)


    xx=top_obj(x0, nelx, nely, nelz, p, phi,Hs,H)[0]
    print(xx)
    
    res = minimize(
        lambda x: top_obj(x, nelx, nely, nelz, p, phi,Hs,H)[0],
        x0, 
        method='trust-constr',#'SLSQP',
        bounds=list(zip(lb, ub)),
        # jac=lambda x: top_obj(x, nelx, nely, nelz, p, phi,Hs,H)[1],
        constraints={'type': 'eq', 'fun': lambda x: np.dot(Aeq, x) - beq},
        options=options,
        callback=lambda x, res: outfun(x,res, nelx, nely, nelz, phi))#(x,res))#, nelx, nely, nelz, p, phi,Hs,H,video))
    
    # Après optimisation, vous pouvez accéder aux résultats dans `res`
    print(f"Optimisation terminée avec {res.fun} itérations.")
    

    # Close video
    video.release()
    
    xloc, yloc, zloc = np.meshgrid(np.arange(nelx + 1), np.arange(nely, -1, -1), np.arange(nelz + 1))
    
    # Réorganiser les axes pour correspondre aux résultats MATLAB
    # xloc = xloc.T  # Transpose pour correspondre à MATLAB
    # yloc = yloc.T
    # zloc = zloc.T    
    
    #     # Si vous souhaitez également obtenir la liste des coordonnées sous forme de vecteurs
    # X_X = np.vstack((xloc.ravel(), yloc.ravel(), zloc.ravel())).T
    
    
    # Plot final layer
    #plot_layer(res.x, X_X)


if __name__ == "__main__":
    top_CFAO(nx=20, ny=6, nz=4, rho0=0.5, theta0=0, p=3, r=1.5)

