# import numpy as np
# import time
# import matplotlib.pyplot as plt
# from plot_outplot import plot_outplot

# def outfun(x, optimValues, state):
#     """
#     Output function for optimization process.

#     Parameters:
#     x (np.ndarray): Current design variable values.
#     optimValues (dict): A dictionary containing optimization values.
#     state (str): The current state of the optimization ('init', 'iter', or 'done').

#     Returns:
#     bool: Indicates whether to stop the optimization (False in this case).
#     """
#     global history, nelx, nely, nelz

#     stop = False

#     if state == 'init':
#         tic = time.time()

#     elif state == 'iter':
#         # Calculate elapsed time
#         t = time.time() - tic
        
#         # Plot design variables
#         plot_outplot(x, t, optimValues, 1)

#         # Display iteration details
#         volume = np.sum(x[:len(x) // 2]) / (nelx * nely * nelz)
#         print(f" Iter: {optimValues['iteration']:4d} "
#               f"Obj: {optimValues['fval']:6.4f} "
#               f"Vol: {volume:6.4f} "
#               f"Step: {optimValues['stepsize']:6.4f} "
#               f"Time Elapsed: {t:6.4f}")

#     elif state == 'done':
#         # Record the iteration
#         history['iter'].append(optimValues['iteration'])

#         t = time.time() - tic
#         plot_outplot(x, t, optimValues)

#     return stop

# def plot_outplot(x, t, optimValues, flag=0):
#     """
#     Placeholder function for plotting the design variables.
#     To be implemented based on specific requirements.

#     Parameters:
#     x (np.ndarray): Current design variable values.
#     t (float): Time elapsed since the last update.
#     optimValues (dict): A dictionary containing optimization values.
#     flag (int): Flag to control plot behavior.
#     """
#     # Example plot implementation (modify as needed)
#     plt.figure()
#     plt.plot(x)
#     plt.title(f"Design Variables at Iteration {optimValues['iteration']}")
#     plt.xlabel('Variable Index')
#     plt.ylabel('Value')
#     plt.grid()
#     plt.show()

# import numpy as np
# import time
# import matplotlib.pyplot as plt
# from plot_outplot import plot_outplot

# def outfun(x, optimValues, state):
#     """
#     Fonction de sortie pour le processus d'optimisation.

#     Paramètres:
#     x (np.ndarray) : Valeurs actuelles des variables de conception.
#     optimValues (dict) : Un dictionnaire contenant les valeurs d'optimisation.
#     state (str) : L'état actuel de l'optimisation ('init', 'iter', 'done').

#     Retourne:
#     bool : Indique si l'optimisation doit être arrêtée (False dans ce cas).
#     """
#     global history, nelx, nely, nelz

#     stop = False  # L'optimisation ne s'arrête pas par défaut

#     if state == 'init':
#         tic = time.time()

#     elif state == 'iter':
#         # Calculer le temps écoulé
#         t = time.time() - tic
        
#         # Afficher les graphiques des variables de conception
#         plot_outplot(x, t, optimValues, 1)

#         # Afficher les détails de l'itération
#         volume = np.sum(x[:len(x) // 2]) / (nelx * nely * nelz)  # Volume des matériaux
#         print(f" Iter: {optimValues['iteration']:4d} "
#               f"Obj: {optimValues['fval']:6.4f} "
#               f"Vol: {volume:6.4f} "
#               f"Step: {optimValues['stepsize']:6.4f} "
#               f"Time Elapsed: {t:6.4f} "
#               f"Feval: {optimValues['feval']:4d} "
#               f"Status: {optimValues['status']}")

#     elif state == 'done':
#         # Enregistrer l'itération finale
#         history['iter'].append(optimValues['iteration'])

#         t = time.time() - tic


import time
import numpy as np
import time
from scipy.optimize import minimize
from plot_outplot import plot_outplot  # Assurez-vous que cette fonction est définie

# Dictionnaire global pour stocker les informations de l'optimisation
history = {}

# Variable globale pour stocker le temps initial (tic)
tic = None  # Initialisation à None pour vérifier si tic est déjà défini

def outfun(x, optimValues, state, nelx, nely, nelz, phi,video):
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
              f"Vol: {volume:6.4f} "
              f"Step: {optimValues['stepsize']:6.4f} "
              f"Time Elapsed: {t:6.4f} "
              f"Feval: {optimValues['feval']:4d} "
              f"Status: {optimValues['status']}")

    elif state == 'done':
        if tic is None:
            tic = time.time()  # Assurez-vous que tic est défini avant de l'utiliser

        # Enregistrer l'itération finale
        history['iter'].append(optimValues['iteration'])
        
        t = time.time() - tic
        print(f"Optimisation terminée. Temps total : {t:.2f} secondes.")
        print(f"Valeur finale de la fonction objectif : {optimValues['fval']:6.4f}")
        
    return stop  # Retourne False pour ne pas arrêter l'optimisation

# Fonction wrapper pour appeler outfun avec les arguments manquants
def outfun_wrapper(x, nelx, nely, nelz, phi,video):
    # Vous devez ici récupérer les informations de `optimValues` et `state` d'une manière adéquate.
    # Vous pouvez soit utiliser un objet global, soit utiliser un dictionnaire.
    
    # Simulation de optimValues (vous devez les récupérer à partir de l'optimiseur dans un vrai cas)
    optimValues = {
        'iteration': 1,  # Remplacez par les vraies valeurs pendant l'optimisation
        'fval': 0.5,
        'stepsize': 0.1,
        'feval': 2,
        'status': 'success'
    }

    # L'état peut être 'init', 'iter' ou 'done'
    state = 'iter'  # Par exemple, ici l'état serait 'iter' à chaque itération

    # Appeler la fonction outfun avec les arguments appropriés
    outfun(x, optimValues, state, nelx, nely, nelz, phi,video)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#     def outfun(x, optimValues, state):
#     """
#     Fonction de rappel pour l'optimisation.

#     Paramètres:
#     x (np.ndarray): Variables de conception.
#     optimValues (dict): Dictionnaire contenant les résultats de l'optimisation.
#     state (str): L'état actuel de l'optimisation, comme 'init', 'iter', 'done'.
#     """

#     # Initialisation
#     if state == 'init':
#         optimValues['iteration'] = 0
#         optimValues['fval'] = None
#         optimValues['volume'] = None
#         optimValues['stepsize'] = None
#         optimValues['status'] = None
#         optimValues['feval'] = 0

#     elif state == 'iter':
#         # Mettre à jour les valeurs à chaque itération
#         optimValues['iteration'] = optimValues.get('iteration', 0) + 1
#         optimValues['fval'] = np.sum(x)  # Exemple : fonction objectif
#         optimValues['volume'] = np.sum(x[:len(x) // 2]) / (nelx * nely * nelz)  # Exemple : volume
#         optimValues['stepsize'] = np.linalg.norm(x)  # Exemple : taille de pas
#         optimValues['status'] = 'success'  # Exemple : statut
#         optimValues['feval'] += 1  # Incrémenter le nombre de fonctions évaluées

#         # Afficher les informations de l'itération
#         print(f"Iter: {optimValues['iteration']} "
#               f"Obj: {optimValues['fval']:.4f} "
#               f"Vol: {optimValues['volume']:.4f} "
#               f"Step: {optimValues['stepsize']:.4f} "
#               f"Feval: {optimValues['feval']} "
#               f"Status: {optimValues['status']}")

#     elif state == 'done':
#         # Enregistrer les résultats finaux
#         print(f"Optimisation terminée: {optimValues['iteration']} itérations")








# from scipy.optimize import minimize
# import numpy as np

# # Définir la fonction objectif
# def cost_function(x):
#     return np.sum(x**2)  # Exemple simple de fonction objectif

# # Initialiser les variables de conception
# x0 = np.ones(10)

# # Définir les bornes pour les variables
# bounds = [(0, 1)] * 10

# # Initialiser un dictionnaire pour stocker les résultats de l'optimisation
# optimValues = {
#     'iteration': 0,
#     'fval': None,
#     'volume': None,
#     'stepsize': None,
#     'status': None,
#     'feval': 0
# }

# # Appeler la fonction minimize avec la fonction de rappel
# res = minimize(cost_function, x0, bounds=bounds, method='SLSQP', callback=lambda x: outfun(x, optimValues, 'iter'))

# # Afficher les résultats finaux
# print("Optimisation terminée")
# print(f"Nombre d'itérations : {optimValues['iteration']}")
# print(f"Valeur de la fonction objectif finale : {optimValues['fval']}")
# $