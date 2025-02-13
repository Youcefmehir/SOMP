%matplotlib qt
import numpy as np
# import matplotlib
# matplotlib.use('TkAgg')  # Ou 'Qt5Agg' selon le backend que vous préférez

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import cv2  # Pour l'enregistrement vidéo
import time
from Top_Obj import top_obj
from fea_input import fea_input
from scipy.sparse import coo_matrix

# Initialisation des variables
tic = None
if tic is None:
    tic = time.time()  # Si tic n'est pas initialisé, l'initialiser ici.

t = time.time() - tic

# Paramètres du modèle
nelx = 6
nely = 4
nelz = 2
rho0 = 0.5
theta0 = 0
p = 3
rmin = 1.5

# Initialisation de la densité et de l'angle
rho0 = rho0 * np.ones((nely, nelx, nelz))
theta0 = theta0 * np.ones((nely, nelx, nelz))

video = cv2.VideoWriter('need_change.avi', cv2.VideoWriter_fourcc(*'XVID'), 5, (640, 480))
nele = nelx * nely * nelz

# Initialisation de l'angle phi
offset = 0
a = -offset / 180 * np.pi
b = offset / 180 * np.pi
phi = a + (b - a) * np.random.rand(nele)
phi = phi.reshape(nely, nelx, nelz)

# Concatenation des variables
x = np.concatenate((rho0.flatten(), theta0.flatten()))   

# Matrices pour l'optimisation
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

# Conversion des listes en arrays numpy pour l'efficacité
iH = np.array(iH, dtype=int)
jH = np.array(jH, dtype=int)
sH = np.array(sH, dtype=float)
H_1 = np.array(H_1, dtype=float)

# Construction de la matrice creuse H
H = coo_matrix((sH, (iH - 1, jH - 1)), shape=(nele, nele))  # Soustraction de 1 pour l'indexation à partir de 0
Hs = np.array(H.sum(axis=1)).flatten()    

# Calcul de l'objet et de ses dérivées
F, dF = top_obj(x, nelx, nely, nelz, p, phi, Hs, H)

# Résultats d'optimisation
optimresult = {
    'iteration': 1,  # Remplacez par les vraies valeurs pendant l'optimisation
    'fval': F,
    'vol': 1,
}

# Extraction de la densité (rho) et de l'angle (theta)
rho = x[:len(x) // 2]
theta = x[len(x) // 2:]

# Reshaping de la densité et des angles
rho2 = np.reshape(rho, (nely, nelx, nelz))
theta2 = np.reshape(theta, (nely, nelx, nelz))

# Calcul des résultats de la simulation FEA
D, ieqn, utot, ptot, F, dcdrho, dcdtheta, youcef, ielem = fea_input(nelx, nely, nelz, rho2, theta2, phi)
XX = youcef

# Filtrage des éléments où rho est supérieur à 0.5
ielem1 = ielem[np.where(rho >= 0.5), :]

# Extraction des coordonnées
xx = XX[:, 0]
yy = XX[:, 1]
zz = XX[:, 2]

theta1 = theta[np.where(rho >= 0.5)]
phi1 = phi.flatten()[np.where(rho >= 0.5)]

# Définition des segments de fibres
center_xx = (np.mean(xx[ielem1], 2)).T
center_yy = (np.mean(yy[ielem1], 2)).T
center_zz = (np.mean(zz[ielem1], 2)).T

tip1x = center_xx + (0.25 * (np.cos(theta1) * np.cos(phi1)).reshape(-1, 1))
tip1y = center_yy + (0.25 * np.sin(theta1).reshape(-1, 1))
tip1z = center_zz - (0.25 * (np.cos(theta1) * np.sin(phi1)).reshape(-1, 1))

tip2x = center_xx - (0.25 * (np.cos(theta1) * np.cos(phi1)).reshape(-1, 1))
tip2y = center_yy - (0.25 * np.sin(theta1).reshape(-1, 1))
tip2z = center_zz + (0.25 * (np.cos(theta1) * np.sin(phi1)).reshape(-1, 1))

# Transformation des coordonnées
if nely >= nelz:
    xxx = np.vstack((center_xx, tip1x, tip2x))
    yyy = np.vstack((center_yy, tip1y, tip2y))
    zzz = np.vstack((center_zz, tip1z, tip2z))
    XXX = np.column_stack((xxx, yyy, zzz))
    XXX[:, [1, 2]] = XXX[:, [2, 1]]
    XXX[:, 1] = -XXX[:, 1]

    xxx = (XXX[:, 0].reshape(3, -1)).T
    yyy = (XXX[:, 1].reshape(3, -1)).T
    zzz = (XXX[:, 2].reshape(3, -1)).T
    center_xx, tip1x, tip2x = xxx[:, 0], xxx[:, 1], xxx[:, 2]
    center_yy, tip1y, tip2y = yyy[:, 0], yyy[:, 1], yyy[:, 2]
    center_zz, tip1z, tip2z = zzz[:, 0], zzz[:, 1], zzz[:, 2]

# Création du graphique 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Définition des faces pour les éléments (cube ici)
faces1 = np.array([
    [0, 1, 2, 3],  # Face 1
    [4, 5, 6, 7],  # Face 2
    [0, 1, 5, 4],  # Face 3
    [1, 2, 6, 5],  # Face 4
    [2, 3, 7, 6],  # Face 5
    [3, 0, 4, 7],  # Face 6
])

# Coordonnées des sommets
coords = np.array([
    [0, 0, 0],  # Sommet 0
    [1, 0, 0],  # Sommet 1
    [1, 1, 0],  # Sommet 2
    [0, 1, 0],  # Sommet 3
    [0, 0, 1],  # Sommet 4
    [1, 0, 1],  # Sommet 5
    [1, 1, 1],  # Sommet 6
    [0, 1, 1]   # Sommet 7
])

# Extraction des faces
faces_coords = [coords[face] for face in faces1]

# Couleurs par densité
rho1 = np.repeat(rho[np.where(rho >= 0.5)], 6, axis=0)
face_colors = plt.cm.viridis(1 - rho1)  # Utilisation de viridis pour les couleurs

# Ajout des faces 3D au graphique
poly3d = Poly3DCollection(faces_coords, facecolors=face_colors, alpha=0.5)
ax.add_collection3d(poly3d)

# Définir les limites des axes
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_zlim([0, 1])

# Affichage du graphique
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

# Enregistrement vidéo (si nécessaire)
video.release()

