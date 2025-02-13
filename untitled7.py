# %matplotlib qt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fea_input import fea_input
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation
from scipy.sparse import coo_matrix, csr_matrix
from scipy.optimize import minimize
import cv2  # For video writing
from Top_Obj import top_obj
from plot_layer import plot_layer
import time
from plot_outplot import plot_outplot  # Assurez-vous que cette fonction est définie


tic = None

if tic is None:
    tic = time.time()  # Si tic n'est pas initialisé, l'initialiser ici.

t = time.time() - tic

nelx=2
nely=3
nelz=2
rho0=0.5
theta0=0
p=3
rmin=1.5

rho0 = rho0 * np.ones((nely, nelx, nelz))
theta0 = theta0 * np.ones((nely, nelx, nelz))


video = cv2.VideoWriter('need_change.avi', cv2.VideoWriter_fourcc(*'XVID'), 5, (640, 480))
nele = nelx * nely * nelz

offset = 0
a = -offset / 180 * np.pi
b = offset / 180 * np.pi
phi = a + (b - a) * np.random.rand(nele)
phi = phi.reshape(nely, nelx, nelz)

x= np.concatenate((rho0.flatten(), theta0.flatten()))   
 

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


F, dF=top_obj(x, nelx, nely, nelz, p, phi,Hs,H)
# a=None



optimresult = {
        'iteration': 1,  # Remplacez par les vraies valeurs pendant l'optimisation
        'fval': F,
        'vol': 1,}

    

   
    # Extract density (rho) and fiber angle (theta)
rho = x[:len(x) // 2]
theta = x[len(x) // 2:]
    

    
    
# phi0 = phi.flatten()
    

rho2 = np.reshape(rho, (nely, nelx, nelz))
theta2 = np.reshape(theta, (nely, nelx, nelz))

    
    
D, ieqn,utot, ptot, F, dcdrho, dcdtheta,youcef,ielem= fea_input(nelx, nely, nelz, rho2, theta2, phi)
XX = youcef

    

# plt.figure(2)
ielem1 = ielem[np.where(rho >= 0.5), :]

    
# Extract coordinates
xx= XX[:, 0]
yy= XX[:, 1]
zz= XX[:, 2]



theta1 = theta[np.where(rho >= 0.5)]
# phi1 = phi0[np.where(rho >= 0.5)]
phi1 = phi.flatten()[np.where(rho >= 0.5)]

xxxxx=xx[ielem1]
yyyyy=yy[ielem1]
zzzzz=zz[ielem1]



    # Define fiber segments
center_xx = (np.mean(xx[ielem1],2)).T
center_yy = (np.mean(yy[ielem1],2)).T
center_zz = (np.mean(zz[ielem1],2)).T
    

tip1x = center_xx + (0.25 * (np.cos(theta1) * np.cos(phi1)).reshape(-1, 1))
tip1y = center_yy + (0.25 * np.sin(theta1).reshape(-1, 1))
tip1z = center_zz - (0.25 * (np.cos(theta1) * np.sin(phi1)).reshape(-1, 1))
    
tip2x = center_xx - (0.25 * (np.cos(theta1) * np.cos(phi1)).reshape(-1, 1))
tip2y = center_yy - (0.25 * np.sin(theta1).reshape(-1, 1))
tip2z = center_zz + (0.25 * (np.cos(theta1) * np.sin(phi1)).reshape(-1, 1))


    
    

    # Coordinate transformation
depth,row, col = ielem1.shape
if nely >= nelz:
    xxx = np.vstack((center_xx, tip1x, tip2x))
    yyy = np.vstack((center_yy, tip1y, tip2y))
    zzz = np.vstack((center_zz, tip1z, tip2z))
    XXX = np.column_stack((xxx, yyy, zzz))
    XXX[:, [1, 2]] = XXX[:, [2, 1]]
    XXX[:, 1] = -XXX[:, 1]


        
    xxx = (XXX[:, 0].reshape(3,row)).T
    yyy = (XXX[:, 1].reshape(3,row)).T
    zzz = (XXX[:, 2].reshape(3,row)).T
    center_xx, tip1x, tip2x = xxx[:, 0], xxx[:, 1], xxx[:, 2]
    center_yy, tip1y, tip2y = yyy[:, 0], yyy[:, 1], yyy[:, 2]
    center_zz, tip1z, tip2z = zzz[:, 0], zzz[:, 1], zzz[:, 2]



# Création du graphique 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# # Plot fiber orientation
fiber_x = np.array([tip2x, center_xx, tip1x])
fiber_y = np.array([tip2y, center_yy, tip1y])
fiber_z = np.array([tip2z, center_zz, tip1z])
ax = plt.axes(projection='3d')
# ax.plot3D(fiber_x.flatten(), fiber_y.flatten(), fiber_z.flatten(), color='red', linewidth=2)

for i in range(fiber_x.shape[1]):  # Parcourir chaque segment
    # Extraire les trois points du segment
    x_coords = fiber_x[:, i]
    y_coords = fiber_y[:, i]
    z_coords = fiber_z[:, i]
    
    # Créer un polygone qui relie les trois points
    verts = [list(zip(x_coords, y_coords, z_coords))]  # Liste de coordonnées des points
    ax.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))

ax.set_xlim([fiber_x.min() - 1, fiber_x.max() + 1])
ax.set_ylim([fiber_y.min() - 1, fiber_y.max() + 1])
ax.set_zlim([fiber_z.min() - 1, fiber_z.max() + 1])

ax.set_title(f'time: {t} sec; iter: {optimresult['iteration']}; '
                  f'compliance: {optimresult['fval']}; Mnd(%): '
                  f'{np.sum(4 * rho * (1 - rho)) / len(rho) * 100:.2f}')






if nely >= nelz:
    XX[:, [1, 2]] = XX[:, [2, 1]]  # Coordinate transformation
    XX[:, 1] = -XX[:, 1]
rho1 = np.repeat(rho[np.where(rho >= 0.5)], 6, axis=0).reshape(-1,1)
#
    # Create faces for elements
faces1 = np.vstack((
        ielem1[0, :, [0, 1, 2, 3]].T,  # front face
        ielem1[0, :, [4, 5, 6, 7]].T,  # back face
        ielem1[0, :, [3, 2, 7, 6]].T,  # right face
        ielem1[0, :, [0, 1, 5, 4]].T,  # left face
        ielem1[0, :, [0, 4, 7, 3]].T,  # bottom face
        ielem1[0, :, [1, 5, 6, 2]].T   # top face
    ))
    
    # coords = np.array([
    #         [0, 0, 0],  # Sommet 0
    #         [1, 0, 0],  # Sommet 1
    #         [1, 1, 0],  # Sommet 2
    #         [0, 1, 0],  # Sommet 3
    #         [0, 0, 1],  # Sommet 4
    #         [1, 0, 1],  # Sommet 5
    #         [1, 1, 1],  # Sommet 6
    #         [0, 1, 1]   # Sommet 7
    #         ])
    
    
    # coords = np.array([[x, y, z] for x in range(6) for y in range(4) for z in range(2)])
# coords=[]
# for k in range (nelz+1):
#         for i in range (nelx+1):
#             for j in range (nely+1):
#                 coords.append(np.array([i, j, k]))
                
                
#     # for face in faces1:
#     #     if np.any(face >= len(coords)):
#     #         print(f"Index out of bounds: {face}")
    
#     # faces_coords = [coords[face] for face in faces1]

#     # Convertir `coords` en un tableau numpy
# coords = np.array(coords)
coords=XX
    
    # Créer une nouvelle matrice pour lier chaque nœud à ses coordonnées
faces_coords = []
    
    # Lier chaque nœud de faces1 à ses coordonnées respectives
for face in faces1:
        node_coords = [coords[node] for node in face]  # Récupérer les coordonnées de chaque nœud dans la face
        faces_coords.append(node_coords)
    
    # Convertir le résultat en un tableau numpy
faces_coords = np.array(faces_coords)

ax.add_collection3d(Poly3DCollection(faces_coords, facecolors=plt.cm.gray(1 - rho1), alpha=0.5))

    # Set the box aspect ratio and hide the axes
ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
ax.axis('off')  # Hide axes
    
    # plt.view_init(elev=30, azim=30)
ax.view_init(elev=30, azim=30)




    # # Plot using Poly3DCollection
    # ax.add_collection3d(Poly3DCollection(faces_coords, facecolors=plt.cm.viridis(1 - rho1), alpha=0.5))

    # # Set axes limits and hide axes
    # ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
    # ax.axis('off')  # Hide axes
    # plt.colormaps('gray')
    # plt.view_init(elev=30, azim=30)


    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    
    #     # Transformation des coordonnées si nely >= nelz
    # if nely >= nelz:
    #     XX[:, [1, 2]] = XX[:, [2, 1]]
    #     XX[:, 1] = -XX[:, 1]
    
    #     # Plot de la densité
    # rho1 = np.tile(rho[rho > 0.1], (6, 1))  # Répéter pour chaque face
        
        
        
        
        

    # if nely >= nelz:
    #     XX[:, [1, 2]] = XX[:, [2, 1]]
    #     XX[:, 1] = -XX[:, 1]

    # # Plot density
    # rho1 = np.repeat(rho[np.where(rho >= 0.5)], 6, axis=0)
    # faces1 = np.vstack((ielem1[:, [0, 1, 2, 3]],
    #                     ielem1[:, [4, 5, 6, 7]],
    #                     ielem1[:, [3, 2, 7, 6]],
    #                     ielem1[:, [0, 1, 5, 4]],
    #                     ielem1[:, [0, 4, 7, 3]],
    #                     ielem1[:, [1, 5, 6, 2]]))
    
    
    
       
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')


    



    # Plot the 3D elements with material distribution
    # ax.add_collection3d(plt.Polygon(faces1, facecolors=1 - rho1, alpha=0.5))
    # ax.add_collection3d(Poly3DCollection(faces1, facecolors=plt.cm.viridis(1 - rho1), alpha=0.5))



    # ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
    # ax.axis('off')
    # plt.colormaps('gray')
    # plt.view_init(elev=30, azim=30)

    # plt.draw()

    # if a is not None:  # If video recording is requested
    #     frame = plt.gcf()  # Get the current figure
    #     video.write(frame)

    # plt.show()



