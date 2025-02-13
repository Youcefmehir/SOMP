import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fea_input import fea_input
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation
# %matplotlib qt

# from scipy.sparse import coo_matrix, csr_matrix






def plot_outplot(x, t, optimresult, nelx, nely, nelz, phi,video,a=None):

    
    a=None
    """
    Plots the fiber orientation and density distribution in 3D.

    Parameters:
    x (np.ndarray): Design variables, where the first half corresponds to density 
                    (rho) and the second half corresponds to fiber angle (theta).
    t (float): Time elapsed in seconds.
    optimresult (object): Contains optimization results such as iteration count and function value.
    a: Additional argument for video recording (optional).
    """
    # Extract density (rho) and fiber angle (theta)
    rho = x[:len(x) // 2]
    theta = x[len(x) // 2:] 
    phi0 = phi.flatten()
    rho2 = np.reshape(rho, (nely, nelx, nelz))
    theta2 = np.reshape(theta, (nely, nelx, nelz))

    D, ieqn,utot, ptot, F, dcdrho, dcdtheta,X_X,ielem= fea_input(nelx, nely, nelz, rho2, theta2, phi)
    XX = X_X

    plt.figure(2)
    ielem1 = ielem[np.where(rho >= 0.5), :]
    

    # Extract coordinates
    xx= XX[:, 0]
    yy= XX[:, 1]
    zz= XX[:, 2]
    
    theta1 = theta[np.where(rho >= 0.5)]
    # phi1 = phi0[np.where(rho >= 0.5)]
    phi1 = phi.flatten()[np.where(rho >= 0.5)]

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

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
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
    
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
    ax.axis('off')
    plt.set_cmap('gray')
    ax.view_init(elev=30, azim=-60)

    plt.draw()

    if a is not None:  # If video recording is requested
        frame = plt.gcf()  # Get the current figure
        video.write(frame)

    plt.show()





 
    
    
    
    
    
    
    








    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    