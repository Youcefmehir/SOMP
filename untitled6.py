import numpy as np
from scipy.sparse import coo_matrix
from elemental_3D import elemental_3D
from numpy.linalg import pinv

nelx, nely, nelz,  rho, theta, phi=6,4,2,0.5,0.1,0
    

D = {}
D['nnode'] = (nelx + 1) * (nely + 1) * (nelz + 1)  # number of nodes in model
D['nel'] = nelx * nely * nelz  # number of elements in model
D['etype'] = 3  # element type (only one per model)
D['nelx'] = nelx
D['nely'] = nely
D['nelz'] = nelz


nnode = (nelx + 1) * (nely + 1) * (nelz + 1)  # number of nodes in model
nel= nelx * nely * nelz  # number of elements in model
etype = 3  # element type (only one per model)


    # Set properties related to element properties
D['elname'] = 'elemental_3D'
D['ndof'] = 3
D['nenode'] = 8
    
#    D['elname'] = 'elemental_3D'
ndof = 3
nenode = 8

#     # Nodal coordinates mesh (standard rectangular)
# xloc, yloc, zloc = np.meshgrid(np.arange(nelx + 1), np.arange(nely, -1, -1), np.arange(nelz + 1))
    
# X_X = np.column_stack((xloc.ravel(), yloc.ravel(), zloc.ravel()))
#     #X_X= np.vstack([xloc.ravel(), yloc.ravel(), zloc.ravel()]).T
    
# Création du meshgrid en Python (corrigé pour correspondre à MATLAB)
xloc1, yloc1, zloc1 = np.meshgrid(np.arange(nelx + 1), np.arange(nely, -1, -1), np.arange(nelz + 1))

# Réorganiser les axes pour correspondre aux résultats MATLAB
xloc = xloc1.T  # Transpose pour correspondre à MATLAB
yloc = yloc1.T
zloc = zloc1.T    

    # Si vous souhaitez également obtenir la liste des coordonnées sous forme de vecteurs
X_X = np.vstack((xloc.ravel(), yloc.ravel(), zloc.ravel())).T

    

    # Créer un vecteur allant de 1 à nnode
node_vect = np.arange(1, nnode + 1)
    
    # Reshaper le vecteur en un tableau 3D
node_mat = node_vect.reshape((nely + 1), (nelx + 1), (nelz + 1), order='F')
    
    # Affichage du résultat
#print(node_mat)

    
    #node_mat = np.reshape(np.arange(1, nnode + 1), (nely + 1, nelx + 1, nelz + 1))
    #node_mat = np.arange(1, nnode + 1).reshape((nely + 1, nelx + 1, nelz + 1))

    # Element data (standard rectangular)
ielem = np.zeros((nel, nenode), dtype=int)
icont = 0

for k in range(nelz):
    for i in range(nelx):
        for j in range(nely):
            icont += 1
            ielem[icont - 1, 0] = (nely + 1) * (i) + (j + 1) + (k) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 1] = (nely + 1) * (i + 1) + (j + 1) + (k) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 2] = (nely + 1) * (i + 1) + (j) + (k) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 3] = (nely + 1) * (i) + (j) + (k) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 4] = (nely + 1) * (i) + (j + 1) + (k + 1) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 5] = (nely + 1) * (i + 1) + (j + 1) + (k + 1) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 6] = (nely + 1) * (i + 1) + (j) + (k + 1) * (nelx + 1) * (nely + 1)
            ielem[icont - 1, 7] = (nely + 1) * (i) + (j) + (k + 1) * (nelx + 1) * (nely + 1)

iprops = np.ones(len(ielem), dtype=int)

    # Cantilever beam (force in y-direction) edge forces
D['loadcase'] = 1  # number of loadcases
D['nforce'] =  1  # number of applied forces per loadcase
D['ndisp'] = 3 * (nely + 1) * (nelz + 1)

    # Read applied loads
fnode = node_mat[-1, -1, np.floor(nelz // 2).astype(int)]
iforce = np.column_stack((fnode.flatten(), np.full(fnode.size, 2)))
force = -np.ones(D['nforce'])

    # Read applied displacements
unode = np.transpose(node_mat[:, 0,:])
idisp = np.column_stack([np.tile(unode.flatten(), 3), np.repeat([1, 2, 3], len(unode.flatten()))])
disp = np.zeros(D['ndisp'])  # Initialisation des déplacements


    # Set total degrees of freedom
ndoftot = nnode * ndof

    # Initialize ieqn
ieqn = np.zeros(ndoftot)

    # Identify fixed degrees of freedom
iessential = 0
 
for i in range(D['ndisp']):
    iessential += 1
        # Ajustement de l'indexation pour correspondre à MATLAB (indexation de 1)
    num = (idisp[i, 0] ) * ndof + (idisp[i, 1] -1) - ndof  # Ajustement -1 pour la base 1
    ieqn[num] = -iessential


    # Identify free degrees of freedom
ifree = 0
for i in range(ndoftot):
    if ieqn[i] == 0:
        ieqn[i] = ifree+1
        ifree += 1

    # Determine number of element free and fixed nodes
ndofff = 0
ndofee = 0
ndoffe = 0
for i in range(nel):
    jj = 0

    ieleqn = np.zeros(nenode* ndof)  # Initialize local equation numbers
    for j1 in range(nenode):
        for j2 in range(ndof):
            jj += 1
            numj = (ielem[i, j1] +1)* ndof + j2 - ndof
               
            ieleqn[jj - 1] = ieqn[numj]  # Store local equation numbers
                
    ndofff += np.sum(ieleqn > 0) ** 2
    ndofee += np.sum(ieleqn < 0) ** 2
    ndoffe += np.sum(ieleqn > 0) * np.sum(ieleqn < 0)
        

    # Initialize global matrices
Kindex = np.array([0, 0, 0])
Kffi = np.zeros(ndofff)
Kffj = np.zeros(ndofff)
Kff = np.zeros(ndofff)
Keei = np.zeros(ndofee)
Keej = np.zeros(ndofee)
Kee = np.zeros(ndofee)
Kfei = np.zeros(ndoffe)
Kfej = np.zeros(ndoffe)
Kfe = np.zeros(ndoffe)
    
Uf = np.zeros((nnode * ndof - D['ndisp'], D['loadcase']))
Ue = np.tile(disp.reshape(-1, 1), (1, D['loadcase']))
Ff = np.zeros((nnode * ndof - D['ndisp'], D['loadcase']))
Fe = np.zeros((D['ndisp'], D['loadcase']))

   





    
    
rho = rho * np.ones((nely, nelx, nelz))
theta = theta * np.ones((nely, nelx, nelz))


offset = 0
a = -offset / 180 * np.pi
b = offset / 180 * np.pi
phi = a + (b - a) * np.random.rand(D['nel'])
phi = phi.reshape(nely, nelx, nelz)
    
    





rho = rho.flatten()
theta = theta.flatten()
phi = phi.flatten()

    # Compute and assemble element stiffness matrix and element load vector
for i in range(D['nel']):
    elemnodes = ielem[i, :]
        # Get element dof number array
        
    ieldof = []  # Initialize the list to hold global dofs

    for j1 in range(nenode):
        for j2 in range(ndof):  # j2 goes from 1 to ndof inclusive

            global_dof = (elemnodes[j1]+1) * ndof + j2 - ndof
            ieldof.append(global_dof)
    #ieldof = get_ieldof(elemnodes, D)

        # Calculate element stiffness matrix and load vector
    eldat_F, eldat_K ,eldat_KK= eval(f"{D['elname']}(rho[i], theta[i], phi[i], elemnodes, D, -1 ,0, X_X)")
        #eldat_F, eldat_K = elemental_3D(rho[i], theta[i], phi[i], elemnodes, D, -1)
        
        # Assemble element stiffness matrix and load vector
    #assemkp_sparse(eldat_F, eldat_K, ieqn[ieldof], D, Kff, Kffi, Kffj, Kee, Keei, Keej, Kfe, Kfei, Kfej, Ff, Kindex)
    ieqn_1 = ieqn[ieldof]
    


    neldoftot = ndof * nenode  # nombre total de degrés de liberté par élément
    free = np.where(ieqn_1 > 0)[0]  # degrés de liberté libres de l'élément
    essential = np.where(ieqn_1 < 0)[0]  # degrés de liberté essentiels de l'élément
    
    # Assemblage des tableaux globaux
    # Assemblage du vecteur de charge
    # Ff[ieleqn[free]] += eldat_F[free]
  
    # Assemblage de la matrice de rigidité
    # Définition des variables d'index
    Lff = len(free)
    Lee = len(essential)
    Lfe = (neldoftot - Lff - Lee) // 2
    
    eldati = np.tile(ieqn_1[:, np.newaxis], (1, neldoftot))
    eldatj = np.tile(ieqn_1, (neldoftot, 1))
     
    # Extraction des composantes de la sous-matrice libre-libre
    Eff = eldat_K[np.ix_(free, free)].flatten()
    Effi = eldati[np.ix_(free, free)].flatten()
    Effj = eldatj[np.ix_(free, free)].flatten()
    
    # Extraction des composantes de la sous-matrice essentielle-essentielle
    Eee = eldat_K[np.ix_(essential, essential)].flatten()
    Eeei = -eldati[np.ix_(essential, essential)].flatten()
    Eeej = -eldatj[np.ix_(essential, essential)].flatten()
    
    # Extraction des composantes de la sous-matrice libre-essentielle
    Efe = eldat_K[np.ix_(free, essential)].flatten()
    Efei = np.abs(eldati[np.ix_(free, essential)].flatten())
    Efej = np.abs(eldatj[np.ix_(free, essential)].flatten())
    
    # Assemblage des composantes et des indices dans la forme sparse
    Kff[Kindex[0]:Kindex[0] + Lff**2] = Eff
    Kffi[Kindex[0]:Kindex[0] + Lff**2] = Effi
    Kffj[Kindex[0]:Kindex[0] + Lff**2] = Effj
    
    Kee[Kindex[1]:Kindex[1] + Lee**2] = Eee
    Keei[Kindex[1]:Kindex[1] + Lee**2] = Eeei
    Keej[Kindex[1]:Kindex[1] + Lee**2] = Eeej
    
    Kfe[Kindex[2]:Kindex[2] + Lff * Lee] = Efe
    Kfei[Kindex[2]:Kindex[2] + Lff * Lee] = Efei
    Kfej[Kindex[2]:Kindex[2] + Lff * Lee] = Efej
    
    # Mise à jour de Kindex
    Kindex += np.array([Lff**2, Lee**2, Lff * Lee])
    # Apply nodal forces
    
  
    
 
k = 0
for j in range(D['loadcase']):
    for i in range(D['nforce']):
        k += 1
        num = iforce[k - 1, 0] * D['ndof'] + iforce[k - 1, 1] - D['ndof']  # MATLAB to Python index adjustment
        Ff[int(ieqn[num-1]-1), int(j)] += force[k - 1]







# Vérification de la taille des indices par rapport à la dimension de la matrice
max_index = max(Kffi.max(), Kffj.max())
matrix_size = D['nnode'] * D['ndof'] - D['ndisp']


Kffi = Kffi - 1
Kffj = Kffj - 1
Keei, Keej=Keei-1, Keej-1
Kfei, Kfej=Kfei-1, Kfej-1

    # Define sparse matrices
Kff = coo_matrix((Kff, (Kffi, Kffj)), shape=(D['nnode'] * D['ndof'] - D['ndisp'], D['nnode'] * D['ndof'] - D['ndisp']))
Kee = coo_matrix((Kee, (Keei, Keej)), shape=(D['ndisp'], D['ndisp']))
Kfe = coo_matrix((Kfe, (Kfei, Kfej)), shape=(D['nnode'] * D['ndof'] - D['ndisp'], D['ndisp']))




Kff_inv = pinv(Kff.toarray())  # Conversion de la matrice sparse en dense pour calculer la pseudo-inverse
rhs = Ff - Kfe.dot(Ue)
Uf = Kff_inv.dot(rhs)


    # Solve system of equations
#Uf = np.linalg.solve(Kff, Ff - Kfe.dot(Ue))

    # Compute nodal reactions
Fe = Kee.dot(Ue) + Kfe.transpose().dot(Uf)

    # Complete global displacement and force vector
utot = np.zeros((D['nnode'] * D['ndof'], D['loadcase']))
ptot = np.zeros((D['nnode'] * D['ndof'], D['loadcase']))

for j in range(D['loadcase']):
    for i in range(D['nnode'] * D['ndof']):
        if ieqn[i] < 0:
            utot[i, j] = Ue[-int(ieqn[i]) - 1, int(j)]  # MATLAB to Python index adjustment
            ptot[i, j] = Fe[-int(ieqn[i]) - 1, int(j)]
        else:
            utot[i, j] = Uf[int(ieqn[i])-1, int(j)]
            ptot[i, j] = Ff[int(ieqn[i])-1, int(j)]







#def Fcalc(utot, ptot, D, rho, theta, phi):
    
rho = np.reshape(rho, (nely, nelx, nelz))
theta = np.reshape(theta, (nely, nelx, nelz))
phi = np.reshape(phi, (nely, nelx, nelz))
    
F = 0
dcdrho = np.zeros((D['nely'], D['nelx'], D['nelz']))
dcdtheta = np.zeros((D['nely'], D['nelx'], D['nelz']))

k = 0
for m in range(D['nelz']):
    for i in range(D['nelx']):
        for j in range(D['nely']):
            k += 1
            elemnodes = ielem[k - 1, :]  # Adjusting for 0-based indexing
                # Get element DOF number array
                
                #ieldof = get_ieldof(elemnodes, D)
            ieldof = []  # Initialize the list to hold global dofs

            for j1 in range(nenode):
                for j2 in range(ndof ):  # j2 goes from 1 to ndof inclusive
                    global_dof = (elemnodes[j1] +1)* ndof + j2 - ndof
                    ieldof.append(global_dof)
            
            c_e, dcdrho_e, dcdtheta_e = eval(f"{D['elname']}(rho[j, i, m], theta[j, i, m], phi[j, i, m], elemnodes, D, 1, utot[ieldof, :], X_X)")
            
            #eldat_F, eldat_K = eval(f"{D['elname']}(rho[i], theta[i], phi[i], elemnodes, D, -1 ,0, X_X)")
            
            # Objective function value
            F += c_e
                # Derivative w.r.t density
            dcdrho[j, i, m] += dcdrho_e
                # Derivative w.r.t theta for gradient check
            dcdtheta[j, i, m] += dcdtheta_e



#return D, ieqn,utot, ptot, F, dcdrho, dcdtheta







