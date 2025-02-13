import numpy as np



def elemental_3D(rho, theta, phi, lnodes, D, icode, ue, X):
   
    
    xloc = X[lnodes, :]
    #xloc = np.array([X[lnodes, 0], X[lnodes, 1], X[lnodes, 2]]).T  # Transposée pour s'assurer qu'elle est de taille (8, 3)
    #p = D['p']  # Assurez-vous que D est un dictionnaire avec la clé 'p'
    p=3
    # Propriétés de l'élément
    Ex = 7.34
    Ey = 3.43
    Gxy = 1.39
    vxy = 0.42
    vyz = 0.47
  

      # Matrice de compliance pour un matériau transitoire isotrope
    SS = np.array([
      [1/Ex, -vxy/Ex, -vxy/Ex, 0, 0, 0],
      [-vxy/Ex, 1/Ey, -vyz/Ey, 0, 0, 0],
      [-vxy/Ex, -vyz/Ey, 1/Ey, 0, 0, 0],
      [0, 0, 0, 2*(1+vyz)/Ey, 0, 0],
      [0, 0, 0, 0, 1/Gxy, 0],
      [0, 0, 0, 0, 0, 1/Gxy]
      ])
    C = np.linalg.inv(SS)  # Inverse de la matrice de compliance


    Tinv = np.array([
      [np.cos(theta)**2 * np.cos(phi)**2, np.sin(theta)**2, np.cos(theta)**2 * np.sin(phi)**2, -2 * np.cos(theta) * np.sin(theta) * np.sin(phi), np.cos(theta)**2 * np.sin(2*phi), -2 * np.cos(theta) * np.cos(phi) * np.sin(theta)],
      [np.cos(phi)**2 * np.sin(theta)**2, np.cos(theta)**2, np.sin(theta)**2 * np.sin(phi)**2, np.sin(2*theta) * np.sin(phi), np.sin(theta)**2 * np.sin(2*phi), np.cos(phi) * np.sin(2*theta)],
      [np.sin(phi)**2, 0, np.cos(phi)**2, 0, -2 * np.cos(phi) * np.sin(phi), 0],
      [-np.cos(phi) * np.sin(theta) * np.sin(phi), 0, np.cos(phi) * np.sin(theta) * np.sin(phi), np.cos(theta) * np.cos(phi), np.cos(2*phi) * np.sin(theta), -np.cos(theta) * np.sin(phi)],
      [-np.cos(theta) * np.cos(phi) * np.sin(phi), 0, np.cos(theta) * np.cos(phi) * np.sin(phi), -np.cos(phi) * np.sin(theta), np.cos(theta) * np.cos(2*phi), np.sin(theta) * np.sin(phi)],
      [np.cos(theta) * np.cos(phi)**2 * np.sin(theta), -np.cos(theta) * np.sin(theta), np.cos(theta) * np.sin(theta) * np.sin(phi)**2, np.cos(2*theta) * np.sin(phi), 2 * np.cos(theta) * np.cos(phi) * np.sin(theta) * np.sin(phi), np.cos(2*theta) * np.cos(phi)]
      ])


      # Matrices de transformation
   
    Tinvtransp = Tinv.T  # Transposée de Tinv
    DTinv = np.array([
          [-2 * np.cos(phi)**2 * np.cos(theta) * np.sin(theta), 2 * np.cos(theta) * np.sin(theta), -2 * np.cos(theta) * np.sin(phi)**2 * np.sin(theta), 2 * np.sin(phi) * np.sin(theta)**2 - 2 * np.cos(theta)**2 * np.sin(phi), -2 * np.sin(2*phi) * np.cos(theta) * np.sin(theta), 2 * np.cos(phi) * np.sin(theta)**2 - 2 * np.cos(phi) * np.cos(theta)**2],
          [2 * np.cos(phi)**2 * np.cos(theta) * np.sin(theta), -2 * np.cos(theta) * np.sin(theta), 2 * np.cos(theta) * np.sin(phi)**2 * np.sin(theta), 2 * np.cos(2*theta) * np.sin(phi), 2 * np.sin(2*phi) * np.cos(theta) * np.sin(theta), 2 * np.cos(2*theta) * np.cos(phi)],
          [0, 0, 0, 0, 0, 0],
          [-np.cos(phi) * np.cos(theta) * np.sin(phi), 0, np.cos(phi) * np.cos(theta) * np.sin(phi), -np.cos(phi) * np.sin(theta), np.cos(2*phi) * np.cos(theta), np.sin(phi) * np.sin(theta)],
          [np.cos(phi) * np.sin(phi) * np.sin(theta), 0, -np.cos(phi) * np.sin(phi) * np.sin(theta), -np.cos(phi) * np.cos(theta), -np.cos(2*phi) * np.sin(theta), np.cos(theta) * np.sin(phi)],
          [np.cos(phi)**2 * np.cos(theta)**2 - np.cos(phi)**2 * np.sin(theta)**2, np.sin(theta)**2 - np.cos(theta)**2, np.cos(theta)**2 * np.sin(phi)**2 - np.sin(phi)**2 * np.sin(theta)**2, -2 * np.sin(2*theta) * np.sin(phi), 2 * np.cos(phi) * np.cos(theta)**2 * np.sin(phi) - 2 * np.cos(phi) * np.sin(phi) * np.sin(theta)**2, -2 * np.sin(2*theta) * np.cos(phi)]
      ])
      
    DTinvtransp = DTinv.T  # Transposée de DTinv
      

      # Initialisation des matrices de l'élément
    if icode == -1:
        eldat_F = np.zeros((D['nenode'] * D['ndof'], 1))
        eldat_K = np.zeros((D['nenode'] * D['ndof'], D['nenode'] * D['ndof']))
        eldat_KK = 0
    elif icode == 1:
        eldat_F = 0
        eldat_K = 0
        eldat_KK = 0
    elif icode == 2:
        eldat_F = 0
        eldat_K = 0
      
      # Boucle sur les points de Gauss
    wij = np.array([[2, 0, 0], [1, 1, 0], [5/9, 8/9, 5/9]]).T
    xij = np.array([[0, 0, 0], [-1/np.sqrt(3), 1/np.sqrt(3), 0], [-np.sqrt(3/5), 0, np.sqrt(3/5)]]).T
  

    ngpt = 2
    for ii in range(ngpt):
        for jj in range(ngpt):
            for kk in range(ngpt):
                w = wij[ii, ngpt-1] * wij[jj, ngpt-1] * wij[kk, ngpt-1]
                xi = xij[ii, ngpt-1]
                eta = xij[jj, ngpt-1]
                zeta = xij[kk, ngpt-1]
                  
                  # Matrice dN
                dN = (1/8) * np.array([
                      [-(eta-1)*(zeta-1), (eta-1)*(zeta-1), -(eta+1)*(zeta-1), (eta+1)*(zeta-1), (eta-1)*(zeta+1), -(eta-1)*(zeta+1), (eta+1)*(zeta+1), -(eta+1)*(zeta+1)],
                      [-(xi-1)*(zeta-1), (xi+1)*(zeta-1), -(xi+1)*(zeta-1), (xi-1)*(zeta-1), (xi-1)*(zeta+1), -(xi+1)*(zeta+1), (xi+1)*(zeta+1), -(xi-1)*(zeta+1)],
                      [-(eta-1)*(xi-1), (eta-1)*(xi+1), -(eta+1)*(xi+1), (eta+1)*(xi-1), (eta-1)*(xi-1), -(eta-1)*(xi+1), (eta+1)*(xi+1), -(eta+1)*(xi-1)]
                  ])
                  
                  
                J = dN @ xloc
                detj = np.linalg.det(J)
                  
                  # Define strain-displacement matrix
                B = np.linalg.solve(J, dN)
                B = np.array([
                          [B[0, 0], 0, 0, B[0, 1], 0, 0, B[0, 2], 0, 0, B[0, 3], 0, 0, B[0, 4], 0, 0, B[0, 5], 0, 0, B[0, 6], 0, 0, B[0, 7], 0, 0],
                          [0, B[1, 0], 0, 0, B[1, 1], 0, 0, B[1, 2], 0, 0, B[1, 3], 0, 0, B[1, 4], 0, 0, B[1, 5], 0, 0, B[1, 6], 0, 0, B[1, 7], 0],
                          [0, 0, B[2, 0], 0, 0, B[2, 1], 0, 0, B[2, 2], 0, 0, B[2, 3], 0, 0, B[2, 4], 0, 0, B[2, 5], 0, 0, B[2, 6], 0, 0, B[2, 7]],
                          [0, B[2, 0], B[1, 0], 0, B[2, 1], B[1, 1], 0, B[2, 2], B[1, 2], 0, B[2, 3], B[1, 3], 0, B[2, 4], B[1, 4], 0, B[2, 5], B[1, 5], 0, B[2, 6], B[1, 6], 0, B[2, 7], B[1, 7]],
                          [B[2, 0], 0, B[0, 0], B[2, 1], 0, B[0, 1], B[2, 2], 0, B[0, 2], B[2, 3], 0, B[0, 3], B[2, 4], 0, B[0, 4], B[2, 5], 0, B[0, 5], B[2, 6], 0, B[0, 6], B[2, 7], 0, B[0, 7]],
                          [B[1, 0], B[0, 0], 0, B[1, 1], B[0, 1], 0, B[1, 2], B[0, 2], 0, B[1, 3], B[0, 3], 0, B[1, 4], B[0, 4], 0, B[1, 5], B[0, 5], 0, B[1, 6], B[0, 6], 0, B[1, 7], B[0, 7], 0]
                          ])

                  # Calcul de la matrice de rigidité ke
                ke = np.dot(np.dot(B.T, np.dot(Tinv, C)), np.dot(Tinvtransp, B))
                ke = (ke + ke.T) / 2
                 
                  # Calcul des matrices selon icode
                if icode == -1:
                    eldat_K += rho**p * ke * detj * w
                elif icode == 1:
                  # Compute c at element level
                    eldat_F += np.trace(rho**p * ue.T @ ke @ ue * detj * w)
                  # Compute dc at element level w.r.t density
                    eldat_K -= np.trace(p * rho**(p-1) * ue.T @ ke @ ue * detj * w)  # dc/drho
                  # 1st derivative w.r.t theta
                    eldat_KK -= np.trace(rho**p * ue.T @ (B.T @ (DTinv @ C @ Tinvtransp + Tinv @ C @ DTinvtransp) @ B) @ ue * detj * w)  # dc/dtheta





  #•elif icode == 2:
                      #eldat_KK -= np.trace(rho**p * np.dot(np.dot(ue.T, B.T), np.dot(np.dot(DTinv, C), Tinvtransp) + np.dot(Tinv, np.dot(C, DTinvtransp)), B) * ue) * detj * w

  #return eldat_F, eldat_K, eldat_KK    
      
    # Matrice de compliance pour un matériau transitoire isotrope
  
    # Boucle sur les points de Gauss
  
    
                
 
    return eldat_F, eldat_K, eldat_KK    
 
    
 
    
 
    
 
  