import numpy as np
import matplotlib.pyplot as plt

def plot_layer(x, X_X):
    """
    Plots layer-by-layer visualization in the print plane.

    Parameters:
    x (np.ndarray): Design variables, where the first half corresponds to density 
                    (rho) and the second half corresponds to fiber angle (theta).
    """
    global ielem, nelx, nely, nelz
    
    # Extracting the element connectivity
    ielem1 = ielem[:nelx * nely, :4]

    # Reshaping density and fiber angle arrays
    rho = x[:len(x) // 2].reshape(nely, nelx, nelz)
    theta = x[len(x) // 2:].reshape(nely, nelx, nelz)

    # Extracting coordinates from global variable X
    XX = X_X
    xx = XX[:, 0]
    yy = XX[:, 1]
    zz = XX[:, 2]

    plt.figure()

    for i in range(nelz):  # Iterate through each layer in the z-direction
        ax = plt.subplot(nelz, 1, i + 1)
        rho_plot = 1 - rho[:, :, i].flatten()
        theta_plot = theta[:, :, i].flatten()

        # Plot density
        ax.tripcolor(xx[ielem1].T, yy[ielem1].T, rho_plot, edgecolor='none', cmap='gray')
        ax.set_xlim([min(xx), max(xx)])
        ax.set_ylim([min(yy), max(yy)])
        ax.set_aspect('equal')
        ax.axis('off')

        # Plot fiber angle
        center_xx = np.mean(xx[ielem1], axis=1)
        center_yy = np.mean(yy[ielem1], axis=1)

        tip1x = center_xx + 0.25 * np.cos(theta_plot)
        tip1y = center_yy + 0.25 * np.sin(theta_plot)
        tip2x = center_xx - 0.25 * np.cos(theta_plot)
        tip2y = center_yy - 0.25 * np.sin(theta_plot)

        fiber_x = np.vstack((tip2x, center_xx, tip1x))
        fiber_y = np.vstack((tip2y, center_yy, tip1y))

        ax.plot(fiber_x.flatten(), fiber_y.flatten(), color='white', linewidth=2)

    plt.show()

    # Save the solution variables into a file (uncomment if needed)
    # np.savetxt('20x10x6_optimized_data.dat', x)
