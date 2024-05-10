import matplotlib.pyplot as plt
import numpy as np

from inspect import currentframe

import os
import shutil

"""
Create Your Own Lattice Boltzmann Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate flow past cylinder
for an isothermal fluid

"""



def main(Nx, Ny, Nt):
    """ Lattice Boltzmann Simulation """
    
    # if folder exists delete it so we start fresh
    if os.path.exists("output/reference"):
        shutil.rmtree("output/reference")
    
    # create output/reference folder
    if not os.path.exists("output/reference"):
        os.mkdir("output/reference")

    # Simulation parameters
    rho0                   = 100    # average density
    tau                    = 0.6    # collision timescale
    plotRealTime = True # switch on for plotting as the simulation goes along
    
    # Lattice speeds / weights
    NL = 9
    idxs = np.arange(NL)
    cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
    cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
    weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1
    
    # Initial Conditions
    F = np.ones((Ny,Nx,NL)) #* rho0 / 
    
    # print("Line: " + str(currentframe().f_lineno) + " size of F: ", str(F.shape))
    
    np.random.seed(42)
    F += 0.01#*np.random.randn(Ny,Nx,NL)

    np.save("output/reference/F_start.npy", F)

    X, Y = np.meshgrid(range(Nx), range(Ny))

    np.save("output/reference/X_start.npy", X)
    np.save("output/reference/Y_start.npy", Y)

    # print(X)
    # print(Y)
    
    # print("Line: " + str(currentframe().f_lineno) + " size of X: ", str(X.shape))
    # print("Line: " + str(currentframe().f_lineno) + " size of Y: ", str(Y.shape))
    
    F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))

    np.save("output/reference/F_3_setup.npy", F)

    rho = np.sum(F,2)

    np.save("output/reference/rho_start.npy", rho)

    # print first row of rho    
    # print("Line: " + str(currentframe().f_lineno) + " rho: ", str(rho[0]))
    
    # print("Line: " + str(currentframe().f_lineno) + " size of rho: ", str(rho.shape))

    for i in idxs:
        F[:,:,i] *= rho0 / rho

    np.save("output/reference/F_normalized_start.npy", F)
    
    # Cylinder boundary
    X, Y = np.meshgrid(range(Nx), range(Ny))
    
    # print("Line: " + str(currentframe().f_lineno) + " size of cylinder X: ", str(X.shape))
    # print("Line: " + str(currentframe().f_lineno) + " size of cylinder Y: ", str(Y.shape))
    
    cylinder = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2

    np.save("output/reference/cylinder_start.npy", cylinder)
    np.save("output/reference/cylinder.npy", cylinder)

    # print first 10 elements of X and Y
    # print("Line: " + str(currentframe().f_lineno) + " X: ", str(X[:10]))
    # print("Line: " + str(currentframe().f_lineno) + " Y: ", str(Y[:10]))

    #plot cylinder
    #plt.imshow(cylinder, cmap='gray')
    #plt.show()

    # wait for 5sec
    #plt.pause(5)
    
    # print("Line: " + str(currentframe().f_lineno) + " size of cylinder: ", str(cylinder.shape))

    
    # Prep figure
    fig = plt.figure(figsize=(4,2), dpi=80)
    
    # Simulation Main Loop
    for it in range(Nt):
        print(f"\r{it}", end='')


        #np.save(f"output/reference/F_before_drift_{it}.npy", F)
        
        # Drift
        for i, cx, cy in zip(idxs, cxs, cys):
            F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
            F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)
        
        #np.save(f"output/reference/F_after_drift_{it}.npy", F)

        
        # Set reflective boundaries
        bndryF = F[cylinder,:]
        #np.save(f"output/reference/bndryF_start_{it}.npy", bndryF)

        #print("Line: " + str(currentframe().f_lineno) + " size of bndryF: ", str(bndryF.shape))
        bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]

        #np.save(f"output/reference/bndryF_reordered_{it}.npy", bndryF)
        #print("Line: " + str(currentframe().f_lineno) + " size of bndryF: ", str(bndryF.shape))
    
        
        # Calculate fluid variables
        rho = np.sum(F,2)

        #np.save(f"output/reference/rho_loop_{it}.npy", rho)

        #print("Line: " + str(currentframe().f_lineno) + " size of rho: ", str(rho.shape))
        ux  = np.sum(F*cxs,2) / rho

        #np.save(f"output/reference/ux_loop_{it}.npy", ux)

        #print("Line: " + str(currentframe().f_lineno) + " size of ux: ", str(ux.shape))
        uy  = np.sum(F*cys,2) / rho

        #np.save(f"output/reference/uy_loop_{it}.npy", uy)

        #print("Line: " + str(currentframe().f_lineno) + " size of uy: ", str(uy.shape))
        
        
        # Apply Collision
        Feq = np.zeros(F.shape)
        #print("Line: " + str(currentframe().f_lineno) + " size of Feq: ", str(Feq.shape))
        for i, cx, cy, w in zip(idxs, cxs, cys, weights):
            Feq[:,:,i] = rho * w * ( 1 + 3*(cx*ux+cy*uy)  + 9*(cx*ux+cy*uy)**2/2 - 3*(ux**2+uy**2)/2 )
        
        #np.save(f"output/reference/Feq_{it}.npy", Feq)

        F += -(1.0/tau) * (F - Feq)
        
        # Apply boundary 
        F[cylinder,:] = bndryF
        
        
        # plot in real time - color 1/2 particles blue, other half red
        #if (plotRealTime and (it % 10) == 0) or (it == Nt-1):
        plt.cla()
        ux[cylinder] = 0
        uy[cylinder] = 0
        vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
        #print("Line: " + str(currentframe().f_lineno) + " size of vorticity: ", str(vorticity.shape))
        
        np.save(f"output/reference/vorticity_{it:05}.npy", vorticity)

        vorticity[cylinder] = np.nan


        vorticity = np.ma.array(vorticity, mask=cylinder)

        plt.title(f"iteration: {it}")

        plt.imshow(vorticity, cmap='bwr')
        plt.imshow(~cylinder, cmap='gray', alpha=0.3)
        plt.clim(-.1, .1)
        ax = plt.gca()
        ax.invert_yaxis()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)    
        ax.set_aspect('equal')    
        plt.pause(0.001)
    print("")
    
    # Save figure
    #plt.savefig('latticeboltzmann.png',dpi=240)
    #plt.show()
        
    return 0



if __name__== "__main__":
    import sys
    if len(sys.argv) == 4:
        Nx = int(sys.argv[1])
        Ny = int(sys.argv[2])
        Nt = int(sys.argv[3])
        main(Nx, Ny, Nt)
    elif len(sys.argv) == 1:
        main(400, 100, 5000)
    else:
        print("Usage: latticeboltzmann.py Nx Ny Nt")
        sys.exit(1)

