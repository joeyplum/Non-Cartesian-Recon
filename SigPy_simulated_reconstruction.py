# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 08:50:28 2022

@author: PLUMD2
"""
# Module install:
# (optional for sigpy) conda install sigpy
# (optional for plot support) conda install matplotlib
# (optional for CUDA support) conda install cupy
# (optional for MPI support) conda install mpi4py

# Import modules:
import numpy as np
import sigpy as sp
import sigpy.plot as pl
import sigpy.mri as mr
import matplotlib.pyplot as plt

# Set up simulated phantom:
phant = sp.shepp_logan(((128,128)))   

# Visualize:
pl.ImagePlot(abs(phant), 
              title = "Simulated Phantom")

# Set up spiral saplting:
n_spirals = 15
traj = mr.spiral(fov = 0.35, # field of view in meters.
                      N = len(phant)/2, # effective matrix shape.
                      f_sampling = 0.25, # undersaplting factor in freq encoding direction.
                      R = 0.25, # undersaplting factor.
                      ninterleaves = n_spirals, # number of spiral interleaves.
                      alpha = 2, # variable density factor.
                      gm = 32.5E-3, # maximum gradient apltitude (T/m).
                      sm = 150, # maximum slew rate (T/m/s).
                      gamma = 73.997E6) # gyromagnetic ratio in rad/T/s.
traj = np.reshape(traj, (n_spirals,int(len(traj)/n_spirals),2))
pl.ScatterPlot(traj, 
                title = "Spiral Sampling Coordinates")


# Account for sample density:
dcf = mr.pipe_menon_dcf(traj)
plt.plot(dcf[1,...], '-')
plt.ylabel('Sample number')
plt.xlabel('DCF')
plt.title('Sample Density Compensation')

# Sample image using NUFFT:
ksamp = sp.nufft(phant,traj,oversamp = 1.25, width = 4)

# Reconstruct image using NUFFT_adjoint:
recon_phant = sp.nufft_adjoint(ksamp*dcf,traj,oversamp = 1.25, width = 4)

# Visualize reconstructed phantom:
pl.ImagePlot(abs(recon_phant), 
             title = "Reconstructed Phantom with Spiral Sampling")

# Reconstruct using linear operators:
NUFFT_linop = sp.linop.NUFFT(np.shape(phant),
                             traj, 
                             oversamp = 1.25, 
                             width = 4)
NUFFT_adjoint_linop = sp.linop.NUFFTAdjoint(np.shape(phant),
                                            traj, 
                                            oversamp = 1.25, 
                                            width = 4)
ksamp = NUFFT_linop * phant
recon_phant = NUFFT_adjoint_linop * (dcf * ksamp)

# Visualize reconstructed phantom:
pl.ImagePlot(abs(recon_phant), 
             title = "Linear Operation Reconstructed Phantom with Spiral Sampling")

# # Summary figure:
# fig = plt.figure(figsize=(15, 3))
# rows = 1
# columns = 3

# fig.add_subplot(rows, columns, 1)
# plt.imshow(abs(phant), cmap="gray", origin="lower")
# plt.axis('off')
# plt.title("Simulated Phantom")

# fig.add_subplot(rows, columns, 2)
# plt.imshow(abs(recon_phant), cmap="gray", origin="lower")
# plt.axis('off')
# plt.title("Reconstructed Phantom with Spiral Sampling")

# fig.add_subplot(rows, columns, 3)
# plt.imshow(abs(recon_phant), cmap="gray", origin="lower")
# plt.axis('off')
# plt.title("Linear Operation Reconstructed Phantom with Spiral Sampling")






# Set up radial sampling: 
traj = mr.radial(coord_shape = (256,96,2), # e.g. (512, 96, 2) = 512 proj, 96 sapltes, 2 dimensions
                  img_shape = (128,128),
                  golden = True)
pl.ScatterPlot(traj,
                title = "Radial Sampling Coordinates")

# Account for sample density:
dcf = mr.pipe_menon_dcf(traj)
plt.plot(dcf[1,...], '-')
plt.ylabel('Sample number')
plt.xlabel('DCF')
plt.title('Sample Density Compensation')

    
# Saplte image using NUFFT:
ksamp = sp.nufft(phant,traj,oversamp = 1.25, width = 4)

# Reconstruct image using NUFFT_adjoint:
recon_phant = sp.nufft_adjoint(ksamp*dcf,traj,oversamp = 1.25, width = 4)

# Visualize reconstructed phantom:
pl.ImagePlot(abs(recon_phant), 
             title = "Reconstructed Phantom with Radial Sampling")

# Reconstruct using linear operators:
NUFFT_linop = sp.linop.NUFFT(np.shape(phant),
                             traj, 
                             oversamp = 1.25, 
                             width = 4)
NUFFT_adjoint_linop = sp.linop.NUFFTAdjoint(np.shape(phant),
                                            traj, 
                                            oversamp = 1.25, 
                                            width = 4)
ksamp = NUFFT_linop * phant
recon_phant = NUFFT_adjoint_linop * (dcf * ksamp)

# Visualize reconstructed phantom:
pl.ImagePlot(abs(recon_phant), 
             title = "Linear Operation Reconstructed Phantom with Radial Sampling")




# # Perform compressed sensing reconstruction:
#     # First, simulate a sensitivity map estimate.
# sens_map = mr.app.EspiritCalib(ksamp).run()
# sens_map = np.ones((100,100))
# pl.ImagePlot(abs(sens_map), 
#              title = "Estimated Sensitivity Map")

