import numpy as np
import scipy.io  # for loading .mat files
import os

nelx = 300                          # number of elements in x-direction
nely = 100                           # number of elements in y-direction
volfrac_max = 0.5                   # volume fraction max
penal = 3                           # penalty for SIMP interpolation
rmin = 2.5                          # filter radius


# Defining the path from which the baseline design is uploaded
pathtopbaseline = f"./baseline/{nelx}X{nely}/vol={volfrac_max}/SIMP={penal}/rmin={rmin}/"

# uplading the file here
mat = scipy.io.loadmat(os.path.join(pathtopbaseline, "baseline.mat"))
xPhys = mat['xPhys']

print(xPhys)

# saving as the npz file
varname = os.path.join(pathtopbaseline, "baseline.npz")

print(varname)

np.savez(
    varname,
    xPhys=xPhys    
)

data = np.load(varname)

xPhys1 = np.asarray(data["xPhys"])

print(xPhys1)
