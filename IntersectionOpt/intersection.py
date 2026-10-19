import numpy as np
import os
from scipy.sparse import lil_matrix
from scipy.sparse import coo_matrix
import scipy.io  # for loading .mat files
from scipy.sparse.linalg import spsolve
from mmasub import mmasub
import pdb
import matplotlib.pyplot as plt

# === Design space and ring domain definition - Fixed ===
# === Structured grid ===
#  ___________________________
# |___|___|___|___|___|___|___| }
# |___|___|___|___|___|___|___| > nely
# |___|___|___|___|___|___|___| }
# <------------nelx----------->

nelx = 180                          # number of elements in x-direction
nely = 60                           # number of elements in y-direction
volfrac_max = 0.5                   # volume fraction max
volfrac_min = 0.0                   # volume fraction min
penal = 3                           # penalty for SIMP interpolation
rmin = 2.5                          # filter radius
ft = 2                              # filtering type (1: none, 2: with filtering)
pI = 1                              # penalization on intersection objective
lambdaN = 0.02                      # normalization factor
ndof = 2 * (nely + 1) * (nelx + 1)  # total number of DOFs

# === Material properties ===
E0 = 210000                   # Young's modulus
Emin = E0 * 1e-9              # void stiffness
nu = 0.3                      # Poisson's ratio

# === Intersection Objective Parameters ===
ratio = 3/12            # define the intersection radius and the density filter radius
rmin1 = rmin / ratio    # Intersection Radius
vt = 1                  # Default parameter === do not change ===
wt = 0.2                # Weight on the compliance and intersection objective - reduce the weight to make the intersection prominent
IepsF = 1               # Default paramter === do not change ===

# === Define path to save the variables === 
pathtop = (
    f"./result/{nelx}X{nely}/"
    f"vol={volfrac_max}/"
    f"SIMP={penal}/"
    f"rmin={rmin}/"
    f"ratio={ratio}/"
    f"vt={vt}/"
    f"wt={wt}/"
    f"IepsF={IepsF}/"
)

# === Create the directory ===
os.makedirs(pathtop, exist_ok=True)

# === Setiing up the finite element analysis ===
# === Define matrices ===
A11 = np.array([[12,  3, -6, -3],
                [ 3, 12,  3,  0],
                [-6,  3, 12, -3],
                [-3,  0, -3, 12]])

A12 = np.array([[-6, -3,  0,  3],
                [-3, -6, -3, -6],
                [ 0, -3, -6,  3],
                [ 3, -6,  3, -6]])

B11 = np.array([[-4,  3, -2,  9],
                [ 3, -4, -9,  4],
                [-2, -9, -4, -3],
                [ 9,  4, -3, -4]])

B12 = np.array([[ 2, -3,  4, -9],
                [-3,  2,  9, -2],
                [ 4,  9,  2,  3],
                [-9, -2,  3,  2]])

# === Define element stiffness matrix ===
KE = 1/(1-nu**2)/24 * (np.block([[A11, A12],
                                  [A12.T, A11]]) + nu * np.block([[B11, B12],
                                                                  [B12.T, B11]]))
KE_vec = np.reshape(KE,(8*8,1),order='F')

# === Element number ====
ElementNumber = np.reshape(np.arange(1, nelx * nely + 1), (nely, nelx), order='F')

# === Node numbering ===
nodenrs = np.reshape(np.arange(1, (1+nelx)*(1+nely)+1),(nely+1, nelx+1), order='F')

# === Edof vector ===
# === Quad Element left lower first dof ====
edofVec = np.reshape((2*nodenrs[:-1, :-1] + 1),(nely*nelx, 1), order='F')

# === offset array generation to generate the other dofs ===
offset = np.array([0, 1, 2*nely + 2, 2*nely + 3, 2*nely, 2*nely + 1, -2, -1])

# === first dof arranged in an array ====
edofVec_tiled = np.tile(edofVec, (1, 8))  # shape (nelx*nely, 8)

# === offset is tiled in an array ====
offset_tiled = np.tile(offset, (nelx*nely, 1))  # shape (nelx*nely, 8)

# === to create the dofs array for each element ===
edofMat = edofVec_tiled + offset_tiled

# ===  iK: repeat each row of edofMat 8 times column-wise ===
iK = np.reshape(np.kron(edofMat, np.ones((8, 1))).T.astype(int),(64*nely*nelx, 1),order='F')

# === jK: repeat each row of edofMat 8 times row-wise ===
jK = np.reshape(np.kron(edofMat, np.ones((1, 8))).T.astype(int),(64*nely*nelx, 1),order='F')

# === Define force vector (sparse) ===
F = lil_matrix((ndof, 1))
F[2*(nelx+1)*(nely+1)-1, 0] = -1000  # MATLAB indexing starts at 1, Python at 0

# === Define displacement vector ===
U = np.zeros((ndof, 1))

# === Define fixed DOFs ===
fixeddofs = np.arange(0, 2*(nely+1))

# === All DOFs ===
alldofs = np.arange(ndof)

# === Free DOFs ===
freedofs = np.setdiff1d(alldofs, fixeddofs)

# === Defining the density filter ===
size = nelx * nely * (2 * (int(np.ceil(rmin)) - 1) + 1)**2

iH = np.zeros(size, dtype=int)
jH = np.zeros(size, dtype=int)
sH = np.zeros(size, dtype=float)
k = 0

for i1 in range(nelx):
    for j1 in range(nely):
        e1 = i1 * nely + j1
        i2_min = max(i1 - (int(np.ceil(rmin)) - 1), 0)
        i2_max = min(i1 + (int(np.ceil(rmin)) - 1), nelx - 1)
        j2_min = max(j1 - (int(np.ceil(rmin)) - 1), 0)
        j2_max = min(j1 + (int(np.ceil(rmin)) - 1), nely - 1)
        
        for i2 in range(i2_min, i2_max + 1):
            for j2 in range(j2_min, j2_max + 1):
                e2 = i2 * nely + j2
                iH[k] = e1
                jH[k] = e2
                sH[k] = max(0, rmin - np.sqrt((i1 - i2)**2 + (j1 - j2)**2))
                k += 1

iH = iH[:k]
jH = jH[:k]
sH = sH[:k]
H = coo_matrix((sH, (iH, jH)), shape=(nelx * nely, nelx * nely))
Hs = np.array(H.sum(axis=1))

# === Defining the intersection filter ===
size = nelx * nely * (2 * (int(np.ceil(rmin1)) - 1) + 1)**2

iDI = np.zeros(size, dtype=int)
jDI = np.zeros(size, dtype=int)
sDI = np.zeros(size, dtype=float)
k = 0

for i1 in range(nelx):
    for j1 in range(nely):
        e1 = i1 * nely + j1
        i2_min = max(i1 - (int(np.ceil(rmin1)) - 1), 0)
        i2_max = min(i1 + (int(np.ceil(rmin1)) - 1), nelx - 1)
        j2_min = max(j1 - (int(np.ceil(rmin1)) - 1), 0)
        j2_max = min(j1 + (int(np.ceil(rmin1)) - 1), nely - 1)
        
        for i2 in range(i2_min, i2_max + 1):
            for j2 in range(j2_min, j2_max + 1):
                e2 = i2 * nely + j2
                iDI[k] = e1
                jDI[k] = e2
                sDI[k] = max(0, rmin1 - np.sqrt((i1 - i2)**2 + (j1 - j2)**2))
                k += 1

iDI = iDI[:k]
jDI = jDI[:k]
sDI = sDI[:k]
DI = coo_matrix((sDI, (iDI, jDI)), shape=(nelx * nely, nelx * nely))
DIs = np.array(DI.sum(axis=1))

# === To calculate the weights for the compliance and intersection objective ===
# === Defining the path from which the baseline design is uploaded ===
pathtopbaseline = f"./baseline/{nelx}X{nely}/vol={volfrac_max}/SIMP={penal}/rmin={rmin}/"

# === Uploading the variables ===
mat = np.load(os.path.join(pathtopbaseline, "baseline.npz"))

# === Defining the design variable ===
xPhys = np.asarray(mat["xPhys"])  # This will be a 2D array or possibly a MATLAB column vector

# === Defining the vector ===
xPhys_vec = (np.reshape(xPhys,(nelx*nely,1),order="F"))  # MATLAB column-major flatten

# === Apply intersection filter ===
xPhys1_vec = (DI.dot(xPhys_vec)) / DIs

# === Intersection size calculation ===
xPhys1 = np.reshape(xPhys1_vec, (nely, nelx), order='F')

# === Creating the system level stiffness matrix ===
sK = np.reshape(np.kron(Emin + (xPhys_vec.T)**penal * (E0 - Emin), KE_vec),(64*nelx*nely,1),order = 'F')
K = coo_matrix((sK.flatten(), (iK.flatten()-1, jK.flatten()-1)), shape=(ndof, ndof))
K = (K + K.T) / 2

# === Solving for free DOFs ===
U[freedofs] = np.reshape(spsolve(K[freedofs, :][:, freedofs], F[freedofs]),(len(freedofs),1),order='F')

# === Calculating compliance ===
ce = np.reshape(np.sum((np.squeeze(U[edofMat-1]) @ KE) * np.squeeze(U[edofMat-1]), axis=1),(nely,nelx),order='F')
c = np.sum((Emin + xPhys**penal * (E0 - Emin)) * ce)

# === B-matrix - Strain-Displacement relation ===
B = np.array([
    [-0.5, 0, 0.5, 0, 0.5, 0, -0.5, 0],
    [0, -0.5, 0, -0.5, 0, 0.5, 0, 0.5],
    [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5]
])

# === D-matrix - Strain-Stress Realtionship - Plane stress condition ===
D = ((E0 - Emin) / (1 - nu**2)) * np.array([
    [1, nu, 0],
    [nu, 1, 0],
    [0, 0, 0.5 * (1 - nu)]
])

# === Intersection Objective calculation procedure ===
# === Stress calculation ===
Ue = np.squeeze(U[edofMat-1]).T
Istress = (xPhys_vec**penal) * (D @ B @ Ue).T

# === Principal stress calculation ===
Istress1 = 0.5 * (Istress[:, 0] + Istress[:, 1]) + \
           np.sqrt(0.25 * (Istress[:, 0] - Istress[:, 1])**2 + Istress[:, 2]**2)
Istress2 = 0.5 * (Istress[:, 0] + Istress[:, 1]) - \
           np.sqrt(0.25 * (Istress[:, 0] - Istress[:, 1])**2 + Istress[:, 2]**2)
Istressp = np.column_stack((Istress1, Istress2))

# === Von Mises stress calculation ===
IstressVM = np.reshape(np.sqrt(0.5 * ((Istress1 - Istress2)**2 + Istress1**2 + Istress2**2)),(nely*nelx,1),order ='F')

# === Threshold stress calculation ===
mask = (xPhys_vec > 0.1) & (xPhys_vec < 0.5)
IstressVMStr = np.sort(IstressVM[mask])[::-1]
Percent90VM = IstressVMStr[int(np.floor(0.1 * len(IstressVMStr)))-1]
Iepsilon = (Percent90VM / IepsF)**2
Iepsilon1 = Iepsilon
Istresspabs = Istressp**2 + Iepsilon1

# === Stress ratio ===
Istresspratio = np.column_stack((
    Istresspabs[:, 0] / Istresspabs[:, 1],
    Istresspabs[:, 1] / Istresspabs[:, 0]
))

# === Intersection Indicator ===
T = np.reshape(1 / np.log2(Istresspratio[:, 0] + Istresspratio[:, 1]),(nelx*nely,1),order='F')
O = (xPhys_vec**penal) * T
Oavg = DI @ (O / DIs)

# === Weight based on the intersection size ===
w = (1 - xPhys1_vec**vt)**(1 / vt)

# === Intersection Objective ===
I = (np.sum((w * O)**pI)**(1 / pI)) / np.mean(xPhys.flatten())

# === Normalization weights for the objective functions ===
wcn = c
win = I

# === Initialize design variables ===
# === wt is less than 1 then  initialization happens ===
if wt < 1:
    x = volfrac_max * np.ones((nely, nelx))
    x_vec = np.reshape(x,(nely*nelx,1),order='F')
else:
    x = baseline_data["xPhys"]
    x_vec = np.reshape(x,(nely*nelx,1),order='F')

# === Apply density filtering ===
xPhys = np.reshape(((H @ x_vec) / Hs),(nely,nelx),order='F')
xPhys_vec = np.reshape(xPhys,(nely*nelx,1),order='F')
xPhys1_vec = ((DI @ xPhys_vec) / DIs)
xPhys1 = np.reshape(xPhys1_vec,(nely,nelx),order='F')

# === Initialize other optimization variables ===
xHeavi = x.copy()
loop = 0
change = 1
m = 1
n = nelx * nely
xmin = np.reshape(1e-6 * np.ones(n),(nelx*nely,1),order='F')
xmax = np.reshape(np.ones(n),(nelx*nely,1),order='F')
a0mma = np.ones(1)
amma = np.reshape(np.zeros(m),(m,1),order='F')
cmma = 100000 * np.reshape(np.ones(m),(m,1),order='F')
dmma = np.reshape(np.zeros(m),(m,1),order='F')
move = 0.2

# === Previous designs ===
xold1 = volfrac_max * np.ones((nely, nelx))
xold1_vec = np.reshape(xold1,(nely*nelx,1),order = 'F')
xold2 = volfrac_max * np.ones((nely, nelx))
xold2_vec = np.reshape(xold2,(nely*nelx,1),order = 'F')

# === MMA bounds ===
low = 1e-6 * np.reshape(np.ones(n),(nelx*nely,1),order='F')
upp = np.reshape(np.ones(n),(nelx*nely,1),order='F')

# === Defining storing variable ===
ObjectiveC = np.empty((0, 2))
ObjectiveI = np.empty((0, 2))
TotalObjective = np.empty((0, 2))
Volume = np.empty((0, 2))

# === Starting the optimization loop ===
while loop < 200:
    
    loop += 1

    # === Computing the global stiffness matrix ===
    sK = np.reshape(np.kron(Emin + (xPhys_vec.T)**penal * (E0 - Emin), KE_vec),(64*nelx*nely,1),order = 'F')
    K = coo_matrix((sK.flatten(), (iK.flatten()-1, jK.flatten()-1)), shape=(ndof, ndof))
    K = (K + K.T) / 2

    # === Solving the system equation K(xPhys)u = F ===
    U[freedofs] = np.reshape(spsolve(K[freedofs, :][:, freedofs], F[freedofs]),(len(freedofs),1),order='F')
    
    # === Calculating compliance objective ===
    ce = np.reshape(np.sum((np.squeeze(U[edofMat-1]) @ KE) * np.squeeze(U[edofMat-1]), axis=1),(nely,nelx),order='F')
    c = np.sum((Emin + xPhys**penal * (E0 - Emin)) * ce)

    # === Sensitivity calculation - compliance objective ===
    dc = -penal * (E0 - Emin) * (xPhys ** (penal - 1)) * ce
    dc_vec = np.reshape(dc, (nely*nelx,1), order='F')

    # === Sensitivity calculation - volume constraint ===
    dv_max = np.ones((nely,nelx))
    dv_max_vec = np.reshape(dv_max, (nely*nelx,1), order='F')

    # === ft = 1 - sensitivity filtering and 2 is density filtering ===
    if ft == 1:
        dc = (H @ (x_vec * dc_vec)) / Hs / np.maximum(1e-3, x_vec)
    elif ft == 2:
        dc_vec = (H @ (dc_vec / Hs))
        dc = np.reshape(dc_vec,(nely,nelx),order = 'F')
        dv_max_vec = (H @ (dv_max_vec / Hs))
        dv_max = np.reshape(dv_max_vec,(nely,nelx),order = 'F')
        
    # === Stress calculation ===
    Ue = np.squeeze(U[edofMat-1]).T              
    Istress = (xPhys_vec**penal) * (D @ B @ Ue).T

    # === Principal stress calculation ===
    Istress1 = 0.5 * (Istress[:, 0] + Istress[:, 1]) + \
           np.sqrt(0.25 * (Istress[:, 0] - Istress[:, 1])**2 + Istress[:, 2]**2)
    Istress2 = 0.5 * (Istress[:, 0] + Istress[:, 1]) - \
           np.sqrt(0.25 * (Istress[:, 0] - Istress[:, 1])**2 + Istress[:, 2]**2)
    Istressp = np.column_stack((Istress1, Istress2))
    
    # === Von Mises stress ===
    IstressVM = np.reshape(np.sqrt(0.5 * ((Istress1 - Istress2)**2 + Istress1**2 + Istress2**2)),(nely*nelx,1),order ='F')

    # === Threshold stress cansideration ===
    Istresspabs = Istressp**2 + Iepsilon1

    # === Stress ratio ===
    Istresspratio = np.column_stack((
        Istresspabs[:, 0] / Istresspabs[:, 1],
        Istresspabs[:, 1] / Istresspabs[:, 0]
        ))

    # === Intersection indicator ===
    T = np.reshape(1 / np.log2(Istresspratio[:, 0] + Istresspratio[:, 1]),(nelx*nely,1),order='F')
    O = (xPhys_vec**penal) * T
    Oavg = DI @ (O / DIs)

    # === Weight calcualtion based on intersection size ===
    w = (1 - xPhys1_vec**vt)**(1 / vt)

    # === Intersection Objective ===
    I = (np.sum((w * O)**pI)**(1 / pI)) / np.mean(xPhys.flatten())

    # === Adjoint sensitivity calculation ===
    
    dTdsp1 = (2 * (Istress1 / (Istress2**2 + Iepsilon))) - ((2 * Istress1 * (Istress2**2 + Iepsilon)) / ((Istress1**2 + Iepsilon)**2))
    dTdsp2 = (2 * (Istress2 / (Istress1**2 + Iepsilon))) - ((2 * Istress2 * (Istress1**2 + Iepsilon)) / ((Istress2**2 + Iepsilon)**2))
    dTdsp = np.column_stack([dTdsp1, dTdsp2])

    diff = Istress[:, 0] - Istress[:, 1]
    denom = np.sqrt(0.25 * diff**2 + Istress[:, 2]**2)

    dsp1ds1 = 0.5 + 0.25 * (diff / denom)
    dsp1ds2 = 0.5 + 0.25 * ((Istress[:, 1] - Istress[:, 0]) / denom)
    dsp1ds12 = Istress[:, 2] / denom

    dsp2ds1 = 0.5 - 0.25 * (diff / denom)
    dsp2ds2 = 0.5 - 0.25 * ((Istress[:, 1] - Istress[:, 0]) / denom)
    dsp2ds12 = - Istress[:, 2] / denom

    dIt1 = penal * (xPhys_vec ** (penal - 1)) * w * np.reshape((1.0 / np.log2(Istresspratio[:, 0] + Istresspratio[:, 1])),(nely*nelx,1),order='F')
    
    sum_part = np.stack([
        dTdsp1 * dsp1ds1 + dTdsp2 * dsp2ds1,
        dTdsp1 * dsp1ds2 + dTdsp2 * dsp2ds2,
        dTdsp1 * dsp1ds12 + dTdsp2 * dsp2ds12
        ], axis=-1)

    DBU = (D @ (B @ np.squeeze(U[edofMat.T-1]))).T 
    multiplied = sum_part * DBU 
    sum_multiplied = np.sum(multiplied, axis=1)

    dIt21 = -(np.log(2)) * np.reshape((1.0 / (Istresspratio[:, 0] + Istresspratio[:, 1])),(nelx*nely,1),order = 'F') * \
        np.reshape(((1.0 / np.log(Istresspratio[:, 0] + Istresspratio[:, 1])) ** 2),(nelx*nely,1),order ='F')* \
        (penal * (xPhys_vec ** (2 * penal - 1)) * w) * \
        np.reshape(sum_multiplied,(nelx*nely,1),order='F')

    adj11 = B.T @ D.T @ np.array([
    (dTdsp1 * dsp1ds1 + dTdsp2 * dsp2ds1),
    (dTdsp1 * dsp1ds2 + dTdsp2 * dsp2ds2),
    (dTdsp1 * dsp1ds12 + dTdsp2 * dsp2ds12)])

    adj12 = (
    np.log(2)
    * np.reshape((1.0 / (Istresspratio[:, 0] + Istresspratio[:, 1])),(nelx*nely,1),order='F')
    * np.reshape((1.0 / np.log(Istresspratio[:, 0] + Istresspratio[:, 1]))**2,(nelx*nely,1),order='F')
    * w
    * (xPhys_vec**(2 * penal))
    * adj11.T)

    adj13 = coo_matrix((adj12.T.reshape(-1), (edofMat.T.reshape(-1)-1, np.ones(len(edofMat.T.reshape(-1))).astype(int)-1)))
    adj13 = adj13.tocsr()

    # === Solving adjoint equation ====
    lambda22 = np.zeros((2 * (nelx + 1) * (nely + 1), 1))
    lambda22[freedofs] = np.reshape(spsolve(K[freedofs, :][:, freedofs], adj13[freedofs]),(len(freedofs),1),order='F')
    dIt22t = (penal * (E0 - Emin) * (xPhys ** (penal - 1)))* np.reshape(np.sum((np.squeeze(lambda22[edofMat-1]) @ KE) * np.squeeze(U[edofMat-1]), axis=1),(nely,nelx),order='F')
    dIt22 = np.reshape(dIt22t,(nelx*nely,1),order='F')

    term_inner = (
    dIt1 + dIt21 + dIt22- (DI* ((((1 - xPhys_vec**vt)**((1 / vt) - 1))
                * (xPhys1_vec**(vt - 1)))* O)/ DIs))

    dI = (nelx * nely) * ((np.sum(xPhys_vec) * term_inner
            - (np.sum((w * O)**pI)**(1 / pI)))/ (np.sum(xPhys_vec)**2))

    # === Apply density filter ====
    dI = np.reshape((H @ (dI / Hs)),(nely,nelx),order='F')  
    # ===== Total calulation ends here =========================

    # === Optimization ===============================================
    # === Multi-Objective definition - Complaince and Intersection ===
    if wt < 2 and loop > 0:
        f0val = (wt) * (c / (wcn * lambdaN)) + (1 - wt) * (I / (win * lambdaN))
    
        df0dx = np.reshape((
            (wt) * (dc.T.reshape(-1) / (wcn * lambdaN))
            + (1 - wt) * (dI.T.reshape(-1) / (win * lambdaN))
            ),(nelx*nely,1),order ='F')
    
        df0dx2 = 0 * np.reshape(df0dx,(nelx*nely,1),order ='F')
    else:
        f0val = c
        df0dx = np.reshape(dc.T.reshape(-1),(nelx*nely,1),order ='F')
        df0dx2 = 0 * np.reshape(df0dx,(nelx*nely,1),order ='F')

    # === Volume contraint definition ===
    fval01 = np.sum(xPhys.T.reshape(-1)) / (n * volfrac_max) - 1
    df01dx = np.reshape(dv_max.T.reshape(-1) / (nelx * nely * volfrac_max),(nely*nelx,1),order='F')
    df01dx2 = 0 * np.reshape(df01dx,(nelx*nely,1),order='F')
    
    # === Preparation for MMA optimizer ===
    fval = np.array([fval01])        # vertical concatenation
    dfdx = np.reshape(df01dx.T.reshape(-1),(1,n),order='F')   # transpose and vertical concatenation
    dfdx2 = 0*dfdx               # array of zeros with same shape as dfdx

    xval = x_vec

    xmax = np.minimum(1, x_vec + move)
    xmin = np.maximum(0, x_vec - move)
    
    xmma, ymma, zmma, lam, xsi, eta, mu1, zet, s1, low, upp = mmasub(
        m, n, loop, xval, xmin, xmax, xold1_vec, xold2_vec,
        f0val, df0dx, df0dx2, fval, dfdx, dfdx2,low, upp, a0mma, amma, cmma, dmma)

    xold2_vec = xold1_vec
    xold2 = np.reshape(xold2_vec,(nely,nelx),order = 'F')
    xold1_vec = xval
    xold1 = np.reshape(xold1_vec,(nely,nelx),order = 'F')
    xnew = np.reshape(xmma,(nely,nelx),order='F')
    xnew_vec = np.reshape(xnew,(nely*nelx,1),order = 'F')
    #print(xold1)
    #print(xval)
    
    #pdb.set_trace()

    if ft == 1:
        xPhys = xnew.copy()
    elif ft == 2:
        xPhys_vec = (H @ xnew_vec) / Hs
        xPhys = np.reshape(xPhys_vec,(nely,nelx),order = 'F')
        xPhys1_vec = (DI @ xPhys_vec) / DIs
        xPhys1 = np.reshape(xPhys1_vec,(nely,nelx),order = 'F')
        xPhys_old_vec = (H @ xval) / Hs
    
    if loop < 0:
        change = 1
    else:
        change = np.max(np.abs(xPhys_vec - xPhys_old_vec))

    x = xnew
    x_vec = np.reshape(x,(nely*nelx,1),order='F')
    
    ObjectiveC = np.vstack([ObjectiveC, [loop, c]])
    ObjectiveI = np.vstack([ObjectiveI, [loop, I]])
    TotalObjective = np.vstack([TotalObjective, [loop, f0val]])
    Volume = np.vstack([Volume, [loop, np.mean(xPhys)]])
    
    print("itr: %.2f, c: %.2f, I: %.2f, Total: %.2f, Vol: %.2f, Ch: %.2f" % (loop, c, I, f0val, np.mean(xPhys),change))

    plt.figure(1)
    plt.imshow(1 - xPhys, cmap='gray', vmin=0, vmax=1, aspect='equal')
    plt.axis('off')  # hides axes
    plt.tight_layout()  # similar to axis tight
    plt.savefig(f"./fig/output_{loop}.png", bbox_inches='tight', pad_inches=0)
