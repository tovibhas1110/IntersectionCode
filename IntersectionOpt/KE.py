import numpy as np


# QUAD Element
###  4 _______________3
###   |               |
###   |               |
###   |               |
###   |               |
### 1 |_______________|2
### Element level displacement array {u1,v1,u2,v2,u3,v3,u4,v4}^T

# === Material properties Definition===
E0 = 1                          # Young's modulus
Emin = E0 * 1e-9                     # void stiffness
nu = 0.3                             # Poisson's ratio


def KE(**kwargs):
    
    EleType = kwargs['EleType']
    IntType = kwargs['IntType']

    if EleType == 'QUAD' and IntType == 'Full':
        
        GP = np.sqrt(1/3)

        KE = (Bmat(-GP, -GP).T @ Dmat(E0, Emin, nu) @ Bmat(-GP, -GP) +
        Bmat(GP, -GP).T  @ Dmat(E0, Emin, nu)  @ Bmat(GP, -GP)  +
        Bmat(GP, GP).T   @ Dmat(E0, Emin, nu) @ Bmat(GP, GP)   +
        Bmat(-GP, GP).T  @ Dmat(E0, Emin, nu) @ Bmat(-GP, GP))

        print(KE)

    return KE 

def Bmat(x, y):
    B = 0.25 * np.array([
        [-(1 - y), 0, (1 - y), 0, (1 + y), 0, -(1 + y), 0],
        [0, -(1 - x), 0, -(1 + x), 0, (1 + x), 0, (1 - x)],
        [-(1 - x), -(1 - y), -(1 + x), (1 - y), (1 + x), (1 + y), (1 - x), -(1 + y)]
    ])
    return B

def Dmat(E0, Emin, nu):
    
    factor = (E0 - Emin) / (1 - nu**2)
    D = factor * np.array([
        [1,    nu,        0],
        [nu,   1,         0],
        [0,    0, 0.5*(1 - nu)]
    ])
    return D


