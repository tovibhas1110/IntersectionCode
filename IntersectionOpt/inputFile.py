# Here you input the variables - Domain Size, Filter radius, Filter type, Material Properties

class User:
    def __init__(self, data):
        for key, value in data.items():
            setattr(self, key, value)

    def inputFile():
        # === Design space definition ===
        nelx = 180                           # Number of elements in x-direction
        nely = 60                            # Number of elements in y-direction

        # === Allowed volume in the design domain ===
        volfrac_max = 0.5                    # Volume fraction max
    
        # === TO settings - SIMP coefficient, filter radius, filter type ===  
        penal = 3                            # Penalty for SIMP interpolation
        rmin = 2.5                           # Filter radius
        ft = 2                               # Filtering type (1: none, 2: with filtering)
    
        #=== Paramters to control the Intersection objective ===
        ratio = 3/12                         # Feature size to Intersection Filter domain size
        rmin1 = rmin / ratio                 # Size of the Intersection Filter Domain
        vt = 1                               # Penalization on the Weight Function
        wt = 0.2                             # Contribution of Compliance Objective to the Total Objective
        pI = 1                               # Penalization on intersection objective
        lambdaN = 0.02                       # Normalization factor for the total Objective

        return { "nelx" : nelx,
                 "nely" : nely, 
                 "vol" : volfrac_max, 
                 "p" : penal,
                 "rminC": rmin, 
                 "ft" : ft, 
                 "rminI": rmin1,
                 "vt": vt,
                 "wt": wt, 
                 "pI": pI, 
                 "normTotalObj":lambdaN}
