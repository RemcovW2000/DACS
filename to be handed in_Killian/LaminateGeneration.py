import numpy as np
import copy
import uncertainties
from uncertainties import ufloat

np.set_printoptions(linewidth=300, precision = 3)

class Material:
    def __init__(self, material, E1, E2, G12, v12, v21, strengths, v21_fiber, E1_fiber):
        self.material = material             # string
        self.E1        = E1                  # MPa
        self.E2        = E2                  # MPa
        self.G12       = G12                 # MPa
        self.v12       = v12                 # -
        self.v21       = self.u_get_v21() #self.get_v21()      # - 
        self.strengths = strengths # np.array([Xt, Xc, Yt, Yc, S]) # MPa 
        self.v21_fibre = v21_fiber
        self.E1_fibre  = E1_fiber 

    def get_v21(self):
        v21 = self.v12*self.E2/self.E1
        return v21
    
    u_get_v21 = uncertainties.wrap(get_v21)
    
class Lamina:
    def __init__(self, material, thickness, orientation, Q, z0, z1, strains, stresses, failure_counter, failure_mode): #, degradation_counter):
        self.material    =  material    # Material object
        self.thickness   =  thickness   # thickness layers
        self.orientation =  orientation # orientation in layup 
        self.Q           =  Q           # generally off axis, on axis for 0deg orientation   
        self.z0            =  z0
        self.z1            =  z1   
        self.strains       =  strains     # contains strains and curvatures
        self.stresses      =  stresses    # axial stresses due to axial loads and bending
        self.failure_counter = failure_counter
        self.failure_mode    = failure_mode
        #self.degradation_counter = degradation_counter

class Layup:
    def __init__(self, lay_up, symmetry, thickness):
        self.lay_up          = lay_up           # np.array as in [0,90] from [0,90]_s, 
        self.symmetry        = symmetry         # 0: no sym, 1: sym, -1: asym, 2: sym with last ply as midplane, -2: asym with last ply as midplane
        self.thickness       = thickness
        self.stack           = self.get_stack() # returns full stack of plies as in [0,90,90,0] from [0,90]_s
        self.ply_thicknesses = self.get_ply_thicknesses()
        self.total_thickness = sum(self.ply_thicknesses)

    def get_stack(self):
        if self.symmetry**2 == 4:
            stack = self.lay_up + [int(self.symmetry/2) * x for x in self.lay_up[::-1]][1:]
        else:
            stack = self.lay_up + self.symmetry**2*[self.symmetry * x for x in self.lay_up[::-1]]
        return stack
    
    def get_ply_thicknesses(self):
        if len(self.thickness) == 1:
            t = self.thickness[0] * np.ones_like(self.stack)
        else:
            t = self.thickness + [self.symmetry * x for x in self.thickness[::-1]]
        return t

    
class Laminate:
    def __init__(self, type, lay_up, laminae, A, B, D, ABD, abd, Ex, Ey, Gxy, vxy, vyx, maxstresses,maxstrains, minstresses,minstrains):
        self.type    = type           # string
        self.lay_up  = lay_up         # list
        self.laminae = laminae        # Lamina objects
        self.A, self.B, self.D        = A, B, D
        self.ABD     = ABD            # strains to force/moment intensities
        self.abd     = abd            # force/moment intensities to strains (useful one)
        self.Ex, self.Ey, self.Gxy, self.vxy, self.vyx     = Ex, Ey, Gxy, vxy, vyx             #EQUIVALENT PROPERTIES
        self.maxstresses = maxstresses      # 
        self.minstresses = minstresses 
        self.maxstrains  = maxstrains       # 
        self.minstrains  = minstrains       # 
        

def GenerateLaminate(lamina, lay_up, laminate): #, laminae):
    laminae = []
    # generating laminae orientations based upon specified lay-up 
    for i, (ply,t) in enumerate(zip(lay_up.stack, lay_up.ply_thicknesses)):
        lamina             = copy.deepcopy(lamina)                   # generates new basic ply
        lamina.orientation = lay_up.stack[i]                         # assigns orientation to new ply
        lamina.thickness   = lay_up.ply_thicknesses[i]               # assigns thickness to new ply
        laminae.append(lamina)                                       # collecting all plies i.e. laminae
    laminate.laminae = laminae

def Q_matrix(E1,E2,G12,v12,v21, theta):
    m = np.cos(np.deg2rad(theta))
    n = np.sin(np.deg2rad(theta))

    if v21 == None:
        v21 = v12*E2/E1

    Q = 1 - v12*v21
    Q11 = E1/Q
    Q22 = E2/Q
    Q12 = v12*E2/Q
    Q66 = G12

    Qxx = Q11*m**4 + 2*(Q12 + 2*Q66)*(m**2)*(n**2) + Q22*(n**4)
    Qxy = (Q11 + Q22 - 4*Q66)*(m**2)*(n**2) + Q12*(m**4 + n**4)
    Qyy = Q11*(n**4) + 2*(Q12 + 2*Q66)*(m**2)*(n**2) + Q22*(m**4)
    Qxs = (Q11 - Q12 - 2*Q66)*n*(m**3) + (Q12 - Q22 + 2*Q66)*(n**3)*m
    Qys = (Q11 - Q12 - 2*Q66)*m*(n**3) + (Q12 - Q22 + 2*Q66)*(m**3)*n
    Qss = (Q11 + Q22 - 2*Q12 - 2*Q66)*(n**2)*(m**2) + Q66*(n**4 + m**4)

    return Qxx, Qyy, Qxy, Qxs, Qys, Qss

def ABD_Matrix(lay_up,laminate):
    # calculating A, B, D using every lamina
    t = lay_up.total_thickness
    laminate.A =  np.zeros((3,3)) #,dtype='O')
    laminate.B =  np.zeros((3,3)) #,dtype='O')
    laminate.D =  np.zeros((3,3)) # ,dtype='O')

    for i, lamina in enumerate(laminate.laminae):
        E1, E2, G12, v12, v21 = lamina.material.E1, lamina.material.E2, lamina.material.G12, lamina.material.v12, lamina.material.v21
        Q = Q_matrix(E1, E2, G12, v12, v21, lamina.orientation)
        lamina.Q  = np.array([[Q[0], Q[2], Q[3]],
                            [Q[2], Q[1], Q[4]],
                            [Q[3], Q[4], Q[5]]])
        z0 = -t/2
        if i==0:
            lamina.z0 = z0
        # define zk-1
        if i != 0:
            z0 = laminate.laminae[i-1].z1
            laminate.laminae[i].z0 = z0

        # add thickness to get zk
        z1 = z0 + lamina.thickness
        lamina.z1 = z1
        #print(type(test[0][0].dtype))
        #print(type(laminate.A[0][0].dtype))
        laminate.A += lamina.Q*(z1 - z0)
        laminate.B += lamina.Q*(z1**2 - z0**2)/2 * (-1)    # SIGN WAS RESERVED SOMEHOW?
        laminate.D += lamina.Q*(z1**3 - z0**3)/3
    
        # equivalent properties
    Axx, Ayy, Axy, Ass = laminate.A[0][0], laminate.A[1][1], laminate.A[0][1], laminate.A[2][2]
    A = Axx*Ayy - Axy**2
    laminate.Ex, laminate.Ey   = A/(t*Ayy), A/(t*Axx)
    laminate.Gxy               = Ass/t
    laminate.vxy, laminate.vyx = Axy/Ayy, Axy/Axx

    # building full ABD matrix
    laminate.ABD = np.block([[laminate.A, laminate.B],[laminate.B, laminate.D]])
    #laminate.ABD = laminate.ABD.astype('float64')
    
def T(theta): # only to be used for in-plane analysis
    m = np.cos(np.deg2rad(theta))
    n = np.sin(np.deg2rad(theta))
    T_sigma = np.array([[m**2, n**2,       2*m*n],
                        [n**2, m**2,      -2*m*n],
                        [-m*n,  m*n, m**2 - n**2]])
    
    T_epsilon = np.array([[  m**2,  n**2,         m*n],
                          [  n**2,  m**2,        -m*n],
                          [-2*m*n, 2*m*n, m**2 - n**2]])

    return T_sigma, T_epsilon

