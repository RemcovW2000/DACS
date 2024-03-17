import numpy as np
import copy

np.set_printoptions(linewidth=300, precision = 1)

class Material:
    def __init__(self, material, E1, E2, G12, v12, v21, strengths, v21_fiber, E1_fiber):
        self.material = material             # string
        self.E1        = E1                  # GPa
        self.E2        = E2                  # GPa
        self.G12       = G12                 # GPa
        self.v12       = v12                 # -
        self.v21       = self.get_v21()      # - 
        self.strengths = strengths # np.array([Xt, Xc, Yt, Yc, S]) # GPa 
        self.v21_fibre = v21_fiber
        self.E1_fibre  = E1_fiber 

    def get_v21(self):
        v21 = self.v12*self.E2/self.E1
        return v21
    
class Lamina:
    def __init__(self, material, thickness, orientation, Q, z0, z1, strains, stresses):
        self.material    =  material    # Material object
        self.thickness   =  thickness   # thickness layers
        self.orientation =  orientation # orientation in layup 
        self.Q           =  Q           # generally off axis, on axis for 0deg orientation   
        self.z0          =  z0
        self.z1          =  z1   
        self.strains     =  strains     # contains strains and curvatures
        self.stresses    =  stresses    # axial stresses due to axial loads and bending

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
    def __init__(self, type, lay_up, laminae, A, B, D, ABD, abd, Ex, Ey, Gxy, vxy, vyx):
        self.type    = type           # string
        self.lay_up  = lay_up         # np.array
        self.laminae = laminae        # Lamina objects
        self.A, self.B, self.D        = A, B, D
        self.ABD     = ABD            # strains to force/moment intensities
        self.abd     = abd            # force/moment intensities to strains (useful one)
        self.Ex, self.Ey, self.Gxy, self.vxy, self.vyx     = Ex, Ey, Gxy, vxy, vyx             #EQUIVALENT PROPERTIES
 


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

# specify lay-up(including symmetry and thicknesses), materials
lay_up    = Layup([0,45,-45,90], 1,  [0.125]) 
print('lay_up:',lay_up.stack)
#print('total thickness:', lay_up.total_thickness)
material  = Material('UD CFRP', 140, 10, 5, 0.3, None, [1.5, 1.2, 0.05, 0.25, 0.07], 0.2, 230) # GPa, [-]

# generate lamina
lamina   = Lamina(material, 0.125, None, None, None, None, None, None)
laminae  = []
laminate = Laminate('cross-ply', lay_up.stack, laminae, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0)

# generating laminae orientations based upon specified lay-up 
for i, (ply,t) in enumerate(zip(lay_up.stack, lay_up.ply_thicknesses)):
    lamina = copy.deepcopy(lamina)                   # generates new basic ply
    lamina.orientation = lay_up.stack[i]             # assigns orientation to new ply
    lamina.thickness   = lay_up.ply_thicknesses[i]   # assigns thickness to new ply
    laminae.append(lamina)                           # collecting all plies i.e. laminae

# calculating A, B, D using every lamina
n_plies = len(lay_up.stack)
t = lay_up.total_thickness

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

    #print(z0, z1,  z1 - z0)
    #print('###########################')
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

laminate.abd = np.linalg.inv(laminate.ABD)

ABD = laminate.ABD
abd = laminate.abd
ABD[np.abs(ABD) < np.finfo(np.float32).eps] = 0
print('ABD:')
print(ABD)

# result of epsilon = abd @ F is strain for full laminate in direction 
# of loading. To ensure continuity, perfect bonding is assumed, therefore
# each ply undergoes equal strain in direction of loading.
# This strain serves as input to calculate strain for each ply in the
# respective fibre axis system, through 'matrix of directionality' i.e.
# a transformation to a different axis system oriented at the ply orientation
# angle around z (through-the-thickness)

def T(theta):
    m = np.cos(np.deg2rad(theta))
    n = np.sin(np.deg2rad(theta))
    T_sigma = np.array([[m**2, n**2, 0, 0,  0,       2*m*n],
                        [n**2, m**2, 0, 0,  0,      -2*m*n],
                        [   0,    0, 1, 0,  0,           0],
                        [   0,    0, 0, m, -n,           0],
                        [   0,    0, 0, n,  m,           0],
                        [-m*n,  m*n, 0, 0,  0, m**2 - n**2]])
    
    T_epsilon = np.array([[  m**2,  n**2, 0, 0,  0,         m*n],
                          [  n**2,  m**2, 0, 0,  0,        -m*n],
                          [     0,     0, 1, 0,  0,           0],
                          [     0,     0, 0, m, -n,           0],
                          [     0,     0, 0, n,  m,           0],
                          [-2*m*n, 2*m*n, 0, 0,  0, m**2 - n**2]])

    return T_sigma, T_epsilon

# generate force vector F, containing force and moment intensities
F = np.zeros((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)

# apply loadcase
F[0] = 1 # N/mm -> results in stresses in MPa (N/mm^2)

# calculate laminate strains
e = abd @ F # midplane strains and curvature full laminate
for i, lamina in enumerate(laminate.laminae):
    epsilon = T(lamina.orientation)[1] @ e # midplane strains [-] and curvatures [1/mm] for the respective lamina orientation
    e_min = epsilon[:3] + min(lamina.z0, lamina.z1) * epsilon[3:]
    e_max = epsilon[:3] + max(lamina.z0, lamina.z1) * epsilon[3:]
    lamina.strains = [e_min, e_max]
    lamina.stresses = lamina.Q @ lamina.strains

# FAILURE CRITERIA ACCORDING TO PUCK
def FibreFailure(stress, material): # stress = np.array with s1,s2,s6 ; mat_properties with [Xt,Xc,Yt,Yc,S]
    print('Fibre Failure')
    Xt, Xc, Yt, Yc, S = material.strengths
    s1, s2, s6        = stress[1]
    if s1 < 0:
        R = -Xc
    else:
        R = Xt
    # m_cfrp = 1.1, m_gfrp = 1.3
    f = 1/(R)*(s1 - (material.v21 - material.v21_fibre*1.1*material.E1/material.E1_fibre)*s2)
    return f

def MatrixFailure(stress, material):
    Xt, Xc, Yt, Yc, S = material.strengths
    s1, s2, s6        = stress[1]
    if s2 > 0: # mode A
        print('mode A')
        f = np.sqrt((s6/S)**2 + (1 - 0.3*Yt/S)**2*(s2/Yt)**2) + 0.3*s2/S
    
    s23A = S/(2*0.2)*(np.sqrt(1+2*0.2*Yc/S) - 1)
    p12_minus = 0.25
    p23_minus  = p12_minus*s23A/S # p12(-)????
    s12c = S*np.sqrt(1 + 2*p23_minus)

    if abs(s2/s6) <= s23A/abs(s12c) and abs(s2/s6) >= 0:
        print('mode B')
        f = (np.sqrt(s6**2 + (p12_minus*s2)**2) + p12_minus*s2)/S

    else:
        print('mode C')
        f = ((s6/(2*(1 + p23_minus*S)))**2 + (s2/Yc)**2)*(Yc/-s2)
    return f

for i,lamina in enumerate(laminate.laminae):
    print('LAMINA', i + 1)
    f_f = FibreFailure(lamina.stresses, lamina.material)
    print(f_f)
    f_m = MatrixFailure(lamina.stresses, lamina.material)
    print(f_m)
    print('#####################')