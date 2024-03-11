import numpy as np
import copy

np.set_printoptions(linewidth=300, precision = 3)

class Material:
    def __init__(self, material, E1, E2, G12, v12, v21):
        self.material = material # string
        self.E1       = E1       # GPa
        self.E2       = E2       # GPa
        self.G12      = G12      # GPa
        self.v12      = v12      # -
        self.v21      = v21      # - 

class Lamina:
    def __init__(self, material, thickness, orientation, Q, z0, z1, epsilon, stresses):
        self.material    =  material    # Material object
        self.thickness   =  thickness   # thickness layers
        self.orientation =  orientation # orientation in layup 
        self.Q           =  Q           # generally off axis, on axis for 0deg orientation   
        self.z0          =  z0
        self.z1          =  z1   
        self.epsilon     =  epsilon     # contains strains and curvatures
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
    def __init__(self, type, lay_up, laminae, A, B, D, ABD, abd):
        self.type    = type           # string
        self.lay_up  = lay_up         # np.array
        self.laminae = laminae        # Lamina objects
        self.A, self.B, self.D        = A, B, D
        self.ABD     = ABD            # strains to force/moment intensities
        self.abd     = abd            # force/moment intensities to strains (useful one)

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
lay_up    = Layup([0,30,67], 0,  [0.2]) 
print('lay_up:',lay_up.stack)
print('total thickness:', lay_up.total_thickness)
material  = Material('UD CFRP', 140, 10, 5, 0.3, None) # GPa, [-]

# generate lamina
lamina   = Lamina(material, 0.125, None, None, None, None, None, None)
laminae  = []
laminate = Laminate('cross-ply', lay_up.stack, laminae, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)))

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
    
    # add thickness to get zk
    z1 = z0 + lamina.thickness
    lamina.z1 = z1

    print(z0, z1,  z1 - z0)
    print('###########################')
    laminate.A += lamina.Q*(z1 - z0)
    laminate.B += lamina.Q*(z1**2 - z0**2)/2 * (-1)    # SIGN WAS RESERVED SOMEHOW?
    laminate.D += lamina.Q*(z1**3 - z0**3)/3

# building full ABD matrix
laminate.ABD[0:3,0:3], laminate.ABD[0:3,3:6], laminate.ABD[3:6,0:3], laminate.ABD[3:6,3:6] = laminate.A, laminate.B, laminate.B, laminate.D

laminate.abd = np.linalg.inv(laminate.ABD)

ABD = laminate.ABD
abd = laminate.abd
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
F[0] = 100 # N/mm -> results in stresses in MPa (N/mm^2)

# calculate laminate strains
e = abd @ F # midplane strains full laminate
for i, lamina in enumerate(laminate.laminae):
    lamina.epsilon = T(lamina.orientation)[1] @ e

    if (i + 1) > n_plies/2:
        zk = ((i + 1) - n_plies/2)*lamina.thickness

    if (i + 1) <= n_plies/2:
        zk = -(n_plies/2 - i - 1)*lamina.thickness
    
    lamina.stresses = lamina.Q @ e[:3] + zk * e[3:]
