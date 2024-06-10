E1 = 142000     # From assignment
E2 = 11200       # From assignment
G12 = 5000      # From assignment
v12 = 0.3      # From assignment

elasticproperties = [E1, E2, G12, v12]

# properties needed for failure analysis
E11f = 6000   # TBD
v21f = 0.1      # TBD
msf = 1.1       # TBD
R11t = 2200     # From assignment
R11c = 1800     # From assignment
yt = 70         # From assignment
yc = 300        # From assignment
S = 100       # From assignment

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]

t = 0.135       # mm ply thickness
rho = 0.00161   # 1g/mm3

E1 = 142000     # From assignment
E2 = 11200       # From assignment
G12 = 5000      # From assignment
v12 = 0.3      # From assignment



HRH128props = {
    'Ez' : 538,
    'Sxz' : 3.31,
    'Gxz' : 110,
    'Syz' : 1.79,
    'Gyz' : 66,
    'Xc' : 12.62,
    'rho' : 128
}

HRH144props = {
    'Ez': 621,
    'Sxz': 3.55,
    'Gxz': 121,
    'Syz': 2.07,
    'Gyz': 76,
    'Xc': 14.48,
    'rho': 144
}

PoplarProps = {
    'Ez': 621,
    'Sxz': 3.55,
    'Gxz': 121,
    'Syz': 2.07,
    'Gyz': 76,
    'Xc': 14.48,
    'rho': 401
}