E1 = 142000     # From assignment
E2 = 11200       # From assignment
G12 = 5000      # From assignment
v12 = 0.3      # From assignment

elasticproperties = [E1, E2, G12, v12]

# properties needed for failure analysis
E11f = 500000   # TBD
v21f = 0.1      # TBD
msf = 1.1       # TBD
R11t = 2200     # From assignment
R11c = 1800     # From assignment
yt = 70         # From assignment
yc = 300        # From assignment
S = 100         # From assignment

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]

t = 0.135       # mm ply thickness
rho = 0.00161   # kg/mm3