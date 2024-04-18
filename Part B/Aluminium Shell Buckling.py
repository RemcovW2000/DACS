import numpy as np
import matplotlib.pyplot as plt



Mcr   = 15e9
L     = 1e3 
v     = 0.29
r     = 3e3
E     = 69e3
gamma = 0.3063970424 # t = 1.326

beta = 2*L/(np.pi*r)

def KnockDownFactor(t):
    gamma = 1 - 0.731*(1 - np.exp(-np.sqrt(r/t)/16))
    return gamma

def ThicknessForBuckling(gamma, L):
    A = 12 * Mcr * (L**2) * (1 - v) / (np.pi**3 * r**2 * E)
    B = 12 * (gamma**2) * (L**4) * (1 - v**2) / (np.pi**4 * (r**2) * (1 + (2*L/(np.pi*r))**2)**2)
    C = (1 + (2*L/(np.pi*r))**2)**2

    t = np.roots([C, 0, B, -A])
    return t[-1].real

def plotting(L):
    t_list     = [0, 1.326]
    gamma_list = [0, gamma]

    while abs(t_list[-2] - t_list[-1]) > 0.00001:
        gamma_0 = gamma_list[-1]
        t_0     = ThicknessForBuckling(gamma_0, L)
        gamma_1 = KnockDownFactor(t_0)
        t_list.append(t_0)
        gamma_list.append(gamma_1)
        #print('t', t_0)

    return t_list[-1]

frame_spacings = np.arange(0.5e3, 5e3, 10)
skin_thicknesses = []

print(plotting(2e3))

#for frame_spacing in frame_spacings:
#    skin_thicknesses.append(plotting(L))


print(skin_thicknesses)
print(frame_spacings)

plt.plot(frame_spacings, skin_thicknesses)
plt.show()
