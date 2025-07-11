import Data.MP

from Toolbox import structural_entity
from Toolbox.laminate import laminate_builder

angles = [0, 0, -10, 10, -20, 20, -30, 30, -40, 40, -50, 50, -60, 60, -70, 70, -80, 80, -90, 90]

laminate = laminate_builder(angles, True, True, 1, 'T700')
print('Vxy: ' , laminate.vxy)
print('Ex, Ey: ', laminate.Ex)

laminate2 = laminate_builder([0], True, True, 1, 'T700')
print(laminate2.vxy)
print('Ez = ', laminate2.Ey)
E3 = laminate2.Ey
G13 = laminate2.Gxy
v13 = (E3 / (2*G13)) -1
print('v13: ', v13)
print('Gxz: ', laminate2.Gxy)
print('Gxy: ', laminate.Gxy)
load = 300
laminate.Loads = [0, 0, -load*laminate.h, 0, 0, 0]
Tx = load/laminate.failure_analysis()
