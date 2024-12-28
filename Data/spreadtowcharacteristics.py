from Toolbox.lamina import Lamina
from Toolbox.laminate import Laminate, laminate_builder

laminate = laminate_builder([0,90], True, True, 1, 'T700')
laminate.Loads = [-100, 0 ,0 ,0, 0, 0]
print(laminate.Ncrit())
print(laminate.vxy)
