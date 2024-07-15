from Toolbox.Lamina import Lamina
from Toolbox.Laminate import Laminate, LaminateBuilder

laminate = LaminateBuilder([0,90], True, True, 1, 'T700')
laminate.Loads = [-100, 0 ,0 ,0, 0, 0]
print(laminate.Ncrit())
print(laminate.vxy)
