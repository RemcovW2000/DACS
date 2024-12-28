from Toolbox.laminate import laminate_builder, Laminate
from Toolbox.sandwich import Sandwich
from Data.MP import corematerials
from Toolbox.core import Core

'''
File containing dictionaries containing LAMINATE - or SANDWICH objects for ease of reference.
'''

LaminateQI = laminate_builder([45, -45, 0, 90], True, True, 1, 'T700')
FacesheetInner = laminate_builder([0], False, True, 1, 'T700')
FacesheetOuter = laminate_builder([0, 45], False, True, 1, 'T700')
core = Core(2, corematerials['ROHACELL31A'])
PanelWingRoot = Sandwich(FacesheetInner, FacesheetOuter, core, None, None)

Laminates = {
    'CFQIChristos' : laminate_builder([45, -45, 0, 90], True, True, 1, 'Christos'),
    'CFUDChristos' : laminate_builder([0, 0, 0, 0], True, True, 1, 'Christos'),
    'CFQIT700' : laminate_builder([45, -45, 0, 90], True, True, 1, 'T700'),
    'CFUDT700' : laminate_builder([0, 0, 0, 0], True, True, 1, 'T700'),
    'CFQIezcomposites_spreadtow' : laminate_builder([45, -45, 0, 90], True, True, 1, 'ezcomposites_spreadtow'),
    'Reinforcement': laminate_builder([0, 0, 0, 0, 45, -45, 0, 0, 0], True, True, 1, 'T700')
}

Sandwiches = {
    'PanelWingRoot' : Sandwich(laminate_builder([0], False, True, 1, 'T700'),
                               laminate_builder([0, 45], False, True, 1, 'T700'),
                               Core(2, corematerials['ROHACELL31A'])
                               ),
    'SparPanels' : Sandwich(laminate_builder([45], False, False, 1, 'ezcomposites_spreadtow'),
                              laminate_builder([45], False, False, 1, 'ezcomposites_spreadtow'),
                              Core(2, corematerials['ROHACELL31A'])
                              )
}