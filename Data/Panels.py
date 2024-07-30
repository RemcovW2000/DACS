from Toolbox.Laminate import LaminateBuilder, Laminate
from Toolbox.Sandwich import Sandwich
from Data.MP import corematerials
from Toolbox.Core import Core

'''
File containing dictionaries containing LAMINATE - or SANDWICH objects for ease of reference.
'''

LaminateQI = LaminateBuilder([45, -45, 0, 90], True, True, 1, 'T700')
FacesheetInner = LaminateBuilder([0], False, True, 1, 'T700')
FacesheetOuter = LaminateBuilder([0, 45], False, True, 1, 'T700')
core = Core(2, corematerials['ROHACELL31A'])
PanelWingRoot = Sandwich(FacesheetInner, FacesheetOuter, core, None, None)

Laminates = {
    'CFQIChristos' : LaminateBuilder([45, -45, 0, 90], True, True, 1, 'Christos'),
    'CFUDChristos' : LaminateBuilder([0, 0, 0, 0], True, True, 1, 'Christos'),
    'CFQIT700' : LaminateBuilder([45, -45, 0, 90], True, True, 1, 'T700'),
    'CFUDT700' : LaminateBuilder([0, 0, 0, 0], True, True, 1, 'T700'),
    'CFQIezcomposites_spreadtow' : LaminateBuilder([45, -45, 0, 90], True, True, 1, 'ezcomposites_spreadtow')
}

Sandwiches = {
    'PanelWingRoot' : Sandwich(LaminateBuilder([0], False, True, 1, 'T700'),
                               LaminateBuilder([0, 45], False, True, 1, 'T700'),
                               Core(2, corematerials['ROHACELL31A'])
                               ),
    # 'PanelWingTip' : Sandwich(LaminateBuilder([0], False, False, 1, 'S2Glass80gsm'),
    #                           LaminateBuilder([0, 45], False, False, 1, 'S2Glass80gsm'),
    #                           Core(2, corematerials['ROHACELL31A'])
    #                           )
}