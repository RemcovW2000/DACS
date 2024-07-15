from Toolbox.Laminate import LaminateBuilder, Laminate
from Toolbox.Sandwich import Sandwich
from Data.MP import corematerials
from Toolbox.Core import Core

LaminateQI = LaminateBuilder([45, -45, 0, 90], True, True, 1, 'T700')
FacesheetInner = LaminateBuilder([0], False, True, 1, 'T700')
FacesheetOuter = LaminateBuilder([0, 45], False, True, 1, 'T700')
core = Core(2, corematerials['ROHACELL31A'])
PanelWingRoot = Sandwich(FacesheetInner, FacesheetOuter, core, None, None)
