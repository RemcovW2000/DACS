import numpy as np
from Laminate import Laminate
from Lamina import Lamina
from stringers import TStringer
from Panel import Panel
import copy


class Fuselage:
    def __init__(self, diameter, standardstringers = [], stringertypes = [], stringerangles = [],
                 stringers = [], standardpanels = [], paneltypes = [], panels = [], loads=None):
        # We want the stringers to have a location

        # Stringer types: A list of types, one for each stringer in the angles list
        self.stringertypes = stringertypes

        # We have to sort the stringer angles:
        sorted_angles = sorted(stringerangles)
        self.stringerangles = sorted_angles         # List of angles where stringers are
        self.standardstringers = standardstringers  # List of standard stringers
        self.stringers = stringers                  # List of stringer objects

        self.loads = loads          # The given loads

        self.diameter = diameter
        self.framedistance = 1             # Distance from frame to frame

        self.standardpanels = standardpanels
        self.paneltypes = paneltypes
        self.panels = panels

        # Upon initialisation we 'build' the fuselage
        self.GenerateStringers()
        self.GeneratePanels()


    def GenerateStringers(self):
        self.stringers = []

        # Generate stringers at angles on the fuselage:
        for count, i in enumerate(self.stringerangles):
            # Define the position:
            x = np.sin(np.deg2rad(i)) * self.diameter * -0.5
            y = np.cos(np.deg2rad(i)) * self.diameter * 0.5
            stringerposition = np.array([x, y])

            # Assign the stringer
            stringertype = self.stringertypes[count]

            # Deepcopy the stringer from the standard stringer list
            stringer = copy.deepcopy(self.standardstringers[stringertype])
            stringer.Position = stringerposition
            self.stringers.append(stringer)
        return

    def GeneratePanels(self):
        # We'll generate panels between the stringers, they will have a start and end coordinate
        self.panels = []

        for i in range(len(self.paneltypes)):
            # Get the start and end positions of the current panel
            start_position = self.stringers[i].Position
            end_position = self.stringers[i + 1].Position

            # Create a 2x2 array to represent the panel's start and end positions
            panel_positions = np.array([start_position, end_position])

            # Find the correct panel type:
            paneltype = self.paneltypes[i]
            print('panel type: ',paneltype)

            panel = copy.deepcopy(self.standardpanels[paneltype])
            panel.Positions = panel_positions
            panel.Length = self.framedistance
            print('Panel positions:')
            print(panel.Positions)

            print('--------------------------')
            self.panels.append(panel)
        return

    def EIEquivalence(self):

        return

    def CalulatePanelEA(self):
        for panel in self.panels:
            panel.CalculateEA()
        return

    def CalculateYbar(self):
        StringerEAysum = 0
        StringerEAsum = 0

        EApreviouspanel = 0
        # First calculate the sum of y*EA of all the stringers
        for count, i in enumerate(self.stringers):
            EAStringer = i.EAtotal
            # For the first or last stringer we must take only half the EA:
            if count == 0 or count == len(self.stringers):
                EAstringer = 0.5 * EAStringer

            EAnextpanel = self.panels[count].EA/2
            EANode = EAstringer + EApreviouspanel + EAnextpanel
            y = i.Position[1]
            StringerEAysum += i.EAtotal*y
            StringerEAsum += i.EAtotal
            EApreviouspanel = EAnextpanel

        PanelEAysum = 0
        PanelEAsum = 0

        self.ybar = (PanelEAysum + StringerEAysum)/(PanelEAsum + StringerEAsum)
        return self.ybar


# now we test the code:
E1 = 142e3     # From assignment
E2 = 11.2e3       # From assignment
G12 = 5e3      # From assignment
v12 = 0.3      # From assignment
elasticproperties = [E1, E2, G12, v12]

# properties needed for failure analysis
E11f = 324e3   # TBD
v21f = 0.2      # TBD
msf = 1.1       # TBD
R11t = 2200     # From assignment
R11c = 1800     # From assignment
yt = 70         # From assignment
yc = 300        # From assignment
S = 100         # From assignment
t = 0.135       # thickness

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]

V0 = Lamina(t, 0, elasticproperties, failureproperties)
V1 = Lamina(t, 0, elasticproperties, failureproperties)
V2 = Lamina(t, 0, elasticproperties, failureproperties)
V3 = Lamina(t, 0, elasticproperties, failureproperties)

H0 = Lamina(t, 0, elasticproperties, failureproperties)
H1 = Lamina(t, 0, elasticproperties, failureproperties)
H2 = Lamina(t, 0, elasticproperties, failureproperties)
H3 = Lamina(t, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
LaminasV = [V0, V1, V2, V3]
LaminasH = [H0, H1, H2, H3]

LaminateH = Laminate(LaminasH)
LaminateV = Laminate(LaminasV)

# We make the stringer objects:
TStringer1 = TStringer(LaminateH, LaminateV, 30, 30)

TStringer2 = copy.deepcopy(TStringer1)

TStringer3 = copy.deepcopy(TStringer1)

standardstringers = [TStringer1, TStringer2, TStringer3]
stringerangles = [0, 20, 60, 80, 100, 120, 140, 160, 180]
stringertypes = [0, 0, 0, 0, 0, 0, 0, 0, 0]

P0 = Lamina(t, 0, elasticproperties, failureproperties)
P1 = Lamina(t, 90, elasticproperties, failureproperties)
P2 = Lamina(t, 45, elasticproperties, failureproperties)
P3 = Lamina(t, -45, elasticproperties, failureproperties)
quasiisotropiclaminas = [P0, P1, P2, P3]
quasiisotropiclaminate = Laminate(quasiisotropiclaminas)

P4 = Lamina(t, 45, elasticproperties, failureproperties)
P5 = Lamina(t, 0, elasticproperties, failureproperties)
P6 = Lamina(t, 0, elasticproperties, failureproperties)
P7 = Lamina(t, -45, elasticproperties, failureproperties)
zerolaminas = [P4, P5, P6, P7]
zerolaminate = Laminate(zerolaminas)

#Define standard panels list:
standardpanels = [Panel(Laminate=quasiisotropiclaminate), Panel(Laminate=zerolaminate)]

paneltypes = [1, 1, 0, 0, 0, 0, 1, 1]

Fuselage = Fuselage(6, standardstringers, stringertypes, stringerangles,
                    standardpanels = standardpanels, paneltypes = paneltypes)
#print('ybar:', Fuselage.CalculateYbar())