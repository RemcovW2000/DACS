from Data.Airfoils import airfoilcoords
from scipy.interpolate import interp1d
from Toolbox.laminate import laminate_builder
import matplotlib.pyplot as plt
from Toolbox.damaged_region import *
from Toolbox.member import Member
from scipy.integrate import simps
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.linalg import lstsq
import math
from Toolbox.section import Section
from Toolbox.helperclasses import boom, segment
from Data.Panels import Laminates
from Toolbox.airfoil import Airfoil
from Toolbox.structural_entity import StructuralEntity

sparlocations = [80, 200]  # Example spar positions
chordlength = 350  # Example chord length
sparmembers = [Member(Laminates['steel_laminate']) for _ in sparlocations]
topmembers = [Member(Laminates['steel_laminate']) for _ in range(len(sparlocations) + 1)]
botmembers = [Member(Laminates['steel_laminate']) for _ in range(len(sparlocations) + 1)]


airfoil = Airfoil(
airfoil_name='NACA2410',
thickness= 0.1,
chord_length=chordlength,
spar_locations=sparlocations,
top_members=topmembers,
bot_members=botmembers,
spar_members=sparmembers,
)
airfoil.segment_length = 0.01
airfoil.plot=True

mx = 100
my = 0
sx = 0
sy = 100
moments = [mx, my]
shears = [sx, sy]
center = [100, 0]  # place it at quarter chord
airfoil.solve_stresses_CSA(moments, shears, center)

airfoil.plot_shear_flow().show()