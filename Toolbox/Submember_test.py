from Data.Airfoils import airfoilcoords
from scipy.interpolate import interp1d
from Toolbox.Laminate import LaminateBuilder
import matplotlib.pyplot as plt
from Toolbox.DamagedRegion import *
from Toolbox.Member import Member
from scipy.integrate import simps
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.linalg import lstsq
import math
from Toolbox.section import section
from Toolbox.helperclasses import boom, segment

mainlaminate = LaminateBuilder([0, 90, 45, -45], True, True, 1)

reinforcedlaminate = LaminateBuilder([0, 90, -45, 45, 0, 0, 0], True, True, 1)

member = Member(mainlaminate, reinforcedlaminate)
member.subpanel = reinforcedlaminate
member.subpanel_start = 0
member.subpanel_end = 20

print(member.Ex(30))