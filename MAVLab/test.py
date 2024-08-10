import numpy as np
from Toolbox.Member import Member
from Toolbox.Laminate import Laminate, LaminateBuilder
from Data.Panels import LaminateQI

Loads = [1000, 0, 200, 0, 0, 0]
member = Member(LaminateQI, Loads)
# member.R = 1000
print(member.BucklingAnalysis())