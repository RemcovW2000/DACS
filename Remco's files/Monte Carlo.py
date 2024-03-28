import numpy as np

class MonteCarlo:
    def __init__(self, plys, elasticproperties, failureproperties, statisticalproperties):
        self.plys = plys
        self.elasticproperties = elasticproperties
        self.failureproperties = failureproperties
        self.statisticalproperties = statisticalproperties

    # make a function that makes the standard layers with orientation
    def Makelaminate(self):
        # we make the standard laminas:
        laminas =
        for lamina in :
            s0 = Lamina(0.25, 0, elasticproperties, failureproperties)
            lamina =

# # Statistical properties:
# self.statisticalproperties = statisticalproperties
# self.E1std = statisticalproperties[0]
# self.E2std = statisticalproperties[1]
# self.v12std = statisticalproperties[2]
# self.G12std = statisticalproperties[3]
# self.Xtstd = statisticalproperties[4]
# self.Xcstd = statisticalproperties[5]
# self.Ytstd = statisticalproperties[6]
# self.Ycstd = statisticalproperties[7]
# self.Sstd = statisticalproperties[8]