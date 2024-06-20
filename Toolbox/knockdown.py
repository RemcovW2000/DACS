import numpy as np
import scipy as sp



def CAI_laminate():

    def scf_5_2(l, v12, E11, E22, G12,a):
        '''
        
        '''
        num = 1 + (l + (1-l*v12**2*(E22/E11))*np.sqrt(2*(np.sqrt(E11/E22) -v12)+(E11/G12))) + ((E11/G12)-v12)*np.sqrt(E22/E11)
        den = 1 + l*(l+(1+np.sqrt(E22/E11)*np.sqrt(2*(np.sqrt(E11/E22) -v12)+(E11/G12)))) +((E11/G12)-2*l*v12)*np.sqrt(E22/E11)-(1-l)**2*v12**2*(E22/E11)
        k_inf = 1 - (1-l)*num/den
        SCF = (2+(1-(2*a/2)**3))/ ( 3 * ( 1 - (2*a/2))) * k_inf 
        return SCF
    
def CAI_strength(strength,SCF):
    s_CAI = strength/SCF
    return s_CAI