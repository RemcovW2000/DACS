import numpy as np



def CAI_laminate(Eave,E11,E22,v12,G12,R,w):
    
    '''
    Eave - average modulus of the damaged laminate 
    E11 -
    R - radius 
    w - width
    Make sure units are consistent
    '''

    l = Eave/E11
  
    num = 1 + (l + (1-l*v12**2*(E22/E11))*np.sqrt(2*(np.sqrt(E11/E22) -v12)+(E11/G12))) + ((E11/G12)-v12)*np.sqrt(E22/E11)
    den = 1 + l*(l+(1+np.sqrt(E22/E11)*np.sqrt(2*(np.sqrt(E11/E22) -v12)+(E11/G12)))) +((E11/G12)-2*l*v12)*np.sqrt(E22/E11)-(1-l)**2*v12**2*(E22/E11)
    k_inf = 1 - (1-l)*num/den
    SCF = (2+(1-(2*R/w)**3))/ ( 3 * ( 1 - (2*R/w))) * k_inf 
    return SCF

def CAI_strength(strength,SCF):
    return strength/SCF 