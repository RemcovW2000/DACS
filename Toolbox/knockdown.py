import numpy as np



def CAI_laminate(Eave,Elaminate,E11,E22,v12,G12,R,w):
    
    '''
    Eave - average modulus of the damaged laminate 
    E11 -
    R - radius
    w - width
    Make sure units are consistent
    '''

    l = Elaminate/Eave
  
    num = 1 + (l + (1-l*v12**2*(E22/E11))*np.sqrt(2*(np.sqrt(E11/E22) -v12)+(E11/G12))) + ((E11/G12)-v12)*np.sqrt(E22/E11)
    den = 1 + l*(l+(1+np.sqrt(E22/E11)*np.sqrt(2*(np.sqrt(E11/E22) -v12)+(E11/G12)))) +((E11/G12)-2*l*v12)*np.sqrt(E22/E11)-(1-l)**2*v12**2*(E22/E11)
    SCF1 = 1 - (1-l)*num/den
    SCF2 = (2+(1-(2*R/w)**3))/ ( 3 * ( 1 - (2*R/w))) * SCF1
    return max(SCF1,SCF2)

def CAI_strength(strength,SCF):
    return strength/SCF 
 