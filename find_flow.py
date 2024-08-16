import math
import numpy as np
import matplotlib.pyplot as plt

def inducerFlowWRTminimumRelativeVelocity(N, rt1arr, Ctheta1arr, Cm1arr):
    U1ti = 2 * math.pi * rt1arr * N / 60                       # Inlet blade tip speed [m/s]                    eqn (60) rapport
    W1ti = (Cm1arr ** 2 + (U1ti - Ctheta1arr) ** 2) ** 0.5        # Todo    : enig i denne linjen, ikke den over
    indice_min = np.argmin(W1ti)
    rt1 = rt1arr[indice_min]           # Replaced by prior loop 
    Ctheta1 = Ctheta1arr[indice_min]     # Get the value of Ctheta1 for minimal value of W1t
    C1 = C1i[indice_min]        
    T1 = T1i[indice_min]
    M1 = M1i[indice_min]
    P1 = P1i[indice_min]
    rho1 = rho1i[indice_min]
    A1 = A1i[indice_min]
    U1t = U1ti[indice_min]
    Cm1 = Cm1i[indice_min]
    beta1 = math.degrees(math.atan((U1t - Ctheta1) / Cm1))      # Inlet relative velocity angle [deg]
    W1t = W1ti[indice_min]
    omega = U1t/rt1                 # =2*np.pi*N/60
    return U1ti,W1ti,rt1,Ctheta1,C1,T1,M1,P1,rho1,A1,U1t,Cm1,beta1,W1t,omega