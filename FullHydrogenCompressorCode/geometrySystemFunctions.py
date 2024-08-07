# import geometryV8.py
import settingsOffDesign
import math


""" The function utilize the slip factor in order to find the impeller outlet velocities. It also finds the work input coefficient"""

def impellerOutletVelocities(slipFactor, beta2B, U2):
    workInputCoeff = slipFactor * settingsOffDesign.lambda2 / (settingsOffDesign.lambda2 - math.tan(beta2B))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 
    Ctheta2m = workInputCoeff * U2                                                # Absolute tangential exit velocity [m/s]         from work coefficient
    Cm2m = Ctheta2m / settingsOffDesign.lambda2                       # Absolute meridional exit velocity [m/s]         
    C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5 

    return workInputCoeff, Ctheta2m, Cm2m, C2


""" The diffuserFlow-function finds the relevant diffuser variables """
def diffuserFlow(P2, P02, rho2, C2):
    P3 = P2 + settingsOffDesign.CpD * (P02 - P2)                    # Diffuser exit static pressure [Pa]
    C3 = C2 / settingsOffDesign.AR                                  # Diffuser exit absolute velocity [m/s]
    P03 = P3 + 0.5 * rho2 * C3 ** 2                                 # Diffuser exit stagnation pressure [Pa]

    return P3, P03, C3


""" systemTotalPerformance determines the overall performance of the compressor. This is done through finding the """
def systemTotalPerformance(P03, T02, U2, T1, workInputCoeff):
    etaIterate = ( (P03 / settingsOffDesign.P00) ** ((settingsOffDesign.k - 1) / settingsOffDesign.k) - 1 ) / ( (T02 / settingsOffDesign.T00) - 1 )        # Iterative stage efficiency [-], isothermal diffuser assumed?
    PressureRatioEstimate = ((etaIterate * U2 ** 2 * workInputCoeff) / (settingsOffDesign.Cp * T1) + 1) ** (settingsOffDesign.k / (settingsOffDesign.k - 1))          # Estimate of the pressure ratio, equation is validated
    PressureTestOuterLoop = (PressureRatioEstimate-settingsOffDesign.Pr)/settingsOffDesign.Pr

    return etaIterate, PressureRatioEstimate, PressureTestOuterLoop 