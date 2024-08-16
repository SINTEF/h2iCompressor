# import geometry.py
import settings
import math


""" The function utilize the slip factor in order to find the impeller outlet velocities. It also finds the work input coefficient"""

def impellerOutletVelocities(slipFactor, beta2B, U2):
    workInputCoeff = slipFactor * settings.lambda2 / (settings.lambda2 - math.tan(beta2B))                # Work input coefficient [-]    MSG: Can formula be adjusted? tan(-x=-tan(x)) 
    Ctheta2m = workInputCoeff * U2                                                # Absolute tangential exit velocity [m/s]         from work coefficient
    Cm2m = Ctheta2m / settings.lambda2                       # Absolute meridional exit velocity [m/s]         
    C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5 

    return workInputCoeff, Ctheta2m, Cm2m, C2


""" The diffuserFlow-function finds the relevant diffuser variables """
def diffuserFlow(P2, P02, rho2, C2):
    P3 = P2 + settings.CpD * (P02 - P2)                    # Diffuser exit static pressure [Pa]
    C3 = C2 / settings.AR                                  # Diffuser exit absolute velocity [m/s]
    P03 = P3 + 0.5 * rho2 * C3 ** 2                                 # Diffuser exit stagnation pressure [Pa]

    return P3, P03, C3


""" systemTotalPerformance determines the overall performance of the compressor. This is done through finding the """
def systemTotalPerformance(P03, T02, U2, T1, workInputCoeff):
    etaIterate = ( (P03 / settings.P00) ** ((settings.k - 1) / settings.k) - 1 ) / ( (T02 / settings.T00) - 1 )        # Iterative stage efficiency [-], isothermal diffuser assumed?
    PressureRatioEstimate = ((etaIterate * U2 ** 2 * workInputCoeff) / (settings.Cp * T1) + 1) ** (settings.k / (settings.k - 1))          # Estimate of the pressure ratio, equation is validated
    PressureTestOuterLoop = (PressureRatioEstimate-settings.Pr)/settings.Pr

    return etaIterate, PressureRatioEstimate, PressureTestOuterLoop 