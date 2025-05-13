def impellerOutletVelocities(slipFactor, beta2B, U2, lambda2):
    """ Utilizes the slip factor in order to find the impeller outlet velocities. It also finds the work input coefficient """
    # Import 
    import math
    workInputCoeff = slipFactor * lambda2 / (lambda2 - math.tan(beta2B))      # Work input coefficient [-]  
    print('d', workInputCoeff)
    Ctheta2m = workInputCoeff * U2                  # Absolute tangential exit velocity [m/s]     from work coefficient
    Cm2m = Ctheta2m / lambda2                       # Absolute meridional exit velocity [m/s]         
    C2 = (Ctheta2m ** 2 + Cm2m ** 2) ** 0.5 

    return workInputCoeff, Ctheta2m, Cm2m, C2


def diffuserFlow(P2, P02, rho2, C2, CpD, AR):
    """ Finds the relevant diffuser variables """
    P3 = P2 + CpD * (P02 - P2)                    # Diffuser exit static pressure [Pa]
    C3 = C2 / AR                                  # Diffuser exit absolute velocity [m/s]
    P03 = P3 + 0.5 * rho2 * C3 ** 2               # Diffuser exit stagnation pressure [Pa]

    return P3, P03, C3


def systemTotalPerformance(P03, T02, U2, T1, workInputCoeff, P00, T00, k, Cp, Pr):
    """ Determines the overall performance of the compressor. This is done through finding the """
    etaIterate = ((P03 / P00) ** ((k - 1) / k) - 1) / ((T02 / T00) - 1)        # Iterative stage efficiency [-], isothermal diffuser assumed?
    PressureRatioEstimate = ((etaIterate * U2 ** 2 * workInputCoeff) / (Cp * T1) + 1) ** (k / (k - 1))          # Estimate of the pressure ratio, equation is validated
    PressureTestOuterLoop = (PressureRatioEstimate - Pr) / Pr

    return etaIterate, PressureRatioEstimate, PressureTestOuterLoop