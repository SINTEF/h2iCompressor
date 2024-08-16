"""
The following Python code takes a set of variables as inputs and calculates the geometry of an inducer, 
impeller and diffuser for a centrifugal compressor.

Reference: Lowe, Mitchell (2016-08-21). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244
Author: Petter Resell (SINTEF Energy Research, 2024)

"""

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# import main as m

# Plotting files:
from plotFlow import plotFlowConditions
from plotCompressor import plotCompressorParam
from findFlow import inducerFlowWRTminimumRelativeVelocity
from pressureTest import pressureOverUnderEstimate



def initialize():

    ### Variables-------------------------------------------------------------------------------------------
    # Inlet flow parameters
    global mdot, R, k, Cp, MolarMass, impellerDensity, impellerTensileStrength, etaStage0
    global T00, P00, T00i, Pr, Cm1i, N0
    global alpha1, B1, B2, AR, lambda2, lambda20
    global etad, CpD, CpDi
    global rhDivr1, r1Divr2
    global etaLowerLimit, etaUpperLimit, bladeVelUpperLimit, bladeVelLowerLimit, betamax, betamin, bladeMin, bladeMax, iterTol


    mdot = 0.4                               # Mass flow rate [kg/s]   , mdot = ( (300* 10**3 )/( 60**2 ) )/6  
    CpAir = 1005                            # Cp air [J/kg/K], 
    CpH2 = 14310                            # Engineering toolbox Cp hydrogen, [J/kg/K]                                       
    Cp = CpAir                              # Choose Air or Hydrogen                                      
    MolarMassAir = 0.02897                  # Molecular weight of air [kg/mol]. 
    MolarMassH2 = 2.01568 * 10**-3          # Molecular weight of H2 [kg/mol]. 
    MolarMass = MolarMassAir                # Molecular weight of H2 [kg/mol]. 
    k = 1.41                                # Basicly similar for air and hydrogen
    R_uni = 8.314                           # Universal gas constant [J/mol/K]
    R = R_uni / MolarMass                   # Specific gas constant [J/kg /K]

    """ ------ inlet conditions ------ """
    P00 = 1 * (10**2) * (10**3)                    # Inlet stagnation pressure [Pa]
    T00 = 293                                       # Inlet stagnation temperature [K]
    Pr = 4.8                                        # Pressure ratio [-]
    Cm1i = np.arange(10, 600.5, 0.5)                # Inlet Absolute meridional velocity [m/s]
    alpha1 = 0                                      # Absolute inlet velocity angle [deg]
    B1 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant
    B2 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant NEED GREATER BLOCKAGE FOR OUTLET
    AR = 2.5                                        # inlet/outlet area ratio of the diffuser [-]         
    T00i = np.full(len(Cm1i), 293)                  # Inlet stagnation temperature [K]

    """ ------ Ti-6Al-2Sn-4Zr-2Mo ------> titanium  properties https://www.azom.com/article.aspx?ArticleID=9298 """
    impellerDensity = 1540                     # [kg/m3]
    impellerTensileStrength = 900* 10**6       # Approximate ultimate tensile strength UTS [Pa]

    """ Impeller exit flow parameters """
    lambda2 = 2                 # Exit swirl parameter                                          , check out
    lambda20 = lambda2          # Exit swirl parameter, used for iteration
    etaStage = 0.6             # Isentropic Stage efficiency [-]           
    etaStage0 = etaStage                               

    # Ndes = 120000                     # Rotational speed [rpm] 
    # Ndes = 60000                      # Rotational speed [rpm]   
    N0 = 120000
                              # Rotational speed [rpm]     
    

    """ -------------- Diffuser Calculation -------------- """
    etad = 0.85                                 # Estimated diffuser efficiency
    CpDi = 1 - 1 / (AR ** 2)                    # Ideal pressure recovery coefficient
    CpD = etad * CpDi    

    """ ------------- Iteration control ----------"""
    etaLowerLimit = 0.35                        # Lowest efficiency allowed
    etaUpperLimit = 0.80                        # Highest efficiency allowed
    bladeVelUpperLimit = 1200                   # Highest blade velocity allowed
    bladeVelLowerLimit = 0                      # lowest blade velocity allowed 
    betamax = -45                               # "Maximum" beta iterated over
    betamin = +1                                 # "Minimum" beta iterated over  
    bladeMin = 1                                # Lowest bladenumber allowed
    bladeMax = 30                               # Highest bladenumber allowed
    iterTol = 0.05                            # loop tolerance condition
    

    """ ------- VARY THESE PARAMETERS ------- """

    rhDivr1=0.364                                 # Ratio commonly given
    r1Divr2 = 0.522                                # rt1 divided by rt2 ,making rt2 increase through array progression
    # r1Divr2 = np.arange(0.8, 0.3, -0.05)        # rt1 divided by rt2 ,making rt2 increase through array progression

    """ ------------------------------------- """

initialize()                            # Initializing values above


# file = open("./geometryV7.py")
# exec(file.read())                       # Running script containing main functionality
# file.close()

