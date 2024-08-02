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
# import logic
# import main as m

# Plotting files:
from plotFlow import plotFlowConditions
from plotCompressor import plotCompressorParam
from findFlow import inducerFlowWRTminimumRelativeVelocity
from pressureTest import pressureOverUnderEstimate


# showGeometryPlots = plt.show()
showGeometryPlots = 0

# def initialize():

### Variables-------------------------------------------------------------------------------------------
# Inlet flow parameters
global mdot, R, k, Cp, MolarMass, impellerDensity, impellerTensileStrength, etaStage0
global T00, P00, T00i, Pr, Cm1i, N0, rho0
global alpha1, B1, B2, AR, lambda2, lambda20
global etad, CpD, CpDi
global rhDivr1, r1Divr2
global etaLowerLimit, etaUpperLimit, bladeVelUpperLimit, bladeVelLowerLimit, beta2Bmax, beta2Bmin, bladeMin, bladeMax, iterTol


mdot = 25                               # Mass flow rate [kg/s]   , mdot = ( (300* 10**3 )/( 60**2 ) )/6  
CpAir = 1006                            # Cp air [J/kg/K], 
CpH2 = 14310                            # Engineering toolbox Cp hydrogen, [J/kg/K]                                       
Cp = CpH2                               # Choose Air or Hydrogen                                      
MolarMassAir = 0.02897                  # Molecular weight of air [kg/mol]. 
MolarMassH2 = 2.01568 * 10**-3          # Molecular weight of H2 [kg/mol]. 
MolarMass = MolarMassH2                 # Molecular weight of H2 [kg/mol]. 
k = 1.41                                # Basicly similar for air and hydrogen
R_uni = 8.314                           # Universal gas constant [J/mol/K]
R = R_uni / MolarMass                   # Specific gas constant [J/kg /K]

""" ------ inlet conditions ------ """
P00 = 24 * (10**2) * (10**3)                    # Inlet stagnation pressure [Pa]
T00 = 293                                       # Inlet stagnation temperature [K]
Pr = 1.24                                        # Pressure ratio [-]
Cm1i = np.arange(10, 600.5, 1)                # Inlet Absolute meridional velocity [m/s]
alpha1 = 0                                      # Absolute inlet velocity angle [deg]
B1 = 0.04                                       # Boundary layer blockage [-],        substance viscocity dependant
AR = 2.5                                        # inlet/outlet area ratio of the diffuser [-]         
T00i = np.full(len(Cm1i), 293)                  # Inlet stagnation temperature [K]
rho0 = P00/(R*T00)

""" ------ Ti-6Al-2Sn-4Zr-2Mo ------> titanium  properties https://www.azom.com/article.aspx?ArticleID=9298 """
impellerDensity = 1540                     # [kg/m3]
impellerTensileStrength = 900* 10**6       # Approximate ultimate tensile strength UTS [Pa]
# impellerDensity = 2710                     # [kg/m3]
# impellerTensileStrength = 690* 10**6       # Approximate ultimate tensile strength UTS [Pa]

""" Impeller exit flow parameters """
lambda2 = 2                 # Exit swirl parameter                                          , check out
lambda20 = lambda2          # Exit swirl parameter, used for iteration
etaStage = 0.6             # Isentropic Stage efficiency [-]           
etaStage0 = etaStage                               

# Ndes = 120000                     # Rotational speed [rpm] 
# Ndes = 60000                      # Rotational speed [rpm]   
N0 = 100000
# N0 = 47000

# Point of interest for off-design evaluation
bladeAngle = np.deg2rad(-40)
bladeNumber = 17




""" -------------- Diffuser Calculation -------------- """
etad = 0.85                                 # Estimated diffuser efficiency
CpDi = 1 - 1 / (AR ** 2)                    # Ideal pressure recovery coefficient
CpD = etad * CpDi    

""" ------------- Iteration control ----------"""
etaLowerLimit = 0.25                        # Lowest efficiency allowed
etaUpperLimit = 0.85                        # Highest efficiency allowed
bladeVelUpperLimit = 1200                   # Highest blade velocity allowed
bladeVelLowerLimit = 0                      # lowest blade velocity allowed 
beta2Bmax = -65                               # "Maximum" beta iterated over
beta2Bmin = 0                                 # "Minimum" beta iterated over  
bladeMin = 2                                # Lowest bladenumber allowed
bladeMax = 50                               # Highest bladenumber allowed
iterTol = 0.01                           # loop tolerance condition


""" ------- VARY THESE PARAMETERS ------- """

rhDivr1=0.3                                 # Ratio commonly given
r1Divr2 = 0.75                                # rt1 divided by rt2 ,making rt2 increase through array progression


    # r1Dr2 = np.arange(0.8, 0.3, -0.05)        # rt1 divided by rt2 ,making rt2 increase through array progression

""" ------------------------------------- """

# initialize()                            # Initializing values above


# file = open("./geometryV7.py")
# exec(file.read())                       # Running script containing main functionality
# file.close()

