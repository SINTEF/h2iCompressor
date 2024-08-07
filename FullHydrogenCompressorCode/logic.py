"""
The following Python code iterates through a range of rotational speeds to deteremine compressor off-design performance.
A geoetrical basis i found in geometry.py which are further built on in off-design_performance.py. Plotting of several variables of interest 
is done by utilizing plotCompessor.py, plotSystem.py, plotText.py, plotVelocities.py and pressureTest.py.
The script called settingsOffDesign is used for all parametrization of for instance fluid properties, inlet flow conditions, diffuser propertis etc. 

Author: Petter Resell (SINTEF Energy Research, 2024)

"""


import math
import numpy as np
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import geometryV8
import off_design_performanceV8 as odp
import settingsOffDesign
from plotFlow import plotFlowConditions
from plotSystem import  plotSystemVariables
from plotCompressor import plotCompressorParam
from plotVelocities import plotVelocities
from plotText import plotText

# Set plot parameters------------------------------------------------------------------------------
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})
plt.close('all')
# ------------------------------

showGeometryPlots = plt.show()


# ------------------- Plotting from geometry.py -------------------
plotCompressorParam(geometryV8.systemVar, geometryV8.Zcompressor, geometryV8.designParam, geometryV8.flowVar)
plotFlowConditions(geometryV8.systemVar, geometryV8.Zflow, geometryV8.designParam, geometryV8.flowVar)
plotVelocities(geometryV8.systemVar, geometryV8.Zvelocities, geometryV8.designParam, geometryV8.flowVar)
plotSystemVariables(geometryV8.systemVar, geometryV8.Zsystem, geometryV8.designParam, geometryV8.flowVar)
plotText(geometryV8.text)

print('r1: ' + str(round(geometryV8.r1, 3)) + 'm')
print('r2: ' + str(round(geometryV8.r2, 3)) + 'm' )

""" Running off-desing_performance.py for rotational speeds of 15000rpm to the critical speed with increments of 5000 """
Narr = np.arange(15000, geometryV8.Ncrit , 5000)      

# ------------------- Plotting from off_design_performance.py -------------------
for iN in range(0, len(Narr)):                                                                      # Running and plotting for all rpm N's
    Pro, P03o, T2oabs, mdoto, etao, U2o = odp.off_design_performance(Narr[iN])                      # Running everything
    Pro = np.array(Pro)                                                                             # Pressure ratio
    P03o = np.array(P03o)                                                                           # Outlet pressure
    T2oabs = np.array(T2oabs)                                                                       # Outlet temperature
    mdoto = np.array(mdoto)                                                                         # mass flow
    etao = np.array(etao)                                                                           # Efficiency
    correctedMass = mdoto*(P03o/settingsOffDesign.P00)/np.sqrt(T2oabs/settingsOffDesign.T00)        # Corrected mass flow parameter

    """ checking if all efficiencies are greater than 0 and less or equal to 1 """
    if all(etao > 0) and all(etao <=1):                                 
        # Plotting pressure ratio curve
        fig1 = plt.figure('Pressure ratio')
        plt.plot(mdoto, Pro, label= str(Narr[iN])  ) # marker = markers[iN],
        plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m}$', fontsize = 15)
        plt.legend()
        plt.plot([0, 50], [1, 1], 'r--')
        plt.xlim([0, 50])
        plt.title(r'Pressure ratio')
        plt.grid(True)

        # Compressor map 
        fig99 = plt.figure('Compressor map')
        plt.plot(correctedMass, Pro, label= str(Narr[iN]) )# marker = markers[iN]
        plt.ylabel( r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m} \frac{P_{03}}{P_{00}} \sqrt{\frac{T_{00}}{T_{03}} }$', fontsize = 15)
        plt.legend()
        plt.plot([0, 50], [1, 1], 'r--')
        plt.xlim([0, 50])
        plt.title(r'Compressor map')
        plt.grid(True)
        
        # Plotting efficiency curve 
        fig2 = plt.figure('Efficiency')
        plt.plot(correctedMass, etao, label= str(Narr[iN])  ) # marker = markers[iN],
        plt.grid()
        plt.ylabel(r'$\eta$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m} \frac{P_{03}}{P_{00}} \sqrt{\frac{T_{00}}{T_{03}} }$', fontsize = 15)
        plt.legend()
        plt.title(r'Efficiency ')
        plt.grid(True)
        
        fig4 = plt.figure('Exit Tip velocity')
        plt.plot(Narr[iN],U2o, 'ko', label= str(Narr[iN])  ) # marker = markers[iN],
        plt.grid()
        plt.xlabel(r'$N $  [rpm]', fontsize = 15)
        plt.ylabel(r'$U_2 $ [m/s]', fontsize = 15)
        plt.grid(True)


        debug = 0
    else:
        print('   -> RPM  = ' + str(Narr[iN]) + ' has efficiency outside range [0, 1] !')


print('Off-design_performance.py successfully run. \n')


plt.show()