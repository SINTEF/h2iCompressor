"""
The following Python code iterates through a range of rotational speeds to deteremine compressor off-design performance.
A geoetrical basis i found in geometry.py which are further built on in off-design_performance.py. Plotting of several variables of interest 
is done by utilizing plot_compessor.py, plot_system.py, plot_text.py, plot_velocities.py and pressure_test.py.
The script called settings is used for all parametrization of for instance fluid properties, inlet flow conditions, diffuser propertis etc. 

Authors: Petter Resell, Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


# Import 
import math
import numpy as np
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import geometry
import off_design_performance as odp
import settings
from plot_scripts.plot_flow import plotFlowConditions
from plot_scripts.plot_system import  plotSystemVariables
from plot_scripts.plot_compressor import plotCompressorParam
from plot_scripts.plot_velocities import plotVelocities
from plot_scripts.plot_text import plotText


# Set plot parameters------------------------------------------------------------------------------
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})
plt.close('all')
# ------------------------------

showGeometryPlots = plt.show()


# ------------------- Plotting from geometry.py -------------------
plotCompressorParam(geometry.systemVar, geometry.Zcompressor, geometry.designParam, geometry.flowVar)
# plotFlowConditions(geometry.systemVar, geometry.Zflow, geometry.designParam, geometry.flowVar)
plotVelocities(geometry.systemVar, geometry.Zvelocities, geometry.designParam, geometry.flowVar)
plotSystemVariables(geometry.systemVar, geometry.Zsystem, geometry.designParam, geometry.flowVar)
plotText(geometry.text)



""" Running off-desing_performance.py for rotational speeds of 15000rpm to the critical speed with increments of 5000 """
Narr = np.arange(15000, geometry.Ncrit , 5000)      

# ------------------- Plotting from off_design_performance.py -------------------
for iN in range(0, len(Narr)):                                                                      # Running and plotting for all rpm N's
    Pro, P03o, T2oabs, mdoto, etao, U2o, M2o = odp.off_design_performance(Narr[iN])                      # Running everything
    Pro = np.array(Pro)                                                                             # Pressure ratio
    P03o = np.array(P03o)                                                                           # Outlet pressure
    T2oabs = np.array(T2oabs)                                                                       # Outlet temperature
    mdoto = np.array(mdoto)                                                                         # mass flow
    etao = np.array(etao)                                                                           # Efficiency
    correctedMass = mdoto*(P03o/settings.P00)/np.sqrt(T2oabs/settings.T00)        # Corrected mass flow parameter

    """ checking if all efficiencies are greater than 0 and less or equal to 1 """
    if all(etao > 0) and all(etao <=1):       
        Nplot = Narr[iN]
        Nplot = str(Nplot *10**-3)+ ' Krpm'
        # Plotting pressure ratio 
        fig1 = plt.figure('Pressure ratio')
        plt.plot(mdoto, Pro, label=Nplot  ) # marker = markers[iN],
        plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 15)
        # plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.plot([0, 50], [1, 1], 'r--')
        plt.title(r'Pressure ratio')
        plt.grid(True)
        

        # Compressor map 
        fig99 = plt.figure('Compressor map')
        plt.plot(correctedMass, Pro, label=Nplot )# marker = markers[iN]
        plt.ylabel( r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m} \frac{P_{03}}{P_{00}} \sqrt{\frac{T_{00}}{T_{03}} }  $' + ' ' +'[kg/s]', fontsize = 15)
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title('Compressor map')
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        # plt.xlim([0, 50])
        plt.grid(True)

        # Mach number
        fig19 = plt.figure('Impeller exit mach number')
        plt.plot(mdoto, M2o, label= Nplot )# marker = markers[iN]
        plt.ylabel( r'${Ma}_2$', fontsize = 15)
        plt.xlabel(r'$\dot{m}$' + ' ' +'[kg/s]', fontsize = 15)
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title('Impeller mach number')

        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.grid(True)
        
        # Plotting efficiency curve 
        fig2 = plt.figure('Efficiency')
        plt.plot(mdoto, etao, label= Nplot  ) # marker = markers[iN],
        plt.grid()
        plt.ylabel(r'$\eta$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m} $'+ ' ' + '[kg/s]', fontsize = 15)
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.title(r'Efficiency ')
        plt.grid(True)
        
        fig4 = plt.figure('Exit Tip velocity')
        plt.plot(Narr[iN],U2o, 'ko', label= Nplot  ) # marker = markers[iN],
        plt.grid()
        plt.xlabel(r'$N $  [rpm]', fontsize = 15)
        plt.ylabel(r'$U_2 $ [m/s]', fontsize = 15)
        plt.title('Impeller exit tip speed')
        plt.grid(True)


        debug = 0
    else:
        print('   -> RPM  = ' + str(Narr[iN]) + ' has efficiency outside range [0, 1] !')

fig1 = plt.figure('Pressure ratio', label='Design pt.')
plt.plot(settings.mdot, settings.Pr, 'ro')


print('off-design_performance.py successfully run. \n')


plt.show()