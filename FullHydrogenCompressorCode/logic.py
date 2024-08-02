### Import---------------------------------------------------------------------------------------------------
# import settings
# import geometryV7     # MSG: This runs the geometryV7.py file. Avoid this by including main function?
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

# import off_design_performance

# Set plot parameters------------------------------------------------------------------------------
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})
plt.close('all')
# ------------------------------

showGeometryPlots = plt.show()
# showGeometryPlots = 0


Npass = 0
Narr = np.arange(15000, geometryV8.Ndes , 5000)

# Specify point of interest for off-design
# beta2B = 8
# ZB = 5

# pressureRatioArray = []
# massFlowRateArray = []
# massFlowRateArray = []
# efficiencyMatrix = []
# markers = ['o','v','s','P', 'D', 'X']
Xs = []
Ys = []
etas = []
xMinArr = []
yMinArr = []
xMaxArr = []
yMaxArr = []

# ------------------- Plotting from geometry.py -------------------
plotCompressorParam(geometryV8.systemVar, geometryV8.Zcompressor, geometryV8.designParam, geometryV8.flowVar, geometryV8.text )
plotFlowConditions(geometryV8.systemVar, geometryV8.Zflow, geometryV8.designParam, geometryV8.flowVar, geometryV8.text)
plotVelocities(geometryV8.systemVar, geometryV8.Zvelocities, geometryV8.designParam, geometryV8.flowVar, geometryV8.text)
plotText(geometryV8.text)
plotSystemVariables(geometryV8.systemVar, geometryV8.Zsystem, geometryV8.designParam, geometryV8.flowVar, geometryV8.text)

print('r1: ' + str(round(geometryV8.r1, 3)) + 'm')
print('r2: ' + str(round(geometryV8.r2, 3)) + 'm' )

# ------------------- Plotting from off_design_performance.py -------------------
for iN in range(0, len(Narr)):
    Pro, P03o, T2oabs, mdoto, etao = odp.off_design_performance(Narr[iN])
    Pro = np.array(Pro)
    P03o = np.array(P03o)
    T2oabs = np.array(T2oabs)
    mdoto = np.array(mdoto)
    etao = np.array(etao)
    correctedMass = mdoto*(P03o/settingsOffDesign.P00)/np.sqrt(T2oabs/settingsOffDesign.T00)
    Xs.append(np.array(correctedMass))
    Ys.append(np.array(Pro))
    etas.append(np.array(etao))
    xMinArr.append(min(np.array(correctedMass)))
    yMinArr.append(min(np.array(Pro)))
    xMaxArr.append(max(np.array(correctedMass)))
    yMaxArr.append(max(np.array(Pro)))


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


        # # efficiency / Pressure ratio 
        # fig4 = plt.figure('Efficiency / Pressure ratio')
        # plt.plot(Pro,etao, label= str(Narr[iN]) )# marker = markers[iN]'
        # plt.ylabel(r'$\eta$', rotation=45, fontsize = 15)
        # plt.xlabel( r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        # plt.legend()
        # plt.plot([1, 1],[0, 50],  'r--')
        # plt.title(r'Efficiency / Pressure ratio')
        # plt.ylim([0, 1])
        # plt.grid(True)

        # # Temperature  
        # fig5 = plt.figure('Temperature @ outlet')
        # plt.plot(mdoto,T2oabs, label= str(Narr[iN]) )# marker = markers[iN]'
        # plt.ylabel(r'$T_{02}$', rotation=45, fontsize = 15)
        # plt.xlabel( r'$\dot{m}$', rotation=45, fontsize = 15)
        # plt.plot([0, 50], [settingsOffDesign.T00, settingsOffDesign.T00], 'k--')#, label='Inlet stang. Temp.')
        # plt.legend()
        # plt.xlim([5, 45])
        # plt.title(r'Stagnation temperature at outlet')
        # plt.grid(True)
    else:
        print('   -> RPM  = ' + str(Narr[iN]) + ' gives invalid efficiency !')



xMinVal = min(xMinArr)
yMinVal = min(yMinArr)
xMaxVal = max(xMaxArr)
yMaxVal = max(yMaxArr)

x = np.concatenate(Xs)
y = np.concatenate(Ys)
z = np.concatenate(etas)

gridX, gridY = np.meshgrid( np.linspace(xMinVal, xMaxVal, 100), np.linspace(yMinVal, yMaxVal, 100) )
gridZ = griddata((x, y), z, (gridX, gridY), method='cubic')

fig99 = plt.figure('Compressor map')
plt.contour(gridX, gridY, gridZ, levels = 14 , colors='black')
plt.plot()
plt.colorbar()









print('Off-design_performance.py successfully run. \n')


plt.show()