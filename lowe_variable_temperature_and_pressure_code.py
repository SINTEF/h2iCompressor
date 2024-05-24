"""
The following Python code was used to calculate the lines of constant power, the effect of controlling inlet
conditions and the plot of the 100kW line of constant power with mass flow rate in kg/s on the x-axis
for Figure 9 and Figure 10. This code was also used to recreate the HE351Ve compressor map.

Reference: Lowe, Mitchell (2016). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244

Notes:
    - When increasing inlet temperature, we do not account for the associated reduction in air density
    - When increasing inlet pressure, we do not account for the associated increase in air density


Author: Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


# Import-------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt


# Variables---------------------------------------------------------------------------------------
cp = 1.006                  # Specific heat at constant pressure [kJ]
y = 1.4                     # Ratio of specific heats [-]
eta = 0.75                  # Assumed efficiency [-]
comp_power = [20, 100]      # Compressor power [kW]
speed = [120]               # Shaft Speed [krpm]
T1 = 298                    # Inlet stagnation temperature [K]
P1 = 0.1                    # Inlet pressure [MPa]
mfp_arr = np.linspace(0.0001, 160, 100)    # Mass flow parameter [10^(- 6) s K^0.5]


# Pressure ratio as a function of mass flow parameter, Equation [8]------------------------------
def pressure_ratio(mfp, eta, comp_power, P1, T1, y, cp):
    return (eta * comp_power / (mfp * cp * P1 * T1 ** 0.5) + 1) ** (y / (y - 1))


# Plot requirement lines for different compressor power--------------------------------------------
plt.figure()
for i in range(len(comp_power)):
    Pr = pressure_ratio(mfp_arr, eta, comp_power[i], P1, T1, y, cp)
    plt.plot(mfp_arr, Pr, linewidth = 2, label = str(comp_power[i]) + ' kW')
plt.xlim(0, 160)
plt.ylim(0, 5.0)
plt.xlabel('Mass flow parameter (MFP)', fontsize = 15)
plt.ylabel('Total to total pressure ratio', fontsize = 15)
plt.title('Requirement lines', fontsize = 15)
plt.legend()
#plt.show()


# Plot pressure ratio for variable temperature (not using defined T1 variable)-------------------
T1_arr = [298, 323, 348, 600]
plt.figure()
for i in range(len(T1_arr)):
    Pr = pressure_ratio(mfp_arr, eta, 100, P1, T1_arr[i], y, cp)
    plt.plot(mfp_arr, Pr, linewidth = 2, label = str(T1_arr[i]) + ' K')
plt.xlim(0, 140)
plt.ylim(0, 5.0)
plt.xlabel('Mass flow parameter (MFP)', fontsize = 15)
plt.ylabel('Total to total pressure ratio', fontsize = 15)
plt.title('Variable inlet temperature', fontsize = 15)
plt.legend()
#plt.show()


# Plot pressure ratio for variable pressure (not using defined P variable)---------------------
P_arr = [0.1, 0.125, 0.15, 0.175]
plt.figure()
for i in range(len(P_arr)):
    Pr = pressure_ratio(mfp_arr, eta, 100, P_arr[i], T1, y, cp)
    plt.plot(mfp_arr, Pr, linewidth = 2, label = str(P_arr[i]) + ' MPa')
plt.xlim(0, 140)
plt.ylim(0, 5.0)
plt.xlabel('Mass flow parameter (MFP)', fontsize = 15)
plt.ylabel('Total to total pressure ratio', fontsize = 15)
plt.title('Variable inlet pressure', fontsize = 15)
plt.legend()
#plt.show()


# Plot pressure ratio as a function of mass flow rate in kg/s----------------------------------
def pressure_ratio_SI_units(mdot, eta, comp_power, T1, y, cp):
    return (eta * comp_power / (mdot * cp * T1) + 1) ** (y / (y - 1))
mdot = np.linspace(0.0001, 1.4, 100)
plt.figure()
for i in range(len(comp_power)):
    Pr = pressure_ratio_SI_units(mdot, eta, comp_power[i], T1, y, cp)
    plt.plot(mdot, Pr, linewidth = 2, label = str(comp_power[i]) + ' kW')
plt.xlim(0, 1.4)
plt.ylim(0, 5.0)
plt.xlabel('Mass flow rate [kg/s]', fontsize = 15)
plt.ylabel('Total to total pressure ratio', fontsize = 15)
plt.legend()
#plt.show()


# Plotting Holset HE351Ve compressor map----------------------------------------------------
plt.figure()
plt.plot([25.10948472, 31.16982327, 38.02921278, 44.41023875, 51.27211801, 57.499736], [2.101782025, 2.074240887, 2.000797853, 1.927354818, 1.789649129, 1.560139647], 'k-o',
         [36.71872019, 43.89559601, 51.55154674, 59.68870643, 67.19089349, 74.22467601], [2.900475022, 2.863753505, 2.808671229, 2.680145919, 2.478177575, 2.019158611], 'k-o', 
         [52.15024017, 56.93529829, 62.99670388, 69.21969805, 75.44696036, 80.90113623], [3.81851295, 3.781791433, 3.726709157, 3.616544606, 3.396215503, 2.55162061], 'k-o',
         [13.01192671, 18.59319036, 23.8569677, 28.96129053, 33.90651452, 39.65078947], [1.560139647, 1.550959268, 1.505057372, 1.459155475, 1.4040732, 1.303089027], 'k-o',
         [56.90506563, 62.00796574, 68.38685764, 73.4922475, 79.5600553, 85.18435604], [4.562123672, 4.552943293, 4.534582534, 4.4611395, 4.240810397, 3.120804125], 'k-o',
         [12.69337336, 24.79093137, 36.55926568, 51.51277779, 56.74561111], [1.550959268, 2.092601646, 2.900475022, 3.809332571, 4.562123672], 'k-o',
         [39.49133496, 57.34028149, 74.22467601, 80.90078055, 85.02419017], [1.303089027, 1.560139647, 2.019158611, 2.560800989, 3.139164883], 'k-o')
plt.annotate('126.86 RPS/K', xy = (56.74561111, 4.562123672), xytext = (50, 4.7), fontsize = 15)
plt.annotate('114.17 RPS/K', xy = (51.51277779, 3.809332571), xytext = (38, 4), fontsize = 15)
plt.annotate('97.58 RPS/K', xy = (36.55926568, 2.900475022), xytext = (27, 3), fontsize = 15)
plt.annotate('78.07 RPS/K', xy = (24.79093137, 2.092601646), xytext = (15, 2.2), fontsize = 15)
plt.annotate('58.55 RPS/K', xy = (12.69337336, 1.550959268), xytext = (10, 1.25), fontsize = 15)
plt.xlim(0, 100)
plt.ylim(0, 5.0)
plt.xlabel('Mass flow parameter (MFP)', fontsize = 15)
plt.ylabel('Total to total pressure ratio', fontsize = 15)
plt.title('Holset HE351Ve compressor map', fontsize = 15)
plt.show()