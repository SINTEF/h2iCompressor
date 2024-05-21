"""
The following Python code was used to calculate the lines of constant power, the effect of controlling inlet
conditions and the plot of the 100kW line of constant power with mass flow rate in kg/s on the x-axis
for Figure 9 and Figure 10. This code was also used to recreate the HE351Ve compressor map. Each
section of the code was enabled and disabled by using three quotation marks (“””) on the top and bottom
of the required section. For example, when plotting the varying pressure over the Holset compressor
map, the requirement lines plot, variable temperature and the standardised mass flow rate sections were
all disabled.


Reference: Lowe, Mitchell (2016). Design of a load dissipation device for a 100 kW supercritical CO2 turbine, https://doi.org/10.14264/uql.2017.244


Author: Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""



# Import----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')


# Variables----------------------------------------------------------------------
comp_power = [20, 100]      # Compressor power [kW]
speed = [120]               # Shaft Speed [krpm]
T1 = 298                    # Inlet stagnation temperature [K]
cp = 1.005                  # Specific heat at constant pressure [kJ]
y = 1.4                     # Ratio of specific heats
P = 0.1                     # Inlet pressure [MPa]
eta = 0.75                  # Assumed efficiency


# Equation [8]-------------------------------------------------------------------
def pressure_ratio(mfp, eta, comp_power, P, T1, y, cp):
    return (eta * comp_power / (mfp * cp * P * T1 ** 0.5) + 1) ** (y / (y - 1))


# Requirement lines plot-------------------------------------------------------
mfp = np.arange(0, 160, 0.5)
plt.figure()
plt.xlim(0, 160)
plt.ylim(1, 5.0)
for i in range(2):
    Pr = pressure_ratio(mfp, eta, comp_power[i], P, T1, y, cp)
plt.plot(mfp, Pr)
plt.legend(('20 kW', '100 kW' ))
plt.show()


# Calculate variable temperature ----------------------------------------------
T0 = [298, 323, 348, 600]
for i in range (4):
    Pr = pressure_ratio(mfp, eta, 100, P, T1[i], y, cp)     # MSG: Change T1 to list or to T0 or convert T0 to T1?
plt.plot(mfp, Pr)
plt.xlim(0, 140)
plt.ylim(1, 5.0)
plt.xlabel('Mass Flow Parameter (MFP)', fontsize = 15)
plt.ylabel('Total to Total pressure ratio', fontsize = 15)
plt.title('Variable Inlet Temperature', fontsize = 20)
plt.legend(('298 K', '323 K', '348 K', '600 K'))