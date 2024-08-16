import numpy as np
import math
import matplotlib.pyplot as plt
U2 = 700
cm2 = 100
T0 = 293
beta = np.arange(0, 50, 5)
eta = 0.65
R = 8.314
Cp = 14310                                             
M = 2.01568 * 10**-3  
gamma = 1.41

Pr = 1 + ( ( (eta*M*(gamma-1)) / (gamma*R*T0 ) ) * U2*(U2-cm2*np.tan(np.deg2rad(beta))) )
print(Pr)

fig, ax = plt.subplots(1)
ax.plot(-beta, Pr, 'k-')
ax.invert_xaxis()
ax.set_yticklabels([])
ax.set_xticks(-beta)
ax.set_ylabel(r'Pressure ratio ', fontsize=12)
ax.set_xlabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
ax.grid()


plt.show()
















