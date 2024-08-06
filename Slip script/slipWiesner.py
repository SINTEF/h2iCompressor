import matplotlib.pyplot as plt
import numpy as np
import math

# plt.rc('text', usetex=True)
plt.rc('font', family='serif')

r1Divr2 = np.arange(0, 1.1, 0.1)

ZB = np.arange(0, 31, 1)
beta2b = np.arange(-20, -41, -10)
sigma = np.zeros((r1Divr2.shape))
sigmaW = np.zeros(ZB.shape)
sigW = np.zeros(( np.shape(ZB)[0], np.shape(beta2b)[0] ) )

limit = sigmaW
ylims = [ [0, 0.47], [0, 1.3], [0, 5.7], [0, 19]  ]
color = ['r-', 'b-', 'k-', 'm-']
colorLimit = ['r--', 'b--', 'k--', 'm--']
fig, axs = plt.subplots(len(beta2b))
fig.set_figwidth(15)
fig.set_figheight(9)


fig.tight_layout(pad=6)

for k in range(0, len(beta2b), 1):
    
    b = beta2b[k]

    for j in range(len(ZB)):
        zb = ZB[j]
        sigmaW[j] = 1 - (math.sqrt(math.cos(np.deg2rad(b)) ) )/(zb **0.7)
        sigW[j, k] = sigmaW[j]
    title = r'$\beta_{2B} =$' + str(b) + r'$ deg$'
    axs[k].plot(ZB, sigmaW, '-bo' )

    axs[k].set_xlabel(r'$Z_B$', fontsize = 13.5)
    axs[k].set_ylabel(r'$\sigma _w$',  fontsize = 13.5)

    # axs[j].axis([min(r1Divr2), max(r1Divr2), yl[0], yl[1] ]) #, -1, 5])
    axs[k].set_xticks(np.arange(1, 25, 1))
    axs[k].set_yticks(np.arange(0.4, 1, 0.1))
    fig.legend([r'$\sigma_{w}$'], fontsize=15, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1))

    axs[k].set_title(title, fontsize = 15)
    axs[k].grid(True)
        
figW ,axW  = plt.subplots(figsize=(10, 6))
plt.grid()
shapes = ['-bo', '-b^', '-bs']
legendW = [r'$\beta_{2B} =$' + str(beta2b[0]) + r'$ deg$', r'$\beta_{2B} =$' + str(beta2b[1]) + r'$ deg$', r'$\beta_{2B} =$' + str(beta2b[2]) + r'$ deg$'] 
plt.legend(legendW)
for k in range(np.shape(sigW)[1]):
    # axW.plot(ZB, sigW[:, k], '-b') #, label=legendW[k])
    axW.plot(ZB, sigW[:, k], shapes[k], markersize=5.6, markerfacecolor = 'red', markeredgecolor='black') 


axW.legend(legendW, fontsize=15, loc='lower right')
axW.set_xlabel(r'$Z_B$', fontsize = 13.5)
axW.set_ylabel(r'$\sigma _w$',  fontsize = 13.5)
axW.axis([.8, 30.2, 0.02, 1])
axW.set_xticks(np.arange(1, np.shape(ZB)[0], 1))
axW.set_yticks(np.arange(0.025, 1, 0.05))
axW.axvspan(5, 25, facecolor='g', alpha=0.1)
axW.annotate(r'Practical range for $Z_B$', fontsize = 13.5, xy=(10.2, 0.229))


plt.show(block=True)

















