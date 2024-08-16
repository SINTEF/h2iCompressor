import matplotlib.pyplot as plt
import numpy as np
import math

r1Divr2 = np.arange(0, 1.1, 0.01)
ZB = np.arange(5, 30, 5)
beta2b = np.arange(-40, -19, 10)
sigma = np.zeros((r1Divr2.shape))
sigmaW = np.zeros(r1Divr2.shape)

limit = 0
ylims = [ [0, 0.9], [0, 2.2], [0, 7.5], [0, 20], [0, 41]  ]
for k in range(0, len(beta2b), 1):
    fig, axs = plt.subplots(len(ZB))
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad=5.0)
    fig.subplots_adjust(top=0.9, bottom=0.05)
    b = beta2b[k]

    axs[0].text(0.0, 1.3, r'$\beta_{2B} = $' + str(b) +' ' r'$ deg$', fontsize = 18)

    for j in range(len(ZB)):
        zb = ZB[j]
        for i in range(len(r1Divr2)):
            sigmaW[i] = 1 - (math.sqrt(math.cos(np.deg2rad(b)) ) )/(zb **0.7)
            limit = math.exp(-8.16 * math.cos(np.deg2rad(b))/zb)
            sigma[i] = sigmaW[i]* ( 1 - ( (r1Divr2[i] - limit)/( 1 - limit ) )**3 )
        title = r'$\underline{\beta_{2B}} =$' + ' ' +str(b) + r'$ deg$'  + r',  $Z_B=$ ' + str(zb)
        title = r'  $Z_B=$ ' + str(zb)
        axs[j].plot(r1Divr2, sigmaW, 'b')
        axs[j].plot(r1Divr2, sigma , 'r')
        axs[j].plot([limit, limit], [-20, 50], 'g--')#[min(min(sigma), min(sigmaW)), max(max(sigma), max(sigmaW))], 'g--')

        fig.legend([r'$\sigma_{w}$',r'$\sigma_{corr}$',  r'$\epsilon_{limit}$'], fontsize=18, ncol=3, loc='upper center',
                    bbox_to_anchor=(0.7, 1), fancybox=True, framealpha=0)



        # axs[j].legend([r'$\sigma$', r'$\sigma_{w}$', r'$\epsilon_{limit}$'], fontsize=12, loc='upper center', bbox_to_anchor=(1.25, 0.5))

        axs[j].set_ylabel(r'slip factor', fontsize = 13.5)
        yl = ylims[j]
        axs[j].axis([min(r1Divr2), max(r1Divr2), 0, 1 ]) #, -1, 5])
        axs[j].set_xticks(np.arange(min(r1Divr2), 1.05, 0.05))
        #axs[j].set_title(title, fontsize = 15) 
        axs[j].set_title(title, fontsize = 15, loc='left') 
        axs[j].grid()
        axs[j].axvspan(0.35, 0.85, facecolor='g', alpha=0.1)

    axs[j].set_xlabel(r'$r_{1} / r_{2}$', fontsize = 13.5)
plt.show(block=True)

















