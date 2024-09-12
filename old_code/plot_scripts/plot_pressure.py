import math
import numpy as np
import matplotlib.pyplot as plt

def pressurePlot4GIF(Zpress, systemVar, text):


    pressErrorMat = Zpress[0]
    PrestMat = Zpress[1]
    VslipMat = Zpress[2]

    
    text1 = text[0]
    text2 = text[1]

    etaStage0 = systemVar[0]
    lambda2 = systemVar[1]
    iterTol = systemVar[2]
    ZBarr = systemVar[3]
    beta2bArr = systemVar[4]
    betamax = np.rad2deg(np.min(beta2bArr))
    betamin = np.rad2deg(np.max(beta2bArr))

    x = ZBarr                  # SAME FOR ALL COUNTOURS
    y = np.rad2deg(beta2bArr)     # SAME FOR ALL COUNTOURS
    
    X, Y = np.meshgrid(x, y)  
    lvls = 30
    colorTheme = 'Reds'

    fig, axs = plt.subplots(2, 2)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad=7.0)
    fig.suptitle(r'Flow properties for proposed  $\eta$ = ' + str(etaStage0) +  r' ,  $\lambda _{2}$ = ' + str(lambda2) + ' and Tolerance set to ' +str(iterTol*100) + '% .', y=0.98, x=0.38)
    fig.subplots_adjust(top=0.9, bottom=0.09)

        # ---------------- Pressure estimate error plot ---------------
    i=0 
    j=1
    Z = pressErrorMat
    con = axs[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs[i, j].invert_yaxis()
    axs[i, j].set_xticks(x[1::4])
    axs[i, j].set_xticklabels(x[1::4], fontsize=10)
    axs[i, j].set_yticks(np.arange(-5, -66, -5))
    axs[i, j].set_yticklabels(np.arange(-5, -66, -5), fontsize=10)
    axs[i, j].set_xlabel(r'Blade number $Z_B$ [deg]' , fontsize=12)
    axs[i, j].set_ylabel(r'$ \beta _{2B}$', fontsize=12)
    axs[i, j].set_title(r'Pressure estimate error contour plot [-]', fontsize=12)
    axs[i, j].grid()


    # -------------- PR plot -------------
    i=0
    j=0
    Z = PrestMat                  
    con = axs[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs[i, j].invert_yaxis()
    axs[i, j].set_xticks(x[1::4])
    axs[i, j].set_xticklabels(x[1::4], fontsize=10)
    axs[i, j].set_yticks(np.arange(betamax, betamin, 4)) #, fontsize=10)
    axs[i, j].set_yticklabels(np.arange(betamax, betamin, 4), fontsize=10)
    axs[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs[i, j].set_title(r'Estimated PR [-] ' , fontsize=12)
    axs[i, j].grid()



    # -------------- Textbox -------------
    i=1
    j=0
    axs[i, j].axis('off')  # Turn off the axis
    axs[i, j].text(0.0, 0.5, text1, ha='left', va='center', fontsize=12, linespacing = 1.8 )

    # -------------- Textbox -------------
    i=1
    j=1           
    axs[i, j].axis('off')  # Turn off the axis
    axs[i, j].text(0.0, 0.5, text2, ha='left', va='center', fontsize=12, linespacing =1.8 )

    