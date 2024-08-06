import math
import numpy as np
import matplotlib.pyplot as plt

def plotVelocities(systemVar, Zvelocities, designParam, flowVar):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zcompressor = [NcritArr, rh1Mat, h2Mat, r1Mat, r2Mat, etaMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """


    c2Mat = Zvelocities[0]
    Ctheta2Mat = Zvelocities[1]
    c2mMat = Zvelocities[2]
    MachExitMat = Zvelocities[3]
    VslipMat = Zvelocities[4]

    etaStage0 = systemVar[0]
    lambda2 = systemVar[1]
    iterTol = systemVar[2]
    ZBarr = systemVar[3]
    beta2bArr = systemVar[4]

    betamax = np.rad2deg(np.min(beta2bArr))
    betamin = np.rad2deg(np.max(beta2bArr))

    rt1 = designParam[0]
    rt2 = designParam[1]
    rh1 = designParam[2]
    rhDivr1 = rh1/rt1
    r1Dr2 = rt1/rt2
    Ncrit = designParam[3]
    N = designParam[4]

    Pr = flowVar[0]
    W1t = flowVar[1]
    Cm1 = flowVar[2]
    U1t = flowVar[3]
    U2 = flowVar[4]
    U2Crit = flowVar[5]

   

    x = ZBarr                  # SAME FOR ALL COUNTOURS
    y = np.rad2deg(beta2bArr)     # SAME FOR ALL COUNTOURS

    X, Y = np.meshgrid(x, y)  
    lvls = 30
    colorTheme = 'Reds'
    
    fig, axs2 = plt.subplots(2, 3)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.suptitle(r'Velocities ', y=0.98, x=0.4)
    fig.tight_layout(pad=7.0)
    fig.subplots_adjust(top=0.9, bottom=0.09)


    # -------------- C2 -------------
    i=0
    j=0
    Z = c2Mat                  
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks([1, 10, 20 , 30, 40 ,50])
    axs2[i, j].set_xticklabels([1, 10, 20 , 30, 40 ,50], fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, -66, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, -66, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Absolute discharge velocity [m/s]', fontsize=12)
    axs2[i, j].grid()





    # -------------- C2m -------------
    i=0
    j=1
    Z = c2mMat
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks([1, 10, 20 , 30, 40 ,50])
    axs2[i, j].set_xticklabels([1, 10, 20 , 30, 40 ,50], fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, -66, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, -66, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Meridonal discharge velocity [m/s] ' , fontsize=12)
    axs2[i, j].grid()

    # -------------- Angular velocity -------------
    i=0
    j=2
    Z = Ctheta2Mat
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks([1, 10, 20 , 30, 40 ,50])
    axs2[i, j].set_xticklabels([1, 10, 20 , 30, 40 ,50], fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, -66, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, -66, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Angular discharge velocity [m/s] ' , fontsize=12)
    axs2[i, j].grid()

    # -------------- Mach number velocity -------------
    i=1
    j=0
    Z = MachExitMat
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks([1, 10, 20 , 30, 40 ,50])
    axs2[i, j].set_xticklabels([1, 10, 20 , 30, 40 ,50], fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, -66, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, -66, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Discharge Mach number [-] ' , fontsize=12)
    axs2[i, j].grid()

    # -------------- Slip velocity -------------
    i=1
    j=1
    Z = VslipMat
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks([1, 10, 20 , 30, 40 ,50])
    axs2[i, j].set_xticklabels([1, 10, 20 , 30, 40 ,50], fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, -66, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, -66, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Slip velocity [m/s] ' , fontsize=12)
    axs2[i, j].grid()

    # -------------- Textbox -------------
    i=1
    j=2          
    axs2[i, j].axis('off')  # Turn off the axis
    axs2[i, j].text(0.0, 0.5, "Vacant plot space", ha='left', va='center', fontsize=12, linespacing = 1.8 )
    # Ready to be filled by plot


  
