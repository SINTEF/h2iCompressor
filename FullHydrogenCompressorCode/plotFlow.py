import math
import numpy as np
import matplotlib.pyplot as plt
import settingsOffDesign
def plotFlowConditions(systemVar, Zflow, designParam, flowVar ):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zflow = [Ctheta2Mat, pressErrorMat, etaMat, MachExitArr, PrestMat, VslipMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """

    Ctheta2Mat = Zflow[0]
    pressErrorMat = Zflow[1]
    etaMat = Zflow[2] 
    MachExitMat = Zflow[3]
    PrestMat = Zflow[4]
    WxMat = Zflow[5]
    dh0SlipCorr =  Zflow[6]

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
    colorTheme = 'viridis'

    fig, axs21 = plt.subplots(2,3 )
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad=7.0)
    fig.canvas.manager.set_window_title('flowResults')
    fig.subplots_adjust(top=0.9, bottom=0.09)

    # -------------- slip corrected compressor work plot -------------

    i=1
    j=2
    Z = dh0SlipCorr*(10**-3)                 
    Z = np.abs(dh0SlipCorr-WxMat)/WxMat
    con = axs21[i, j].contourf(X, Y, Z, levels = 20, cmap=colorTheme)

    contour5percent = axs21[i, j].contour(X, Y, Z,[0.03], colors=('k',),linestyles=('-',),linewidths=(1))
    contourData = contour5percent.allsegs[0]
    x_coords = [segment[:, 0] for segment in contourData]
    y_coords = [segment[:, 1] for segment in contourData]
    axs21[i, j].clabel(contour5percent, inline=True, fontsize=8)

    # cbar = fig.colorbar(con, ax=axs21[i, j])
    # cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, settingsOffDesign.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, settingsOffDesign.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'Slip corrected compressor work [KJ/kg]' , fontsize=12)
    axs21[i, j].grid()
    axs21[i, j].cla()  # Turn off the axis

    # -------------- PR plot -------------
    i=0
    j=1
    Z = PrestMat                  
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, settingsOffDesign.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, settingsOffDesign.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'Estimated PR [-] ' , fontsize=12)
    axs21[i, j].grid()
    for xx, yy in zip(x_coords, y_coords):
        axs21[i, j].plot(xx, yy, 'k-')


    # ---------------- Pressure estimate error plot ---------------
    i=0 
    j=2
    Z = pressErrorMat
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, settingsOffDesign.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, settingsOffDesign.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ [deg]' , fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$', fontsize=12)
    axs21[i, j].set_title(r'Pressure estimate error contour plot [-]', fontsize=12)
    axs21[i, j].grid()
    for xx, yy in zip(x_coords, y_coords):
        axs21[i, j].plot(xx, yy, 'k-')


    # -------------- Efficiency plot -------------
    i=0
    j=0
    Z = Zflow[2]                  
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, settingsOffDesign.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, settingsOffDesign.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, settingsOffDesign.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'Efficiency [-]', fontsize=12)
    axs21[i, j].grid()
    for xx, yy in zip(x_coords, y_coords):
        axs21[i, j].plot(xx, yy, 'k-')






