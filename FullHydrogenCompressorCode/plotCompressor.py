import math
import numpy as np
import matplotlib.pyplot as plt
import settings

def plotCompressorParam(systemVar, Zcompressor, designParam, flowVar):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zcompressor = [NcritArr, rh1Mat, h2Mat, r1Mat, r2Mat, etaMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """



    h2Mat = Zcompressor[0]
    VslipMat = Zcompressor[1]
    beta2FlowMat = Zcompressor[2]


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
    
    fig, axs2 = plt.subplots(2, 3)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.canvas.manager.set_window_title('compressor')
    fig.tight_layout(pad=7.0)
    fig.subplots_adjust(top=0.9, bottom=0.09)


    # -------------- Outlet exit height plot -------------
    i=0
    j=1
    Z = h2Mat                  
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks(np.arange(0, settings.bladeMax+1, 5))
    axs2[i, j].set_xticklabels(np.arange(0, settings.bladeMax+1, 5), fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, settings.beta2Bmax+1, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, settings.beta2Bmax+1, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Impeller cylinder height [m]', fontsize=12)
    axs2[i, j].grid()





    # -------------- Slip velocity -------------
    i=0
    j=0
    Z = VslipMat                  
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks(np.arange(0, settings.bladeMax+1, 5))
    axs2[i, j].set_xticklabels(np.arange(0, settings.bladeMax+1, 5), fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, settings.beta2Bmax+1, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, settings.beta2Bmax+1, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Slip velocity [m/s] ' , fontsize=12)
    axs2[i, j].grid()

    # -------------- Outlet flow angle -------------
    i=0
    j=2
    Z = beta2FlowMat                  
    con = axs2[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs2[i, j])
    cbar.ax.tick_params(labelsize=10)

    contour10s = axs2[i, j].contour(X, Y, Z,[10, 20, 30, 40, 50, 60, 70], colors=('k',),linestyles=('-',),linewidths=(1))
    contourData = contour10s.allsegs[0]
    x_coords = [segment[:, 0] for segment in contourData]
    y_coords = [segment[:, 1] for segment in contourData]
    axs2[i, j].clabel(contour10s, inline=True, fontsize=8)

    axs2[i, j].invert_yaxis()
    axs2[i, j].set_xticks(np.arange(0, settings.bladeMax+1, 5))
    axs2[i, j].set_xticklabels(np.arange(0, settings.bladeMax+1, 5), fontsize=10)
    axs2[i, j].set_yticks(np.arange(-5, settings.beta2Bmax+1, -10))
    axs2[i, j].set_yticklabels(np.arange(-5, settings.beta2Bmax+1, -10), fontsize=10)
    axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs2[i, j].set_title(r'Outlet flow angle $ \beta _{2}$ [m/s] ' , fontsize=12)
    axs2[i, j].grid()


