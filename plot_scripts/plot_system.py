def plotSystemVariables(Compressor, IterationMatrix):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zsystem = [Ctheta2Mat, pressErrorMat, etaMat, MachExitArr, PrestMat, VslipMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """
    # Import
    import numpy as np
    import matplotlib.pyplot as plt

    systemVar = [Compressor.etaStage0, Compressor.lambda2, Compressor.iterTol, IterationMatrix.ZBarr, IterationMatrix.beta2BArr]
    designParam = [Compressor.r1, Compressor.r2, Compressor.rh1, Compressor.Ncrit, Compressor.Ndes]
    flowVar = [Compressor.Pr, Compressor.W1t, Compressor.Cm1, Compressor.U1t, Compressor.U2, Compressor.U2Crit]
    Zsystem = [IterationMatrix.WxMat, IterationMatrix.dh0SlipCorrMAt, IterationMatrix.sigmaMat, IterationMatrix.etaMat, IterationMatrix.PrestMat]                                                 # Goes to plotSystem.py


    WxMat = Zsystem[0]
    dh0SlipCorr = Zsystem[1]
    sigmaMat = Zsystem[2]
    etaMat = Zsystem[3]
    PressureEstimate = Zsystem[4]

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
    rhDivr1 = rh1 / rt1
    r1Dr2 = rt1 / rt2
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
    lvls = 20
    colorTheme = 'viridis'

    fig, axs21 = plt.subplots(2, 3)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad=7.0)
    fig.canvas.manager.set_window_title('System results')
    fig.subplots_adjust(top=0.9, bottom=0.09)

    # -------------- slip corrected compressor work plot -------------

    i=0
    j=2
    Z = (np.abs(dh0SlipCorr-WxMat)/WxMat)*100
    con = axs21[i, j].contourf(X, Y, Z, levels = 20, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    
    # contour5percent = axs21[i, j].contour(X, Y, Z,[0.03], colors=('k',),linestyles=('-',),linewidths=(1))
    # contourData = contour5percent.allsegs[0]
    # x_coords = [segment[:, 0] for segment in contourData]
    # y_coords = [segment[:, 1] for segment in contourData]
    # axs21[i, j].clabel(contour5percent, inline=True, fontsize=8)
    
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'Work deviation [%]' , fontsize=12)
    axs21[i, j].grid()

    Wmin = min(np.nanmin(WxMat), np.nanmin(dh0SlipCorr))
    Wmax = max(np.nanmax(WxMat), np.nanmax(dh0SlipCorr))


    # -------------- Work plot -------------
    i=0
    j=0
    Z = WxMat* 10**(-3)                  
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme, vmin=Wmin *10**-3, vmax=Wmax *10**-3)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'Work by efficiency corr. [KJ/kg] ' , fontsize=12)
    axs21[i, j].grid()
    # for xx, yy in zip(x_coords, y_coords):
    #     axs21[i, j].plot(xx, yy, 'k-')

    # -------------- Work plot -------------
    i=0
    j=1
    Z = dh0SlipCorr*10**-3              
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme, vmin=Wmin *10**-3, vmax=Wmax *10**-3)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r' $ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'Work by slip velocity corr. [KJ/kg]', fontsize=12)
    axs21[i, j].grid()
    # for xx, yy in zip(x_coords, y_coords):
    #     axs21[i, j].plot(xx, yy, 'k-')



        
    # -------------- Slip velocity plot -------------
    i=1
    j=1
    Z = sigmaMat         
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r'slip factor [-] ' , fontsize=12)
    axs21[i, j].grid()
    # for xx, yy in zip(x_coords, y_coords):
    #     axs21[i, j].plot(xx, yy, 'k-')

    # -------------- Textbox -------------
    i=1
    j=0
    Z = etaMat         
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r' Efficiency demand' , fontsize=12)
    axs21[i, j].grid()
    contourEfficiency = axs21[i, j].contour(X, Y, Z,np.arange(0, 1, 0.1), colors=('w',),linestyles=('-',),linewidths=(1))
    contourData = contourEfficiency.allsegs[0]
    axs21[i, j].clabel(contourEfficiency, inline=True, fontsize=8)

    # x_coords = [segment[:, 0] for segment in contourData]
    # y_coords = [segment[:, 1] for segment in contourData]
    # for xx, yy in zip(x_coords, y_coords):
    #     axs21[i, j].plot(xx, yy, 'k-')


    # -------------- Pressure ratio estimate -------------
    i=1
    j=2
    Z = PressureEstimate         
    con = axs21[i, j].contourf(X, Y, Z, levels = lvls, cmap=colorTheme)
    cbar = fig.colorbar(con, ax=axs21[i, j])
    cbar.ax.tick_params(labelsize=10)
    axs21[i, j].invert_yaxis()
    axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    axs21[i, j].set_title(r' Pressure ratio [-] ' , fontsize=12)
    axs21[i, j].grid()
    contourEfficiency = axs21[i, j].contour(X, Y, etaMat,np.arange(0, 1, 0.1), colors=('w',),linestyles=('-',),linewidths=(1))
    contourData = contourEfficiency.allsegs[0]
    axs21[i, j].clabel(contourEfficiency, inline=True, fontsize=8)
    contourPr124 = axs21[i, j].contour(X, Y, Z,[Compressor.Pr], colors=('w',),linestyles=('--',),linewidths=(1))
    contourData = contourPr124.allsegs[0]
    axs21[i, j].clabel(contourPr124, inline=True, fontsize=8)