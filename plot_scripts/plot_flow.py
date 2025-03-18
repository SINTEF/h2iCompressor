import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

def plotFlowConditions(Compressor, IterationMatrix):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zflow = [Ctheta2Mat, pressErrorMat, etaMat, MachExitArr, PrestMat, VslipMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """
    
    systemVar = [Compressor.etaStage0, Compressor.lambda2, Compressor.iterTol, IterationMatrix.ZBarr, IterationMatrix.beta2BArr]
    designParam = [Compressor.r1, Compressor.r2, Compressor.rh1, Compressor.Ncrit, Compressor.Ndes]
    flowVar = [Compressor.Pr, Compressor.W1t, Compressor.Cm1, Compressor.U1t, Compressor.U2, Compressor.U2Crit]
    Zflow = [IterationMatrix.Ctheta2Mat, IterationMatrix.pressErrorMat, IterationMatrix.etaMat, IterationMatrix.MachExitMat, IterationMatrix.PrestMat, IterationMatrix.WxMat, IterationMatrix.dh0SlipCorrMAt]
    
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

    x = ZBarr
    y = np.rad2deg(beta2bArr)

    X, Y = np.meshgrid(x, y)  
    lvls = 30
    colorTheme = 'viridis'

    fig, axs21 = plt.subplots(2, 3)
    fig.set_figwidth(30)
    fig.set_figheight(15)
    fig.tight_layout(pad = 7.0)
    fig.canvas.manager.set_window_title('flowResults')
    fig.subplots_adjust(top = 0.9, bottom = 0.09)

    def format_coord_factory(X, Y, Z):
        def format_coord(x, y):
            return f'x={x:.2f}, y={y:.2f}, z={get_z_value(X, Y, Z, x, y):.2f}'
        return format_coord

    plot_configs = [
        (1, 2, (np.abs(dh0SlipCorr-WxMat)/WxMat)*100, r'Slip corrected compressor work [KJ/kg]', None),
        (0, 1, PrestMat, r'Estimated PR [-]', None),
        (0, 2, pressErrorMat, r'Pressure estimate error contour plot [-]', None),
        (0, 0, etaMat, r'Efficiency [-]', None)
    ]

    for i, j, Z, title, zlabel in plot_configs:
        con = axs21[i, j].contourf(X, Y, Z, levels=20, cmap=colorTheme)
        cbar = fig.colorbar(con, ax=axs21[i, j])
        cbar.ax.tick_params(labelsize=10)
        axs21[i, j].invert_yaxis()
        axs21[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
        axs21[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
        axs21[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
        axs21[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
        axs21[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
        axs21[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
        axs21[i, j].set_title(title, fontsize=12)
        axs21[i, j].grid()

        cursor = Cursor(axs21[i, j], useblit=True, color='red', linewidth=1)
        axs21[i, j].format_coord = format_coord_factory(X, Y, Z)


def get_z_value(X, Y, Z, x, y):
    if X.min() <= x <= X.max() and Y.min() <= y <= Y.max():
        x_idx = np.abs(X[0] - x).argmin()
        y_idx = np.abs(Y[:, 0] - y).argmin()
        return Z[y_idx, x_idx]
    return np.nan