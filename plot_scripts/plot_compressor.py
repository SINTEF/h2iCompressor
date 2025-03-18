import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

def plotCompressorParam(Compressor, IterationMatrix):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zcompressor = [NcritArr, rh1Mat, h2Mat, r1Mat, r2Mat, etaMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """
    systemVar = [Compressor.etaStage0, Compressor.lambda2, Compressor.iterTol, IterationMatrix.ZBarr, IterationMatrix.beta2BArr]
    designParam = [Compressor.r1, Compressor.r2, Compressor.rh1, Compressor.Ncrit, Compressor.Ndes]
    flowVar = [Compressor.Pr, Compressor.W1t, Compressor.Cm1, Compressor.U1t, Compressor.U2, Compressor.U2Crit]
    Zcompressor = [IterationMatrix.b2Mat, IterationMatrix.VslipMat, IterationMatrix.beta2flowMat]

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
    
    fig, axs2 = plt.subplots(2, 3)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.canvas.manager.set_window_title('compressor')
    fig.tight_layout(pad = 7.0)
    fig.subplots_adjust(top = 0.9, bottom = 0.09)

    def format_coord_factory(X, Y, Z):
        def format_coord(x, y):
            return f'x={x:.2f}, y={y:.2f}, z={get_z_value(X, Y, Z, x, y):.2f}'
        return format_coord

    plot_configs = [
        (0, 1, h2Mat, r'Impeller cylinder height [m]'),
        (0, 0, VslipMat, r'Slip velocity [m/s]'),
        (0, 2, beta2FlowMat, r'Outlet flow angle $ \beta _{2}$ [deg]')
    ]

    for i, j, Z, title in plot_configs:
        con = axs2[i, j].contourf(X, Y, Z, levels=lvls, cmap=colorTheme)
        cbar = fig.colorbar(con, ax=axs2[i, j])
        cbar.ax.tick_params(labelsize=10)
        axs2[i, j].invert_yaxis()
        axs2[i, j].set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
        axs2[i, j].set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
        axs2[i, j].set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
        axs2[i, j].set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
        axs2[i, j].set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
        axs2[i, j].set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
        axs2[i, j].set_title(title, fontsize=12)
        axs2[i, j].grid()

        cursor = Cursor(axs2[i, j], useblit=True, color='red', linewidth=1)
        axs2[i, j].format_coord = format_coord_factory(X, Y, Z)



def get_z_value(X, Y, Z, x, y):
    if X.min() <= x <= X.max() and Y.min() <= y <= Y.max():
        x_idx = np.abs(X[0] - x).argmin()
        y_idx = np.abs(Y[:, 0] - y).argmin()
        return Z[y_idx, x_idx]
    return np.nan