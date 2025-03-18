import numpy as np
import matplotlib
matplotlib.use('tkagg')     # Use tkagg to avoid crashing when using X11-forwarding for plotting
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

def plotSystemVariables(Compressor, IterationMatrix):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zsystem = [Ctheta2Mat, pressErrorMat, etaMat, MachExitArr, PrestMat, VslipMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """

    systemVar = [Compressor.etaStage0, Compressor.lambda2, Compressor.iterTol, IterationMatrix.ZBarr, IterationMatrix.beta2BArr]
    designParam = [Compressor.r1, Compressor.r2, Compressor.rh1, Compressor.Ncrit, Compressor.Ndes]
    flowVar = [Compressor.Pr, Compressor.W1t, Compressor.Cm1, Compressor.U1t, Compressor.U2, Compressor.U2Crit]
    Zsystem = [IterationMatrix.WxMat, IterationMatrix.dh0SlipCorrMAt, IterationMatrix.sigmaMat, IterationMatrix.etaMat, IterationMatrix.PrestMat]

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

    x = ZBarr
    y = np.rad2deg(beta2bArr)
    X, Y = np.meshgrid(x, y)

    fig, axs = plt.subplots(2, 3, figsize=(15, 15))
    fig.tight_layout(pad=7.0)
    fig.canvas.manager.set_window_title('System results')
    fig.subplots_adjust(top=0.9, bottom=0.09)

    Wmin = min(np.nanmin(WxMat), np.nanmin(dh0SlipCorr))
    Wmax = max(np.nanmax(WxMat), np.nanmax(dh0SlipCorr))

    plot_configs = [
        # (row, col, Z, title, zlabel, vmin, vmax)
        (0, 2, (np.abs(dh0SlipCorr-WxMat)/WxMat)*100, 'Work deviation [%]', None, None, None),
        (0, 0, WxMat * 1e-3, 'Work by efficiency corr. [KJ/kg]', None, Wmin * 1e-3, Wmax * 1e-3),
        (0, 1, dh0SlipCorr * 1e-3, 'Work by slip velocity corr. [KJ/kg]', None, Wmin * 1e-3, Wmax * 1e-3),
        (1, 1, sigmaMat, 'Slip factor [-]', None, None, None),
        (1, 0, etaMat, 'Efficiency demand (not enthalpy based)', None, None, None),
        (1, 2, PressureEstimate, 'Pressure ratio [-]', None, None, None)
    ]

    def format_coord_factory(X, Y, Z):
        def format_coord(x, y):
            return f'x={x:.2f}, y={y:.2f}, z={get_z_value(X, Y, Z, x, y):.2f}'
        return format_coord

    for i, j, Z, title, zlabel, vmin, vmax in plot_configs:
        im = plot_contour(axs[i, j], X, Y, Z, title, zlabel, Compressor, vmin, vmax)
        
        if (i, j) == (1, 0) or (i, j) == (1, 2):
            contour = axs[i, j].contour(X, Y, etaMat, np.arange(0, 1, 0.1), colors='w', linewidths=1)
            axs[i, j].clabel(contour, inline=True, fontsize=8)
        
        if (i, j) == (1, 2):
            pr_contour = axs[i, j].contour(X, Y, Z, [Compressor.Pr], colors='w', linestyles='--', linewidths=1)
            axs[i, j].clabel(pr_contour, inline=True, fontsize=8)
        
        # Add cursor to each subplot
        cursor = Cursor(axs[i, j], useblit=True, color='red', linewidth=1)
        
        # Add hover function to display values
        axs[i, j].format_coord = format_coord_factory(X, Y, Z)


def plot_contour(ax, X, Y, Z, title, zlabel, Compressor, vmin=None, vmax=None):
    con = ax.contourf(X, Y, Z, levels=20, cmap='viridis', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(con, ax=ax)
    cbar.ax.tick_params(labelsize=10)
    if zlabel:
        cbar.set_label(zlabel, fontsize=10)
    
    ax.invert_yaxis()
    ax.set_xticks(np.arange(0, Compressor.bladeMax+1, 5))
    ax.set_xticklabels(np.arange(0, Compressor.bladeMax+1, 5), fontsize=10)
    ax.set_yticks(np.arange(-5, Compressor.beta2Bmax+1, -10))
    ax.set_yticklabels(np.arange(-5, Compressor.beta2Bmax+1, -10), fontsize=10)
    ax.set_xlabel(r'Blade number $Z_B$ ', fontsize=12)
    ax.set_ylabel(r'$ \beta _{2B}$ [deg]', fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.grid()
    return con

def get_z_value(X, Y, Z, x, y):
    """Get the z-value at a given x, y coordinate"""
    if X.min() <= x <= X.max() and Y.min() <= y <= Y.max():
        x_idx = np.abs(X[0] - x).argmin()
        y_idx = np.abs(Y[:, 0] - y).argmin()
        return Z[y_idx, x_idx]
    return np.nan