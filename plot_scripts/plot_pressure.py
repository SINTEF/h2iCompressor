import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

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

    x = ZBarr
    y = np.rad2deg(beta2bArr)
    
    X, Y = np.meshgrid(x, y)  
    lvls = 30
    colorTheme = 'Reds'

    fig, axs = plt.subplots(2, 2)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad=7.0)
    fig.suptitle(r'Flow properties for proposed  $\eta$ = ' + str(etaStage0) +  r' ,  $\lambda _{2}$ = ' + str(lambda2) + ' and Tolerance set to ' +str(iterTol*100) + '% .', y=0.98, x=0.38)
    fig.subplots_adjust(top=0.9, bottom=0.09)

    def format_coord_factory(X, Y, Z):
        def format_coord(x, y):
            return f'x={x:.2f}, y={y:.2f}, z={get_z_value(X, Y, Z, x, y):.2f}'
        return format_coord

    plot_configs = [
        (0, 1, pressErrorMat, r'Pressure estimate error contour plot [-]'),
        (0, 0, PrestMat, r'Estimated PR [-]')
    ]

    for i, j, Z, title in plot_configs:
        con = axs[i, j].contourf(X, Y, Z, levels=lvls, cmap=colorTheme)
        cbar = fig.colorbar(con, ax=axs[i, j])
        cbar.ax.tick_params(labelsize=10)
        axs[i, j].invert_yaxis()
        axs[i, j].set_xticks(x[1::4])
        axs[i, j].set_xticklabels(x[1::4], fontsize=10)
        axs[i, j].set_yticks(np.arange(-5, -66, -5))
        axs[i, j].set_yticklabels(np.arange(-5, -66, -5), fontsize=10)
        axs[i, j].set_xlabel(r'Blade number $Z_B$ [deg]', fontsize=12)
        axs[i, j].set_ylabel(r'$ \beta _{2B}$', fontsize=12)
        axs[i, j].set_title(title, fontsize=12)
        axs[i, j].grid()

        cursor = Cursor(axs[i, j], useblit=True, color='red', linewidth=1)
        axs[i, j].format_coord = format_coord_factory(X, Y, Z)

    # -------------- Textbox -------------
    axs[1, 0].axis('off')
    axs[1, 0].text(0.0, 0.5, text1, ha='left', va='center', fontsize=12, linespacing=1.8)

    # -------------- Textbox -------------
    axs[1, 1].axis('off')
    axs[1, 1].text(0.0, 0.5, text2, ha='left', va='center', fontsize=12, linespacing=1.8)

    plt.show()
def get_z_value(X, Y, Z, x, y):
    if X.min() <= x <= X.max() and Y.min() <= y <= Y.max():
        x_idx = np.abs(X[0] - x).argmin()
        y_idx = np.abs(Y[:, 0] - y).argmin()
        return Z[y_idx, x_idx]
    return np.nan