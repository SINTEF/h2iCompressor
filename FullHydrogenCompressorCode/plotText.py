import math
import numpy as np
import matplotlib.pyplot as plt

def plotText(text ):
    # ------------- PLOTTING ------------ 
    """
    systemVar = [etaStage0, lambda2, iterTol, Zbarr, beta2bArr]
    Zsystem = [Ctheta2Mat, pressErrorMat, etaMat, MachExitArr, PrestMat, VslipMat, trueFalseMat]
    designParam = [r1, rt2, rh, Ncrit, N]
    flowVar = [PR, W1t, Cm1, U1t, U2, U2Crit]
    """

    text1 = text[0]
    text2 = text[1]
    text3 = text[2]

   
    fig, axs21 = plt.subplots(1, 3)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad=7.0)
    fig.suptitle(r'System output data', y=0.98, x=0.38)
    fig.subplots_adjust(top=0.9, bottom=0.09)

  
    # -------------- Textbox -------------
    i=0
    j=0
    axs21[j].axis('off')  # Turn off the axis
    axs21[j].text(0.0, 0.5, text3, ha='left', va='bottom', fontsize=12, linespacing = 1.8 )
    # -------------- Textbox -------------
    i=0
    j=1
    axs21[j].axis('off')  # Turn off the axis
    axs21[j].text(0.0, 0.5, text1, ha='left', va='center', fontsize=12, linespacing = 1.8 )

    # -------------- Textbox -------------
    i=0
    j=2           
    axs21[j].axis('off')  # Turn off the axis
    axs21[j].text(0.0, 0.5, text2, ha='left', va='center', fontsize=12, linespacing =1.8 )
