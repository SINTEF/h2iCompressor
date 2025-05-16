"""
Python script for plotting the conditions after optimization of the geometry for different blade numbers
and blade angles. 

Author(s): Martin Spillum Grønli (SINTEF Energy Research, 2025)
"""

# Import
import numpy as np
import matplotlib.pyplot as plt


def plot_dimensions(Compressor, InletConditions, IterationMatrix):
    """ Plot dimensions of compressor stage after geometry optimization """
    
    ZB = IterationMatrix.ZB
    beta2B = np.rad2deg(IterationMatrix.beta2B)
    r1Divr2 = Compressor.r1 / IterationMatrix.r2
    b2 = IterationMatrix.b2
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6))

    # Plot 1
    im1 = ax1.imshow(r1Divr2, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    ax1.set_xlabel(r'$Z_\mathrm{B}$ [-]')
    ax1.set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im1, ax = ax1, label = r'$r_1/r_2$ [-]')

    # Plot 2
    im2 = ax2.imshow(b2 * 1e3, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    ax2.set_xlabel(r'$Z_\mathrm{B}$ [-]')
    ax2.set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im2, ax = ax2, label = r'$b_2$ [mm]')

    # Set a common title for the figure
    fig.suptitle('Compressor stage dimensions\n' + r'$r_\mathrm{h} = $' + str(round(Compressor.rh, 3) * 1000) + ' mm \n' + r'$r_\mathrm{h}/r_1$ = ' + str(round(Compressor.rh / Compressor.r1, 3)))

    # Adjust the spacing between subplots
    fig.subplots_adjust(top = 0.82, wspace = 0.4)
    #plt.show()


def plot_temperatures(Compressor, InletConditions, IterationMatrix):
    """ Plot temperatures of compressor stage after geometry optimization """
    
    ZB = IterationMatrix.ZB
    beta2B = np.rad2deg(IterationMatrix.beta2B)    
    T02m = IterationMatrix.T02m - 273.15
    T2m = IterationMatrix.T2m - 273.15
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6))

    # Plot 1
    im1 = ax1.imshow(T02m, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    ax1.set_xlabel(r'$Z_\mathrm{B}$ [-]')
    ax1.set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im1, ax = ax1, label = r'$T_\mathrm{02m}\,\, [^\circ\mathrm{C}]$')

    # Plot 2
    im2 = ax2.imshow(T2m, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    ax2.set_xlabel(r'$Z_\mathrm{B}$ [-]')
    ax2.set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im2, ax = ax2, label = r'$T_\mathrm{2m}\,\, [^\circ\mathrm{C}]$')

    # Set a common title for the figure
    fig.suptitle('Compressor stage temperatures\n' + r'$T_{00} = $' + str(round(InletConditions.T00 - 273.15, 3)) + '\n' + r'$T_1 = $' + str(round(Compressor.T1 - 273.15, 3)))

    # Adjust the spacing between subplots
    fig.subplots_adjust(top = 0.82, wspace = 0.4)
    #plt.show()


def plot_efficiency_shaft_work(Compressor, InletConditions, IterationMatrix):
    """ Plot stage compressor efficiency and total shaft work after geometry optimization """
    
    ZB = IterationMatrix.ZB
    beta2B = np.rad2deg(IterationMatrix.beta2B)   
    etaStage = IterationMatrix.eta
    Wx = IterationMatrix.Wx
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6))

    # Plot 1
    im1 = ax1.imshow(etaStage, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    ax1.set_xlabel(r'$Z_\mathrm{B}$ [-]')
    ax1.set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im1, ax = ax1, label = r'Required $\eta_\mathrm{stage}$ [-]')

    # Plot 2
    im2 = ax2.imshow(Wx / 1e3, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    ax2.set_xlabel(r'$Z_\mathrm{B}$ [-]')
    ax2.set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im2, ax = ax2, label = r'$W_\mathrm{x}$ [kJ/kg]')

    # Set a common title for the figure
    fig.suptitle('Compressor stage efficiency and shaft work \n' + r'$\eta_\mathrm{rotor} = $' + str(round(Compressor.eta_rotor, 3)) + '\n' + r'$\Delta h_\mathrm{0s} = $' + str(round(Compressor.dh0s / 1e3, 1)) + ' kJ/kg')

    # Adjust the spacing between subplots
    fig.subplots_adjust(top = 0.82, wspace = 0.4)
    #plt.show()


def plot_pressure(Compressor, InletConditions, IterationMatrix):
    """ Plot stage compressor pressures after geometry optimization """
    
    ZB = IterationMatrix.ZB
    beta2B = np.rad2deg(IterationMatrix.beta2B)  
    P02m = IterationMatrix.P02m
    P2m = IterationMatrix.P2m
    p5 = IterationMatrix.p5
    Pr = IterationMatrix.Pr
    
    # Create a figure with two subplots
    fig, axs = plt.subplots(2, 2, figsize = (10, 10))

    # Plot 1
    im1 = axs[0, 0].imshow(Pr, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()], vmin = 1)
    axs[0, 0].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[0, 0].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im1, ax = axs[0, 0], label = r'Pr [-]')

    # Plot 2
    im2 = axs[0, 1].imshow(p5 / 1e5, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[0, 1].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[0, 1].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im2, ax = axs[0, 1], label = r'$p_5$ [bar]')

    # Plot 3
    im3 = axs[1, 0].imshow(P02m / 1e5, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[1, 0].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[1, 0].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im3, ax = axs[1, 0], label = r'$p_\mathrm{02m}$ [bar]')

    # Plot 4
    im4 = axs[1, 1].imshow(P2m / 1e5, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[1, 1].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[1, 1].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im4, ax = axs[1, 1], label = r'$p_\mathrm{2m}$ [bar]')

    # Set a common title for the figure
    fig.suptitle('Compressor stage pressures \n' + r'$P_{00} = $' + str(round(InletConditions.P00 / 1e5, 3)) + ' bar')

    # Adjust the spacing between subplots
    fig.subplots_adjust(top = 0.9, wspace = 0.4)
    #plt.show()


def plot_velocity(Compressor, InletConditions, IterationMatrix):
    """ Plot stage compressor velocities after geometry optimization """
    
    ZB = IterationMatrix.ZB
    beta2B = np.rad2deg(IterationMatrix.beta2B)  
    U2 = IterationMatrix.U2
    M2 = IterationMatrix.M2
    C2 = IterationMatrix.c2
    Cm2m = IterationMatrix.cm2m
    
    # Create a figure with two subplots
    fig, axs = plt.subplots(2, 2, figsize = (10, 10))

    # Plot 1
    im1 = axs[0, 0].imshow(U2, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[0, 0].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[0, 0].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im1, ax = axs[0, 0], label = r'$U_2$ [m/s]')

    # Plot 2
    im2 = axs[0, 1].imshow(M2, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[0, 1].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[0, 1].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im2, ax = axs[0, 1], label = r'$M_2$ [-]')

    # Plot 3
    im3 = axs[1, 0].imshow(Cm2m, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[1, 0].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[1, 0].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im3, ax = axs[1, 0], label = r'$C_\mathrm{m2m}$ [m/s]')

    # Plot 4
    im4 = axs[1, 1].imshow(C2, aspect = 'auto', extent = [ZB.min(), ZB.max(), beta2B.max(), beta2B.min()])
    axs[1, 1].set_xlabel(r'$Z_\mathrm{B}$ [-]')
    axs[1, 1].set_ylabel(r'$\beta_\mathrm{2B}$ [deg]')
    fig.colorbar(im4, ax = axs[1, 1], label = r'$C_2$ [m/s]')

    # Set a common title for the figure
    fig.suptitle('Compressor stage velocities \n' + r'$C_\mathrm{m1} = $' + str(round(Compressor.Cm1, 3)) + ' m/s \n' + r'$U_\mathrm{1t} = $' + str(round(Compressor.U1t, 3)) + ' m/s')

    # Adjust the spacing between subplots
    fig.subplots_adjust(top = 0.9, wspace = 0.4)
    #plt.show()


def plot_text(Compressor):
    """ Plot a summary of some important design point conditions """
    
    # Text for summarizing design point conditions
    text1 = r"$\bf{Design\,\, point:}$" + "\n" \
            "Pressure ratio = " + str(Compressor.Pr) +  "\n" \
            "Rotational speed = " + str(round(Compressor.Ndes, 1)) + " rpm \n" \
            "Blade number = " + str(Compressor.bladeNumber) + "\n" \
            "Blade angle = " + str(np.rad2deg(Compressor.bladeAngle)) + r"$^\circ$" + "\n" \
            "Stage efficiency = " + str(round(Compressor.etaStage, 3)) + "\n"
    text2 = "\n" + r"$\bf{Velocities:}$" + "\n" \
            r"$C_\mathrm{m1} = $" + str(round(Compressor.Cm1, 3)) + " m/s \n" \
            r"$U_\mathrm{1t} = $" + str(round(Compressor.U1t, 3)) + " m/s \n" \
            r"$W_\mathrm{1t} = $" + str(round(Compressor.W1t, 3)) + " m/s \n" \
            r"$M_1 = $" + str(round(Compressor.M1, 3)) + "\n" \
            r"$U_2 = $" + str(round(Compressor.U2, 3)) + " m/s \n" \
            r"$U_\mathrm{2t, crit} = $" + str(round(Compressor.U2Crit, 3)) + " m/s \n"
    text3 = "\n\n" + r"$\bf{Compressor\,\, dimensions:}$" + "\n" \
            r"$r_\mathrm{h}$ = " + str(round(Compressor.rh * 1000, 3)) + " mm \n" \
            r"$r_\mathrm{t1}$ = " + str(round(Compressor.r1 * 1000, 3)) + " mm \n" \
            r"$r_\mathrm{t2}$ = " + str(round(Compressor.r2 * 1000, 3)) + " mm \n" \
            r"$b_2$ = " + str(round(Compressor.b2 * 1000, 3)) + " mm \n" \
            r"$\frac{r_\mathrm{h}}{r_1}$ = " + str(round(Compressor.rh / Compressor.r1, 3)) + "\n" \
            r"$\frac{r_1}{r_2}$ = " + str(round(Compressor.r1 / Compressor.r2, 3)) + "\n" 
    text = [text1, text2, text3]     
    text1 = text[0]
    text2 = text[1]
    text3 = text[2]
   
    fig, axs21 = plt.subplots(1, 3)
    fig.set_figwidth(15)
    fig.set_figheight(15)
    fig.tight_layout(pad = 7.0)
    fig.subplots_adjust(top = 0.9, bottom = 0.09)
    fig.canvas.manager.set_window_title('Design point conditions from geometry optimization')

    # Textbox 1
    j = 0
    axs21[j].axis('off')  # Turn off the axis
    axs21[j].text(0.0, 0.5, text1, ha = 'left', va = 'center', fontsize = 12, linespacing = 1.8)
    
    # Textbox 2
    j = 1
    axs21[j].axis('off')  # Turn off the axis
    axs21[j].text(0.0, 0.5, text2, ha = 'left', va = 'center', fontsize = 12, linespacing = 1.8)

    # Textbox 3
    j = 2           
    axs21[j].axis('off')  # Turn off the axis
    axs21[j].text(0.0, 0.5, text3, ha = 'left', va = 'center', fontsize = 12, linespacing = 1.8)
    #plt.show()


def plot_min_relative_velocity(Compressor, InletConditions):
    """ Plot relative velocity as a function of inlet meridional velocity """

    fig, ax11 = plt.subplots()
    plt.plot(InletConditions.Cm1i, Compressor.W1ti)
    plt.plot(Compressor.Cm1, Compressor.W1t, 'ro')
    ax11.tick_params(axis = 'both', which = 'major', labelsize = 10)
    plt.xlabel(r'$C_{\mathrm{m1}}$ [m/s]', fontsize = 12)
    plt.ylabel(r'$W_{\mathrm{1t}}$ [m/s]', fontsize = 12)
    plt.title(r'Minimization of relative velocity')
    plt.grid()
    #plt.show()


def plot_geometry(Compressor, InletConditions, IterationMatrix):
    plot_min_relative_velocity(Compressor, InletConditions)
    plot_text(Compressor)
    plot_dimensions(Compressor, InletConditions, IterationMatrix)
    plot_temperatures(Compressor, InletConditions, IterationMatrix)
    plot_efficiency_shaft_work(Compressor, InletConditions, IterationMatrix)
    plot_pressure(Compressor, InletConditions, IterationMatrix)
    plot_velocity(Compressor, InletConditions, IterationMatrix)
    #plt.show()


if __name__ == '__main__':
    pass