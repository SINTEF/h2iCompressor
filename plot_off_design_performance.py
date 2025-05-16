"""
Python script for plotting the off-design performance computed in off_design_performance.py

Author(s): Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2025)
"""


# Import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_off_design(results_off_design, Compressor, InletConditions):
    """ Plot the results from the off-design performance calculations """
    
    # Plot pressure ratio
    fig1 = plt.figure('Total-to-static pressure ratio')
    
    plot_experimental_values = True
    if plot_experimental_values:        
        df = pd.read_csv('properties/lut_compressor_pressure_ratios.csv')
        x_values_1 = df.iloc[:, 0]
        y_values_curve_1 = df.iloc[:, 1]
        x_values_2 = df.iloc[:, 2]
        y_values_curve_2 = df.iloc[:, 3]
        x_values_3 = df.iloc[:, 4]
        y_values_curve_3 = df.iloc[:, 5]
        x_values_4 = df.iloc[:, 6]
        y_values_curve_4 = df.iloc[:, 7]
        x_values_5 = df.iloc[:, 8]
        y_values_curve_5 = df.iloc[:, 9]

    for result in results_off_design:
        Nplot = f"{result['N'] * 1e-3:.0f} krpm"
        plt.plot(result['mdoto'], result['Pro'], label = Nplot)
    plt.ylabel(r'$\frac{P_{3}}{P_{00}}$', rotation = 45, fontsize = 12)
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)      
    plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    if plot_experimental_values:
        plt.plot(np.array(x_values_1), np.array(y_values_curve_1), 'k-', label = 'Design curve 1')
        plt.plot(np.array(x_values_2), np.array(y_values_curve_2), 'k--', label = 'Design curve 2')
        plt.plot(np.array(x_values_3), np.array(y_values_curve_3), 'k-.', label = 'Design curve 3')
        plt.plot(np.array(x_values_4), np.array(y_values_curve_4), 'k:', label = 'Design curve 4')
        plt.plot(np.array(x_values_5), np.array(y_values_curve_5), 'k:', label = 'Design curve 5')
    plt.tight_layout(rect = [0, 0, 0.98, 1])
    #plt.plot([0, 50], [1, 1], 'r--')
    plt.title(r'Total-to-static pressure ratio')
    #plt.grid(True)
    plt.ylim(1, )
    plt.plot(InletConditions.mdot, Compressor.Pr, 'ro', label = 'Design pt.')

    # Plot Mach number
    plt.figure('Impeller exit Mach number')
    for result in results_off_design:
        Nplot = f"{result['N'] * 1e-3:.0f} krpm"
        plt.plot(result['mdoto'], result['M2o'], label = Nplot)
    plt.axvline(x = InletConditions.mdot, color = 'r', linewidth = 0.75, label = r'$\dot{m}_\mathrm{des}$')
    plt.ylabel(r'${Ma}_2$', fontsize = 12)
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)
    plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.title('Impeller Mach number')
    plt.tight_layout(rect = [0, 0, 0.98, 1])
    plt.grid(True)

    # Plot efficiency curve
    plt.figure('Efficiency')
    for result in results_off_design:
        Nplot = f"{result['N'] * 1e-3:.0f} krpm"
        plt.plot(result['mdoto'], result['etao'], label = Nplot)
    plt.axvline(x = InletConditions.mdot, color = 'r', linewidth = 0.75, label = r'$\dot{m}_\mathrm{des}$')
    plt.ylabel(r'$\eta$', rotation = 45, fontsize = 12)
    plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)
    plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.tight_layout(rect = [0, 0, 0.98, 1])
    plt.title(r'Efficiency ')
    plt.grid(True)

    # Plot exit tip velocity
    plt.figure('Exit tip velocity')
    for result in results_off_design:
        plt.plot(result['N'], result['U2o'], 'ko')
    plt.xlabel(r'$N$ [rpm]', fontsize = 12)
    plt.ylabel(r'$U_2 $ [m/s]', fontsize = 12)
    plt.title('Impeller exit tip speed')
    plt.grid(True)

    # Plot enthalpy rise/loss
    for result in results_off_design:
        total_loss = result['enthalpy_rise_loss'][1] + result['enthalpy_rise_loss'][2] + result['enthalpy_rise_loss'][3] + result['enthalpy_rise_loss'][4] + result['enthalpy_rise_loss'][5] + result['enthalpy_rise_loss'][6]      #MSG: Does not match the definition of efficiency where RC and DF are placed in the denominator
        
        plt.figure(r'Enthalpy rise/loss N = ' + f"{result['N'] * 1e-3:.0f} krpm")
        #plt.plot(result['mdoto'], result['enthalpy_rise_loss'][0], label = 'Aerodynamic rise')
        plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][1] / total_loss), label = 'Incidence loss')
        plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][2] / total_loss), label = 'Blade loading loss')
        plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][3] / total_loss), label = 'Skin friction loss')
        plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][4] / total_loss), label = 'Vaneless diffuser loss')
        plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][5] / total_loss), label = 'Recirculation loss')
        plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][6] / total_loss), label = 'Impeller disk friction loss')
        #plt.plot(result['mdoto'], np.array(result['enthalpy_rise_loss'][0] ), label = 'Enthalpy rise')
        plt.plot(result['mdoto'], total_loss / total_loss, label = 'Total loss')
        plt.axvline(x = InletConditions.mdot, color = 'r', linewidth = 0.75, label = r'$\dot{m}_\mathrm{des}$')
        plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 12)
        plt.ylabel(r'Share of total loss [%]', fontsize = 12)
        plt.title(r'Enthalpy rise/loss N = ' + f"{result['N'] * 1e-3:.0f} krpm")
        plt.legend()
        plt.ylim(0,)
        plt.xlim(min(result['mdoto']), max(result['mdoto']))

    # Plot performance map
    plt.figure('Performance map')
    
    # Obtain the lines of constant efficiency for plotting the compressor performance map
    constant_eff_lines = Compressor.constant_eff_lines
    constant_eff_line_mdot = [[[[], []] for i in range(len(constant_eff_lines))] for _ in range(len(Compressor.N_off_design_arr))]
    constant_eff_line_pr = [[[[], []] for i in range(len(constant_eff_lines))] for _ in range(len(Compressor.N_off_design_arr))]
    
    for iN in range(0, len(Compressor.N_off_design_arr)):
        mdoto = results_off_design[iN]['mdoto']
        Pro = results_off_design[iN]['Pro']
        etao = results_off_design[iN]['etao']
        
        # Finding the constant efficiency lines
        for i in range(len(constant_eff_lines)):
            if any(etao >= constant_eff_lines[i]):
                if etao[0] <= constant_eff_lines[i]:
                    index_0 = np.argmax(etao >= constant_eff_lines[i])
                    constant_eff_line_mdot[iN][i][0] = mdoto[index_0]
                    constant_eff_line_pr[iN][i][0] = Pro[index_0]
                    if etao[- 1] <= constant_eff_lines[i]:
                        index_1 = np.argmax(etao[index_0:] < constant_eff_lines[i])
                        constant_eff_line_mdot[iN][i][1] = mdoto[index_0 + index_1]
                        constant_eff_line_pr[iN][i][1] = Pro[index_0 + index_1]
                    else:
                        constant_eff_line_mdot[iN][i][1] = np.nan
                        constant_eff_line_pr[iN][i][1] = np.nan
                else:
                    if etao[- 1] <= constant_eff_lines[i]:
                        index_1 = np.argmax(etao < constant_eff_lines[i])
                        constant_eff_line_mdot[iN][i][1] = mdoto[index_1]
                        constant_eff_line_pr[iN][i][1] = Pro[index_1]
                        constant_eff_line_mdot[iN][i][0] = np.nan
                        constant_eff_line_pr[iN][i][0] = np.nan
                    else:
                        constant_eff_line_mdot[iN][i][0] = np.nan
                        constant_eff_line_pr[iN][i][0] = np.nan
                        constant_eff_line_mdot[iN][i][1] = np.nan
                        constant_eff_line_pr[iN][i][1] = np.nan
            else:
                constant_eff_line_mdot[iN][i][0] = np.nan
                constant_eff_line_pr[iN][i][0] = np.nan
                constant_eff_line_mdot[iN][i][1] = np.nan
                constant_eff_line_pr[iN][i][1] = np.nan    
    
    colors = ['r', 'g', 'b', 'c', 'm', 'y']  # Colors for the constant efficiency lines    
    
    for result in results_off_design:
        Nplot = f"RPM = {result['N']}"
        plt.plot(result['mdoto'], result['Pro'], label = Nplot)
    
    for i in range(len(constant_eff_lines)):
        line_1_mdot = [constant_eff_line_mdot[j][i][0] for j in range(len(constant_eff_line_mdot))]
        line_2_mdot = [constant_eff_line_mdot[j][i][1] for j in range(len(constant_eff_line_mdot))]
        line_1_pr = [constant_eff_line_pr[j][i][0] for j in range(len(constant_eff_line_pr))]
        line_2_pr = [constant_eff_line_pr[j][i][1] for j in range(len(constant_eff_line_pr))]
        plt.plot(line_1_mdot, line_1_pr, linestyle = '--', color = colors[i], label = f'Efficiency = {constant_eff_lines[i]}')
        plt.plot(line_2_mdot, line_2_pr, linestyle = '--', color = colors[i])
    plt.plot(InletConditions.mdot, Compressor.Pr, 'ro', label = 'Design point')  
    plt.xlabel(r'$\dot{m}$' + ' [kg/s]', fontsize = 12)      
    plt.ylabel(r'$\frac{P_{3}}{P_{00}}$', rotation = 45, fontsize = 12)
    plt.title('Performance map')
    plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.tight_layout(rect = [0, 0, 0.98, 1])
    plt.ylim(1, )
    plt.grid(True)

    # Plot exit velocity V2
    plt.figure('Exit velocity V2')
    for result in results_off_design:
        plt.plot(result['mdoto'], result['V2'], label = f"{result['N'] * 1e-3:.0f} krpm")
    plt.axvline(x = InletConditions.mdot, color = 'r', linewidth = 0.75, label = r'$\dot{m}_\mathrm{des}$')
    plt.xlabel(r'$\dot{m}$' + ' [kg/s]', fontsize = 12)
    plt.ylabel(r'$V_2$' + ' [m/s]', fontsize = 12)
    plt.title('Exit velocity V2')   
    plt.legend()

    #plt.show()


if __name__ == '__main__':
    pass