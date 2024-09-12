"""
The following Python code iterates through a range of rotational speeds to deteremine compressor off-design performance.
A geometrical basis is found in geometry.py which is further built on in off-design_performance.py. Plotting of several variables of interest 
is done by utilizing plot_compessor.py, plot_system.py, plot_text.py, plot_velocities.py and pressure_test.py.
The script called settings.py is used for all parametrization of for instance fluid properties, inlet flow conditions, diffuser propertis etc. 

Authors: Petter Resell (summer intern, 2024), Martin Spillum Grønli (SINTEF Energy Research, 2024)
"""


def main():
    fluid_name = 'h2'           # Select working fluid, 'h2' or 'air'
    import settings
    import geometry
    import matplotlib       
    matplotlib.use('tkagg')     # Use tkagg to avoid crashing when using X11-forwarding for plotting
    from matplotlib import pyplot as plt    

    # Load and make instances of classes for fluid, inlet conditions, compressor and iteration matrices
    Fluid = settings.Fluid(path_to_fluid_properties_toml = './properties/' + fluid_name + '.toml')
    InletConditions = settings.InletConditions(fluid_instance = Fluid, path_to_inlet_toml = './properties/' + fluid_name + '.toml')
    Compressor = settings.Compressor(fluid_instance = Fluid, inlet_conditions_instance = InletConditions, path_to_compressor_toml = './properties/compressor.toml')
    IterationMatrix = settings.IterationMatrix(compressor_instance = Compressor)
    
    print('\nCalculating geometry...') 
    geometry.inducer_and_impeller_calculations(Fluid, InletConditions, Compressor)
    geometry.iterate_blade_number_and_blade_angle(Compressor, InletConditions, Fluid, IterationMatrix)
    geometry.print_and_plot_geometry(Compressor, InletConditions, IterationMatrix)
    plt.show()

    print('\nCalculating off-design performance...')


if __name__ == '__main__':
    main()



""" Running off-desing_performance.py for rotational speeds of 15000rpm to the critical speed with increments of 5000 """
Narr = np.arange(15000, geometry.Ncrit , 5000)      

# ------------------- Plotting from off_design_performance.py -------------------
for iN in range(0, len(Narr)):                                                                      # Running and plotting for all rpm N's
    Pro, P03o, T2oabs, mdoto, etao, U2o, M2o = odp.off_design_performance(Narr[iN])                      # Running everything
    Pro = np.array(Pro)                                                                             # Pressure ratio
    P03o = np.array(P03o)                                                                           # Outlet pressure
    T2oabs = np.array(T2oabs)                                                                       # Outlet temperature
    mdoto = np.array(mdoto)                                                                         # mass flow
    etao = np.array(etao)                                                                           # Efficiency
    correctedMass = mdoto*(P03o/settings.P00)/np.sqrt(T2oabs/settings.T00)        # Corrected mass flow parameter

    """ checking if all efficiencies are greater than 0 and less or equal to 1 """
    if all(etao > 0) and all(etao <=1):       
        Nplot = Narr[iN]
        Nplot = str(Nplot *10**-3)+ ' Krpm'
        
        # Plotting pressure ratio 
        fig1 = plt.figure('Pressure ratio')
        plt.plot(mdoto, Pro, label=Nplot  ) # marker = markers[iN],
        plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m}$' + ' ' + '[kg/s]', fontsize = 15)
        # plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.plot([0, 50], [1, 1], 'r--')
        plt.title(r'Pressure ratio')
        plt.grid(True)
        

        # Compressor map 
        fig99 = plt.figure('Compressor map')
        plt.plot(correctedMass, Pro, label=Nplot )# marker = markers[iN]
        plt.ylabel(r'$\frac{P_{03}}{P_{00}}$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m} \frac{P_{03}}{P_{00}} \sqrt{\frac{T_{00}}{T_{03}} }  $' + ' ' +'[kg/s]', fontsize = 15)
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title('Compressor map')
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        # plt.xlim([0, 50])
        plt.grid(True)

        # Mach number
        fig19 = plt.figure('Impeller exit mach number')
        plt.plot(mdoto, M2o, label= Nplot )# marker = markers[iN]
        plt.ylabel( r'${Ma}_2$', fontsize = 15)
        plt.xlabel(r'$\dot{m}$' + ' ' +'[kg/s]', fontsize = 15)
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title('Impeller mach number')

        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.grid(True)
        
        # Plotting efficiency curve 
        fig2 = plt.figure('Efficiency')
        plt.plot(mdoto, etao, label= Nplot  ) # marker = markers[iN],
        plt.grid()
        plt.ylabel(r'$\eta$', rotation=45, fontsize = 15)
        plt.xlabel(r'$\dot{m} $'+ ' ' + '[kg/s]', fontsize = 15)
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.title(r'Efficiency ')
        plt.grid(True)
        
        fig4 = plt.figure('Exit Tip velocity')
        plt.plot(Narr[iN],U2o, 'ko', label= Nplot  ) # marker = markers[iN],
        plt.grid()
        plt.xlabel(r'$N $  [rpm]', fontsize = 15)
        plt.ylabel(r'$U_2 $ [m/s]', fontsize = 15)
        plt.title('Impeller exit tip speed')
        plt.grid(True)


        debug = 0
    else:
        print('   -> RPM  = ' + str(Narr[iN]) + ' has efficiency outside range [0, 1] !')

fig1 = plt.figure('Pressure ratio', label='Design pt.')
plt.plot(settings.mdot, settings.Pr, 'ro')


print('off-design_performance.py successfully run. \n')


#plt.show()