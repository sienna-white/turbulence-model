# turbulence-model

This directory contains a simple 1D vertical water column model that simulates turbulence, temperature, and phytoplankton dynamics in response to wind forcing, heat flux, and light availability. The model is implemented in Python and uses numerical methods to solve the governing equations.

The main code you will use is `watercolumn-model.py`, which sets up the model parameters, initializes the state variables, and runs the simulation loop. The model relies on several helper functions and classes defined in the `model_code/` directory. 

Summary of files in model_code/ 

    advance_variables.py: 
        Contains functions to advance the velocity and scalar fields in time using finite difference methods.
        Includes the following functions.
            1. TDMA: Tri Diagonal Matrix Algorithm solver for tridiagonal systems.
            2. advance_velocity: Advances the velocity profile using the TDMA solver.
            3. advance_algae: Advances the phytoplankton concentration profile.
            4. advance_scalar: Advances scalar fields like temperature or nutrients.
            5. advance_Q2: Advances the turbulent kinetic energy profile.
            6. advance_Q2L: Advances the product of turbulent kinetic energy and length scale

    calculate_physical_variables.py:
        Contains functions to calculate physical variables such as density, Brunt-Vaisala frequency, the turbulent length scale, and turbulent diffusivities. Includes the following functions:
            1. calculate_rho: Calculates density based on temperature.
            2. calculate_Gh: Calculates the stability function Gh (MY scheme)
            3. calculate_ustar: Calculates friction velocity from wind speed.
            4. calculate_turbulent_functions: Calculates turbulent diffusivities and viscosities.
            5. calculate_sm: Calculates the stability function Sm.
            6. calculate_sh: Calculates the stability function Sh.
            7. calculate_brunt_vaisala: Calculates the Brunt-Vaisala frequency (N^2).
            8. calculate_photic_depth: Determines the photic depth based on light intensity.
            9. add_noise_floor: Adds a noise floor to avoid negative values.
           10. calculate_lengthscale: Calculates the turbulent length scale.
    

    constants.py
        contains the constants used in the model such as densities, specific heat capacities, drag coefficients, and parameters for the MY2.5 scheme.
    
    forcings.py
        Contains functions to generate and read external forcings such as wind speed, heat flux, pressure gradients, and light availability. Includes the following functions:
            1. diurnal_light: Generates a diurnal light cycle. I made these equations up!!! 
            2. get_pressure_at_timestep: Calculates pressure gradient forcing at a given timestep.
            3. wind_speed: Generates a diurnal wind speed profile using a cosine function.
            4. temperature: Generates a diurnal temperature profile using a cosine function.
            5. analytical_temperature_profile: Provides an analytical temperature profile based on a hyperbolic tangent function.
            6. read_forcings_from_file: Reads forcings from an external CSV file.
            7. air_temperature: Generates air temperature based on diurnal cycle.
            8. [deprecated] wind: Calculates wind stress based on wind speed.
            9. [deprecated] generate_forcings: Reads external forcings if provided as a csv.

    initial_conditions.py
        Contains functions to initialize the model state variables such as temperature profile, turbulent kinetic energy, length scale, and Brunt-Vaisala frequency. Includes the following functions:
            1. temp_profile: Initializes the temperature profile based on stratification option.
            2. initialize_arrays: Initializes arrays for turbulent kinetic energy, length scale, and diffusivities.
            3. initialize_N_BV: Initializes the Brunt-Vaisala frequency based on density profile.
            4. initialize_turbulent_functions: Initializes turbulent functions based on Brunt-Vaisala frequency.
            5. check_initial_condition: Checks for existing initial conditions from a CSV file.

    phytoplankton.py
        Phytoplankton are coded as individual Algae_Species "classes" with their own growth, loss, and light shading parameters. This file contains the Algae_Species class and functions to calculate light attenuation due to self-shading and background turbidity. Includes the following classes and functions:
            1. Algae_Species: Class definition representing a phytoplankton species with methods for growth, loss, and light calculation.
            2. SelfShading: Class to calculate light attenuation due to self-shading.
            3. self_shading: Function to calculate light intensity at each depth using Lambert-Beer's Law.
    
    save_output.py
        Functions to save output as designated time steps. Contains:
            1. save_run_info: Saves metadata about the model run.
            2. save_1d_data: Saves 1D data that varies only with time.
            3. save_2d_data: Saves 2D data that varies with depth and time.
            4. save_dataset: Saves the dataset to a NetCDF file.

    watercolumn_lib.py
        Main library file that combines all the functions above into one accessible place. Most likely, you will not need to edit this file.



Old notes:
Because the indexing is a little confusing in Python vs. Matlab (0 is the bed, N-1 is 
the top of the water column), and then the point below the top is N-2, when indexing 
the top of the water column, I defined a variable called "top" to be N-1. This is just 
to make the code a little more readable, and hopefully less confusing. Hopefully the concept
of U[0] = bed velocity is a little more intuitive.

There's also a separate python file called watercolumn_lib.py that contains some functions
as well as a class definition called "SavedProfiles". This class is used to store the
profiles at each requested time step, and then plot them at the end. Mostly this was easier
than passing a bunch of arrays around between functions.