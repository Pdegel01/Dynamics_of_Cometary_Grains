import matplotlib.pyplot as plt       # Importing the plotting module
import numpy as np                    # Importing the numerical computation module                       
import scipy.integrate as integrate   # Importing the numerical integration module
import sys                            # Importing the system module
import importlib

# Importing constants, inputs, and source code
import inputs
import constants
import code_source

from constants import *
from inputs import *
from code_source import *

# Initialize conditions for the first grain
first_grain_initial_conditions = initial_conditions[0]
position = first_grain_initial_conditions['position']
velocity = first_grain_initial_conditions['velocity']

# Calculate altitude in kilometers, useful for the legend of the plot
altitude = np.linalg.norm(position) - earth_radius
altitude_km = int(altitude * 10**-3)

# List of grain radii and masses to test
grain_radii = np.logspace(-9, 1, 100)  
grain_masses = np.logspace(-22,-5,100)  

radius_grain = grain_radius
mass_grain = grain_mass

# Initialize lists to store the magnitudes of different forces
magnitude_earth_gravity = []
magnitude_solar_gravity = []
magnitude_solar_pressure = []
magnitude_gas_drag = []
magnitude_poynting = []
magnitude_lorentz = []
magnitude_comet_gravity = []
magnitude_radial_tidal = []

# Loop over different grain masses
for grain_mass in grain_masses:   
    
    if environment == "comet":
        # Calculate comet gravity force
        comet_force = comet_gravity_force(comet_mass, grain_mass) * (position) / np.linalg.norm(position)**3
        magnitude_comet_gravity.append(np.linalg.norm(comet_force))
        
        # Calculate radial tidal force
        radial_tidal_force = radial_tidal_force_magnitude(coeff_disable_tidal, solar_mass, grain_mass) * (position / np.linalg.norm(position)) / np.linalg.norm(solar_position - position)**3
        magnitude_radial_tidal.append(np.linalg.norm(radial_tidal_force))
        
    else:
        # Calculate Earth gravity force
        gravitational_force = earth_gravity_force(grain_mass, position)
        magnitude_earth_gravity.append(np.linalg.norm(gravitational_force))
    
    # Calculate solar gravity force
    solar_gravity_force = solar_gravity_magnitude(coeff_disable_solar, grain_mass, solar_mass) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3 - (solar_position / np.linalg.norm(solar_position)**3))
    magnitude_solar_gravity.append(np.linalg.norm(solar_gravity_force))
    
    # Calculate solar radiation pressure force
    solar_radiation_pressure_force = solar_radiation_pressure_force_magnitude(Q_pr, grain_radius) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3) 
    magnitude_solar_pressure.append(np.linalg.norm(solar_radiation_pressure_force))
    
    # Calculate Poyting-Robertson drag force
    poynting_robertson_drag_force = (poynting_robertson_drag_magnitude(coeff_disable_poynting, grain_radius)) * (velocity / np.linalg.norm(velocity)) * (1 / np.sqrt(np.linalg.norm(solar_position - position)**5)) 
    magnitude_poynting.append(np.linalg.norm(poynting_robertson_drag_force))
    
    # Calculate gas drag force
    gas_drag_force = gas_drag_force_magnitude(drag_coeff, gas_mol_mass, gas_number_density, grain_radius) * np.linalg.norm(abs(velocity - gas_velocity)) * abs(velocity - gas_velocity)
    magnitude_gas_drag.append(np.linalg.norm(gas_drag_force))

# Convert lists to numpy arrays for easier plotting
earth_gravity_list = np.array(magnitude_earth_gravity)
comet_gravity_list = np.array(magnitude_comet_gravity)
radial_tidal_list = np.array(magnitude_radial_tidal) 
solar_gravity_list = np.array(magnitude_solar_gravity)
solar_radiation_pressure_list = np.array(magnitude_solar_pressure)
gas_drag_list = np.array(magnitude_gas_drag)
poynting_robertson_list = np.array(magnitude_poynting)

''' FORCES AS A FUNCTION OF THE RADIUS '''

# Initialize lists to store the magnitudes of different forces as a function of radius
magnitude_earth_gravity_radius = []
magnitude_solar_gravity_radius = []
magnitude_solar_pressure_radius = []
magnitude_gas_drag_radius = []
magnitude_poynting_radius = []
magnitude_lorentz_radius = []
magnitude_comet_gravity_radius = []
magnitude_radial_tidal_radius = []

# Loop over different grain radii
for grain_radius in grain_radii:   
    
    if environment == "comet":
        # Calculate comet gravity force
        comet_force = comet_gravity_force(comet_mass, grain_mass) * (position) / np.linalg.norm(position)**3
        magnitude_comet_gravity_radius.append(np.linalg.norm(comet_force))
        
        # Calculate radial tidal force
        radial_tidal_force = radial_tidal_force_magnitude(coeff_disable_tidal, solar_mass, grain_mass) * (position / np.linalg.norm(position)) / np.linalg.norm(solar_position - position)**3
        magnitude_radial_tidal_radius.append(np.linalg.norm(radial_tidal_force))
        
    else:
        # Calculate Earth gravity force
        gravitational_force = earth_gravity_force(grain_mass, position)
        magnitude_earth_gravity_radius.append(np.linalg.norm(gravitational_force))
    
    # Calculate solar gravity force
    solar_gravity_force = solar_gravity_magnitude(coeff_disable_solar, grain_mass, solar_mass) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3 - (solar_position / np.linalg.norm(solar_position)**3))
    magnitude_solar_gravity_radius.append(np.linalg.norm(solar_gravity_force))
    
    # Calculate solar radiation pressure force
    solar_radiation_pressure_force = solar_radiation_pressure_force_magnitude(Q_pr, grain_radius) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3) 
    magnitude_solar_pressure_radius.append(np.linalg.norm(solar_radiation_pressure_force))
    
    # Calculate Poyting-Robertson drag force
    poynting_robertson_drag_force = (poynting_robertson_drag_magnitude(coeff_disable_poynting, grain_radius)) * (velocity / np.linalg.norm(velocity)) * (1 / np.sqrt(np.linalg.norm(solar_position - position)**5)) 
    magnitude_poynting_radius.append(np.linalg.norm(poynting_robertson_drag_force))
    
    # Calculate gas drag force
    gas_drag_force = gas_drag_force_magnitude(drag_coeff, gas_mol_mass, gas_number_density, grain_radius) * np.linalg.norm(abs(velocity - gas_velocity)) * abs(velocity - gas_velocity)
    magnitude_gas_drag_radius.append(np.linalg.norm(gas_drag_force))

# Convert lists to numpy arrays for easier plotting
earth_gravity_list_radius = np.array(magnitude_earth_gravity_radius)
comet_gravity_list_radius = np.array(magnitude_comet_gravity_radius)
radial_tidal_list_radius = np.array(magnitude_radial_tidal_radius) 
solar_gravity_list_radius = np.array(magnitude_solar_gravity_radius)
solar_radiation_pressure_list_radius = np.array(magnitude_solar_pressure_radius)
gas_drag_list_radius = np.array(magnitude_gas_drag_radius)
poynting_robertson_list_radius = np.array(magnitude_poynting_radius)

# Plotting the results
if environment == "debris":
    
    plt.figure(figsize=(10, 8))
    plt.loglog(grain_masses, solar_gravity_list, label='Solar Gravity')
    plt.loglog(grain_masses, solar_radiation_pressure_list, label='Radiation Pressure')
    plt.loglog(grain_masses, gas_drag_list, label='Atmospheric Drag')
    plt.loglog(grain_masses, poynting_robertson_list, label='Poynting-Robertson')
    plt.loglog(grain_masses, earth_gravity_list, label='Earth Gravity')

    # Set log scale for x-axis
    plt.xscale('log')
    plt.xlabel('Debris Mass (kg)', fontsize=16)
    plt.ylabel('Magnitude', fontsize=16)
    plt.title(f"Force Ratios by mass, radius = {radius_grain} m, altitude = {altitude_km}km", fontsize=20)
    plt.legend(fontsize=16)
    plt.grid(True)

    # Increase tick size on axes
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
    
    # Plotting forces as a function of grain radius
    plt.figure(figsize=(10, 8))
    plt.loglog(grain_radii, solar_gravity_list_radius, label='Solar Gravity')
    plt.loglog(grain_radii, solar_radiation_pressure_list_radius, label='Radiation Pressure')
    plt.loglog(grain_radii, gas_drag_list_radius, label='Atmospheric Drag')
    plt.loglog(grain_radii, poynting_robertson_list_radius, label='Poynting-Robertson')
    plt.loglog(grain_radii, earth_gravity_list_radius, label='Earth Gravity')

    # Set log scale for x-axis
    plt.xscale('log')
    plt.xlabel('Debris Radius (m)', fontsize=16)
    plt.ylabel('Magnitude', fontsize=16)
    plt.title(f"Force Ratios by radius, mass = {mass_grain}kg, altitude = {altitude_km}km", fontsize=20)
    plt.legend(fontsize=16)
    plt.grid(True)

    # Increase tick size on axes
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()

elif environment == "comet":
    
    plt.figure(figsize=(10, 8))
    plt.loglog(grain_masses, solar_gravity_list, label='Solar Gravity')
    plt.loglog(grain_masses, solar_radiation_pressure_list, label='Radiation Pressure')
    plt.loglog(grain_masses, gas_drag_list, label='Gas Drag')
    plt.loglog(grain_masses, poynting_robertson_list, label='Poynting-Robertson')
    plt.loglog(grain_masses, comet_gravity_list, label='Comet Gravity')
    plt.loglog(grain_masses, radial_tidal_list, label='Radial Tidal')

    # Set log scale for x-axis
    plt.xscale('log')
    plt.xlabel('Debris Mass (kg)', fontsize=16)
    plt.ylabel('Magnitude', fontsize=16)
    plt.title(f"Force Ratios by Mass, radius = {radius_grain}m", fontsize=20)
    plt.legend(fontsize=16)
    plt.grid(True)

    # Increase tick size on axes
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
    
    # Plotting forces as a function of grain radius
    plt.figure(figsize=(10, 8))
    plt.loglog(grain_radii, solar_gravity_list_radius, label='Solar Gravity')
    plt.loglog(grain_radii, solar_radiation_pressure_list_radius, label='Radiation Pressure')
    plt.loglog(grain_radii, gas_drag_list_radius, label='Gas Drag')
    plt.loglog(grain_radii, poynting_robertson_list_radius, label='Poynting-Robertson')
    plt.loglog(grain_radii, comet_gravity_list_radius, label='Comet Gravity')
    plt.loglog(grain_radii, radial_tidal_list_radius, label='Radial Tidal')

    # Set log scale for x-axis
    plt.xscale('log')
    plt.xlabel('Debris Radius (m)', fontsize=16)
    plt.ylabel('Magnitude', fontsize=16)
    plt.title(f"Force Ratios by Radius, mass = {mass_grain}kg", fontsize=20)
    plt.legend(fontsize=16)
    plt.grid(True)

    # Increase tick size on axes
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
