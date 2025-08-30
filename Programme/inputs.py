import numpy as np
import sys
import random

from constants import e, me, mi, kb, c, epsilon_0, Ar, G, earth_mass, solar_mass, earth_radius, moon_radius, J_2, sun_luminosity, AU




'''' INPUTS - FOR ELECTRON AND ION CURRENTS '''

Te = 10*11604                         #Electron temperature (K)
Ti = 10*11604                         #Ion temperature (K)
Tph = 2*11604                         #Photoelectron temperature (K)
Tsec = 2.5*11604                      #Secondary electron temperature (K)


ne = 5e6                              #Electrons density in plasma (m^-3)
ni = ne                               #Ions density in plasma (m^-3)

z_number = 1                          #atomic number of target
drag_coeff = 2                        #Drag coefficient
gas_mol_mass = 1e-3                   #Mass of a gas molecule (kg)
gas_number_density = 1e-12            #Gas number density (m^-3)
gas_velocity = [0, 465, 0]            #Gas velocity (m/s)


#Example for the plasma parameters : [Horanyi 1996] [5] takes Te = Ti = 10eV ; Tsec = 2.5 eV ;  ne = ni = 5 cm^-3 ; grain_radius = 1^-6 m ; gamma = 0.1 ; delta_max = 1 ; E_max = 300 eV




''' INPUTS - CHARACTERISTICS OF THE GRAIN - AND THE COMET (IF COMETOCENTRIC ENVIRONMENT) '''

initial_potential = 0                 #Surface potential of the spherical grain (V)
grain_radius = 1e-6                   #Radius of the sphere (m)
grain_temperature = 200               #Temperature at the surface of the grain (K)
grain_mass = 1e-12                    #Mass of the grain (kg)
 
                       

comet_mass = 1e12                                                   #Mass of the comet (kg)
comet_radius = 1000                                                 # Radius of the comet, in meters
comet_position = np.array([0.0, 100.0, 0.0], dtype=np.float64)      #Comet position within the grid
comet_velocity = np.array([-10, 10, 10.0], dtype=np.float64)        # Define the comet velocity in m/s


''' INPUTS - FOR PHOTOEMISSION, SECONDARY ELECTRON CURRENT AND SOLAR RADIATION (MATERIAL DEPENDENT + DISTANCE FROM THE SUN) '''


Wf = 4.8                              #Material dependent Work function (eV), used in thermionic current
gamma = 0.1                           #The quantum yield of the grain material, i.e. the number of photoelectrons emitted per incident photon.
delta_emax = 2.4                      #Secondary electron yield ( number of secondaries electrons per primary )
E_max = 400                           #Primary energy, energy at maximum yield (eV)
Q_pr = 1                              #Scattering efficiency for radiation pressure (adim)



'''Here there are two options : you can either provide inputs for the grain speed relative to the plasma or you can resume 
without it => choose a value for grain_motion '''

grain_motion = 'no'                   #'Yes' to take in account the grain speed relative to the plasma, 'no' to do the opposit. This impact on the ion current.
rel_velocity = 400e3                  # The dust-to-plasma relative velocity w (m/s)




'This section asks the user for the environment : comet or Earth. The gravitational force expression will be chosen accordingly.'


#Please specify if you want to work with the cometary dust or debris environment
environment = input("Please specify the environment : choose between comet/debris) : ").strip().lower() 




''' INITIAL VALUES OF POSITIONS, VELOCITIES AND ELECTROMAGNETIC FIELDS '''



# If the environment is a comet
if environment == "comet":
    
    #Definition of the distance from the sun and the comet
    distance_sun = 2e8
    distance_comet = 1e3
    
    # Calculate the comet's orbital velocity
    comet_velocity_orbit = np.sqrt(G * solar_mass / distance_sun)
    
    # Number of dust grains to simulate
    num_grains = 100
    initial_conditions = []

    # Loop to initialize each dust grain
    for _ in range(num_grains):
        
        # Random position around the comet within a spherical shell
        r = np.random.uniform(comet_radius, (comet_radius + 10))
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2 * np.pi)
        
        # Calculate the 3D position vector
        position = np.array([
            r * np.sin(theta) * np.cos(phi),
            r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta)
        ])
        
        # Random initial velocity for each grain within the range [-1000, 1000] m/s for each component
        velocity = np.random.uniform(-1e3, 1e3, 3)
        
        # Store initial position and velocity in the conditions list
        initial_conditions.append({'position': position, 'velocity': velocity})
        
        
        

else:  # If the environment is around the Earth
    
    # Number of debris particles to simulate
    num_debris = 10
    initial_conditions = []
    
    
    #Definition of the distance from the sun
    distance_sun = 1.5e11

    # Loop to initialize each debris particle
    for _ in range(num_debris):
        
        
        # Random position around the Earth within a spherical shell
        # Random altitude from 1000 km to 10000 km above Earth's surface
        r = np.random.uniform(earth_radius + 1e6, earth_radius + 1e7)
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2 * np.pi)
        
        # Calculate the 3D position vector
        initial_position = np.array([
            r * np.sin(theta) * np.cos(phi),
            r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta)
        ])
        
        # Calculate orbital velocity for circular orbit at this altitude
        orbital_velocity = np.sqrt(G * earth_mass / r)
        
        # Base vectors representing directions in 3D space
        base_vectors = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
        
        # Randomly choose one of the base vectors
        chosen_vector = random.choice(base_vectors)
        
        # Calculate the tangential velocity using the chosen vector and the orbital velocity
        tangential_velocity = np.cross(chosen_vector, initial_position / np.linalg.norm(initial_position)) * orbital_velocity
        
        # Store initial position and velocity in the conditions list
        initial_conditions.append({'position': initial_position, 'velocity': tangential_velocity})
        


magnetic_field = np.array([0, 0, 0])               # Magnetic field (T)
electric_field = np.array([0, 0, 0])               # Electric field (V/m)

# Define the initial charge of the dust according to the initial surface potential
initial_charge = 4 * np.pi * epsilon_0 * grain_radius * initial_potential

#Convert the distance in AU, used in photoemission expression
distance_sun_AU = distance_sun / AU




''' SET TIME STEPS PARAMETERS FOR SIMULATION '''


if environment == "comet":
    
    # Total simulation time in seconds for the Boris algorithm
    total_time_boris = 2  
    
    # Time step in seconds for the Boris algorithm
    dt_boris = 1e-2  
    
else:
    
    #Orbital frequency of the debris around the Earth, used to determine the total time of simulation
    orbital_frequency_Hz = ((G * earth_mass) / (4 * np.pi**2 * np.linalg.norm(initial_position)**3))**(1/2)
    
    # Time step in seconds for the Boris algorithm
    dt_boris = 1  
    
    # Total simulation time in seconds for the Boris algorithm
    total_time_boris = 100 

# Calculate the number of time steps for the Boris algorithm
num_steps_boris = int(total_time_boris / dt_boris) 

# Initialize a counter for tracking progress
counter = 0

