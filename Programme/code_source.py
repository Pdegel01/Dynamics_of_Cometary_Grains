import numpy as np
import math
import scipy.integrate as integrate
import sys


from constants import e, me, mi, kb, c, epsilon_0, Ar, G, earth_mass, solar_mass, earth_radius, moon_radius, J_2, sun_luminosity, AU

from inputs import *



''' DEFINITION OF THE FUNCTIONS FOR THE FORCES  '''

#Definition of the unit position vectors, the idea is to write only the magnitude of the force and then multiply it by the unit vectors
unit_vector_x = np.array([1, 0, 0])
unit_vector_y = np.array([0, 1, 0])
unit_vector_z = np.array([0, 0, 1])

#Definition of a vector to track the solar position, takes the same direction as the axis x according to GEI coordinates
solar_position = distance_sun * unit_vector_x




#Solar gravity force magnitude expression, reference : [Agarwal et al. 2023] [20]
def solar_gravity_magnitude(coeff_disable_solar, grain_mass, solar_mass):
    """
    Calculate the force of the solar gravity force acting on a grain of a given mass at a given distance from the Sun.

    Args:
    grain_mass (float): The mass of the grain, in kilograms.
    solar_position (array-like): The position of the sun, in Cartesian coordinates (x, y, z) in meters.

    Returns:
    array-like: The sun's gravity force, in Newtons, in Cartesian coordinates (F_x, F_y, F_z).
    """

    
    return coeff_disable_solar * (G * solar_mass * grain_mass) 


#Dimensionless coeff that can manipulated to zero the forces when in the validation process
coeff_disable_solar = 1



#Comet gravity force magnitude expression, reference : [Agarwal et al. 2023] [20]
def comet_gravity_force(comet_mass, grain_mass):

    """
    Calculate the force of the gravity force between a comet and a grain of a given mass at a given position.

    Args:
    comet_mass (float): The mass of the comet, in kilograms.
    grain_mass (float): The mass of the grain, in kilograms.
    positions (array-like): The position of the grain, in Cartesian coordinates (x, y, z) in meters.

    Returns:
    array-like: The Comet's gravity force, in Newtons, in Cartesian coordinates (F_x, F_y, F_z).
    """
    

    return -(G * comet_mass * grain_mass) 




def radial_tidal_force_magnitude(coeff_disable_tidal, distance_comet, grain_mass):
    
    return coeff_disable_tidal * (2 * G * solar_mass / distance_sun**3) * grain_mass * distance_comet

#Dimensionless coeff that can manipulated to zero the forces when in the validation process
coeff_disable_tidal = 1


#Solar radiation pressure force magnitude, reference : [Horanyi and Mendis, 1985] [19] 
def solar_radiation_pressure_force_magnitude(Q_pr, grain_radius):
    """
    Calculate the magnitude of the solar radiation pressure force acting on a grain of a given radius at a given distance from the Sun.

    Args:
    grain_radius (float): The radius of the grain, in meters.
    distance_sun (float): The distance between the grain and the Sun, in meters.

    Returns:
    float: The magnitude of the solar radiation pressure force, in Newtons.
    """
    
    return - ( Q_pr * sun_luminosity * grain_radius**2) / c




#Poynting-Robertson drag magnitude, reference : [Horanyi and Mendis, 1985] [19] and  [Agarwal et al. 2023] [20]
def poynting_robertson_drag_magnitude(coeff_disable_poynting, grain_radius):
    """
    Calculate the magnitude of the Poynting-Robertson drag force acting on a grain of a given radius at a given distance from the Sun.

    Args:
    grain_radius (float): The radius of the grain, in meters.
    distance_sun (float): The distance between the grain and the Sun, in meters.

    Returns:
    float: The magnitude of the Poynting-Robertson drag force, in Newtons.
    """

    return - coeff_disable_poynting * ((grain_radius**2 * sun_luminosity) / (4 * c**2)) * np.sqrt(G * solar_mass)

#Dimensionless coeff that can be manipulated to zero the P-R contribution
coeff_disable_poynting = 1


#Gas drag force magnitude, reference : [Agarwal et al. 2023] [20]
def gas_drag_force_magnitude(drag_coeff, gas_mol_mass, gas_number_density, grain_radius):
    """
    Calculate the magnitude of the gas drag force acting on a grain of a given radius and velocity relative to the gas.

    Args:
    gas_mol_mass (float): The molar mass of the gas, in kilograms per mole.
    gas_number_density (float): The number density of the gas, in particles per cubic meter.
    grain_radius (float): The radius of the grain, in meters.
    rel_velocity (float): The velocity of the grain relative to the gas, in meters per second.

    Returns:
    float: The magnitude of the gas drag force, in Newtons.
    """

    return - 0.5 * drag_coeff * gas_mol_mass * gas_number_density * np.pi * grain_radius**2 





#Earth gravitational force including the Moon's perturbation, reference : [JUHASZ et HORANYI,1997] [7]
def earth_gravity_force(grain_mass, position):
    """
    Calculate the Earth's gravity force acting on a grain of a given mass at a given position.

    Args:
    grain_mass (float): The mass of the grain, in kilograms.
    positions (array-like): The position of the grain, in Cartesian coordinates (x, y, z) in meters.

    Returns:
    array-like: The Earth's gravity force, in Newtons, in Cartesian coordinates (F_x, F_y, F_z).
    """


    # Extract the x, y, and z coordinates from the positions array
    x = position[0]
    y = position[1]
    z = position[2]

    # Calculate the distance between the grain and the center of the Earth
    r = np.linalg.norm(position)

    # Calculate the first term of the Earth's gravity force
    term1 = -G * earth_mass / r**3

    # Calculate the second term of the Earth's gravity force
    term2 = 1 - (3/2) * J_2 * (earth_radius / r)**2 * (3 * (z**2 / r**2) - 1)

    # Calculate the product of the first and second terms
    F1 = term1 * term2

    # Calculate the third term of the Earth's gravity force
    term3 = -3 * G * earth_mass * J_2 * (moon_radius / r)**2 / r**5

    # Calculate the x, y, and z components of the third term
    Fx = term3 * x * z**2
    Fy = term3 * y * z**2
    Fz = term3 * (-x**2 - y**2)

    # Calculate the x, y, and z components of the Earth's gravity force
    F_x = F1 * x + Fx
    F_y = F1 * y + Fy
    F_z = F1 * z + Fz


    return  np.array([F_x, F_y, F_z]) * grain_mass 



''' CHARGING CURRENTS '''


capacitance = 4 * np.pi * epsilon_0 * grain_radius          #Capacitance of the spherical grain (F)



# Expression of the electron current to a spherical grain (A), reference :  [Mann et al. 2014] [1]  
def electron_current(Te, ne, grain_radius, surface_potential):
    """
    Calculate the electron current on a spherical grain of a given radius and surface potential in a plasma with a given electron temperature and density.

    Args:
    Te (float): The electron temperature in the plasma, in Kelvin.
    ne (float): The electron density in the plasma, in particles per cubic meter.
    grain_radius (float): The radius of the spherical grain, in meters.
    surface_potential (float): The electric potential at the surface of the grain, in Volts.

    Returns:
    float: The electron current on the grain, in Amperes.
    """
    
    #For a negative potential
    if surface_potential < 0 :
        
        I_e = (4 * np.pi * grain_radius**2 * -e * ne) * (((kb * Te) / (2 * np.pi * me)) ** (1/2)) * np.exp((- e * np.abs(surface_potential)) / (kb * Te))
    
    #For a zero potential
    elif surface_potential == 0 :    
        
        I_e = (4 * np.pi * grain_radius**2 * -e * ne) * (((kb * Te) / (2 * np.pi * me)) ** (1/2)) 
    
    #For a positive potential
    else :                           
        
        I_e = (4 * np.pi * grain_radius**2 * -e * ne) * (((kb * Te) / (2 * np.pi * me)) ** (1/2)) * (1+((e * surface_potential) / (kb * Te))) 
            
    return I_e 




# Expression of the ion current (moving and static) to a spherical grain (A), reference : [Mann et al.] 2014 [1]  
def ion_current(Ti, ni, grain_radius, surface_potential, rel_velocity):
    """
    Calculate the ion current on a spherical grain of a given radius and surface potential in a plasma with a given ion temperature and density, and a given relative velocity between the grain and the plasma.
    Taking in account the relative grain-to-plasma velocity is optional. If grain_motion = 'yes' the input velocity is used and the expression adapted. 
    
    Args:
    Ti (float): The ion temperature in the plasma, in Kelvin.
    ni (float): The ion density in the plasma, in particles per cubic meter.
    grain_radius (float): The radius of the spherical grain, in meters.
    surface_potential (float): The electric potential at the surface of the grain, in Volts.
    rel_velocity (float): The relative velocity between the grain and the plasma, in meters per second.

    Returns:
    float: The ion current on the grain, in Amperes.
    """

    
    #Relative Mach number : the ratio of the dust-to-plasma relative velocity w over the ion thermal speed
    mach_number = rel_velocity / ((2 * kb * Ti) / mi)**(1/2)  

    # If you choose to take in account the motion of the grain relative to the plasma (grain_motion = 'yes')
    if grain_motion == 'yes':

        mach_factor = 1/2 * (((mach_number**2 + 1/2 - (e * surface_potential)/(kb * Ti)) * ((np.pi)**(1/2) / mach_number) * math.erf(mach_number)) + np.exp(-mach_number**2))

        # For a positive potential
        if surface_potential > 0:

            I_i = mach_factor * ( 4 * np.pi * grain_radius**2 * e * ni) * (((kb * Ti) / (2 * np.pi * mi)) ** (1/2)) * np.exp((-e * surface_potential) / (kb * Ti))

        # For a zero potential
        elif surface_potential == 0:

            I_i = mach_factor * (4 * np.pi * grain_radius**2 * e * ni) * (((kb * Ti) / (2 * np.pi * mi)) ** (1/2))

        # For a negative potential
        else:

            I_i = mach_factor * (4 * np.pi * grain_radius**2 * e * ni) * (((kb * Ti) / (2 * np.pi * mi)) ** (1/2)) * (1-((e * surface_potential) / (kb * Ti)))

        return I_i

    # If you choose to NOT take in account the motion of the grain relative to the plasma (grain_motion = 'no'), a default value will is computed for the script to work, but is not used
    else:

        # For a positive potential
        if surface_potential > 0:

            I_i =  (4 * np.pi * grain_radius**2 * e * ni) * (((kb * Ti) / (2 * np.pi * mi)) ** (1/2)) * np.exp((-e * surface_potential) / (kb * Ti))

        # For a zero potential
        elif surface_potential == 0:

            I_i =  (4 * np.pi * grain_radius**2 * e * ni) * (((kb * Ti) / (2 * np.pi * mi)) ** (1/2))

        # For a negative potential
        else:

            I_i = (4 * np.pi * grain_radius**2 * e * ni) * (((kb * Ti) / (2 * np.pi * mi)) ** (1/2)) * (1-((e * surface_potential) / (kb * Ti)))

        return I_i





#Expression of the photoemission current (A), reference : [Horanyi 1996] [5]

#f is the photoelectron constant flux
f = (2.5 * 10**14 *gamma) / distance_sun_AU**2  

def photoemission_current(Tph, surface_potential, grain_radius, gamma, distance_sun):
    """
    Calculate the photoelectron current on a spherical grain of a given radius and surface potential in a plasma with a given photoelectron temperature and a given incident solar flux.

    Args:
    Tph (float): The photoelectron temperature, in Kelvin.
    surface_potential (float): The electric potential at the surface of the grain, in Volts.
    grain_radius (float): The radius of the spherical grain, in meters.
    gamma (float): The quantum yield of the grain material, i.e. the number of photoelectrons emitted per incident photon.
    distance_sun (float): The distance between the grain and the Sun, in meters.

    Returns:
    float: The photoelectron current on the grain, in Amperes.
    """

    # For negative potential
    if surface_potential < 0:

        I_ph = np.pi * grain_radius**2 * e * f

    # For zero potential
    elif surface_potential == 0:

        I_ph = np.pi * grain_radius**2 * e * f

    # For positive potential
    else:
        
        I_ph = np.pi * grain_radius**2 * e * f  * np.exp(- e * surface_potential / (kb * Tph)) * (1 + ( e * surface_potential / (kb * Tph)))

    return I_ph





#Expression of the secondary electron emission (A), reference : [Horanyi 1996] [5]     

#Electron current for a zero potential, creation of the variable I_e to simplify the expression
I_e0 = - electron_current(Te, ne, grain_radius, surface_potential=0)     



#Function F5(x) to solve in the expression of the secondary electron current
def F_n_int(x, n, lower_limit):
    """
    Calculate the integral for the function F5.

    Args:
    x (float): The value of x for which to calculate the integral.
    n (int): The exponent of u in the integrand.
    lower_limit (float): The lower limit of integration.

    Returns:
    float: The value of the integral.
    """

    def integrand(u, x, n):
        """
        The integrand for the function F5(x).

        Args:
        u (float): The value of u for which to calculate the integrand.
        x (float): The value of x for which to calculate the integral.
        n (int): The exponent of u in the integrand.

        Returns:
        float: The value of the integrand.
        """

        return u**n * np.exp(-(x * u**2 + u))

    result, error = integrate.quad(integrand, lower_limit, np.inf, args=(x, n))
    
    return x**2 * result

#Argument of the F5 function and B the lower limit of the integral when the potential is positive
arg = E_max / ((4 * kb * Te)/e)
B_5 = ((e * np.abs(initial_potential) / (kb * Te)) / arg)**(1/2)


#F5 is integrated from 0 to infinity if the potential is negative
F_5_0 = F_n_int(arg, 4, 0)

#F5 is integrated from B to infinity if the potential is positive
F_5_B  = F_n_int(arg, 4, B_5)


def secondary_emission_current(Tsec, Te, delta_emax, E_max, grain_radius, surface_potential):
    """
    Calculate the secondary electron current on a spherical grain of a given radius and surface potential in a plasma with a given temperature and a given incident electron flux.

    Args:
    Tsec (float): The temperature of the secondary electrons, in Kelvin.
    Te (float): The temperature of the incident electrons, in Kelvin.
    delta_emax (float): Secondary electron yield ( number of secondaries electrons per primary ).
    E_max (float): Primary energy, energy at maximum yield in eV.
    grain_radius (float): The radius of the spherical grain, in meters.
    surface_potential (float): The electric potential at the surface of the grain, in Volts.

    Returns:
    float: The secondary electron current on the grain, in Amperes.
    """

    # For negative potential
    if surface_potential < 0:

        I_sec = 3.7 * I_e0 * delta_emax * F_5_0 * np.exp((e * surface_potential) / (kb * Te))

    # For zero potential
    elif surface_potential == 0:

        I_sec = 3.7 * I_e0 * delta_emax * F_5_0

    # For positive potential
    else:

        I_sec = 3.7 * I_e0 * delta_emax * F_5_B * np.exp((- e * surface_potential) / (kb * Tsec)) * (1 + ((e * surface_potential) / (kb * Tsec)))

    return I_sec





#Expression of the thermionic current (A),  reference : [Zhuang Liu et al., 2017 ] [21]
def thermionic_emission_current(Te, grain_temperature, Wf, grain_radius):
    """
    Calculate the thermionic current on a spherical grain of a given radius and temperature in a plasma with a given temperature.

    Args:
    Te (float): The temperature of the plasma, in Kelvin.
    grain_temperature (float): The temperature of the spherical grain, in Kelvin.
    Wf (float): The work function of the grain material, in eV.
    grain_radius (float): The radius of the spherical grain, in meters.

    Returns:
    float: The thermionic current on the grain, in Amperes.
    """

    I_th = 4 * np.pi * grain_radius**2 * Ar * grain_temperature**2 * np.exp(- Wf * e / (kb * grain_temperature))

    return I_th




#Expression of the backscattered electrons current (A), reference : [Horanyi,1989] [17]
def backscattered_emission_current(Te, ne, z_number, surface_potential):
    """
    Calculate the backscattered electron current on a spherical grain of a given radius and surface potential in a plasma with a given temperature and density.

    Args:
    Te (float): The temperature of the plasma, in Kelvin.
    ne (float): The density of the plasma electrons, in m^-3.
    z_number (float): The atomic number of the grain material.
    surface_potential (float): The electric potential at the surface of the grain, in Volts.

    Returns:
    float: The backscattered electron current on the grain, in Amperes.
    """

    #E_bs is the resulting energy of the backscattered electrons in eV
    E_bs = ((0.45 + 2*10**-3 * z_number) * (kb * Te + (e * surface_potential)))/ e
    
    #E_min is the energy distribution of the backscattered electrons in eV
    E_min = 100 
    
    # Look again for the condition
    if E_bs < E_min:

        I_bs = 0

    else:

        I_bs = 0.3 * electron_current(Te, ne, grain_radius, surface_potential)

    return I_bs




#Sum of the currents, we will need this when solving the surface potential at equilibrium and when calculating the charge
def total_current(surface_potential):
    """
    Calculate the total current on a spherical grain of a given radius and surface potential in a plasma with a given temperature and density.

    Args:
    surface_potential (float): The electric potential at the surface of the grain, in Volts.
    All of the previous arguments are included implicitly within the currents. 
        
    Returns:
    float: The total current on the grain, in Amperes.
    """

    total_current = electron_current(Te, ne, grain_radius, surface_potential) + ion_current(Ti, ni, grain_radius, surface_potential, rel_velocity) + photoemission_current(Tph, surface_potential, grain_radius, gamma, distance_sun) + secondary_emission_current(Tsec, Te, delta_emax, E_max, grain_radius, surface_potential) + backscattered_emission_current(Te, ne, z_number, surface_potential) + thermionic_emission_current(Te, grain_temperature, Wf, grain_radius)

    return total_current




''' INTEGRATION FUNCTIONS '''



'Estimation of the surface potential at equilibrium with the dichotomy method '

#Initialization of the initial interval bounds a and b : the interval needs to contain the solution we are looking for
a = -1000.0
b = 1000.0
tolerance = 1e-6                  #Defines how close we want to get to the actual solution, when the numerical solution has the same value with a 1e-6 we consider that equilibrium has been reached
max_iterations = 1000             #Maximum number of iterations


# Function to be solved : sum of the currents = 0, this function is created to later evaluate the sum of the currents for diferent values of potential 
def equation_to_solve(surface_potential):
    
    'The goal is to solve for the potential when the sum of the currents is zero'
    
    return total_current(surface_potential)




# Bisection method to solve the equation
def solve_with_bisection(a, b):
    """
    INPUT : a, b = bounds of the interval
    OUTPUT : surface_potential_solution = potential at equilibrium
    """
    # Check if the initial bounds are valid (root should be between them)
    if equation_to_solve(a) * equation_to_solve(b) >= 0:
       
        print("Les bornes initiales ne sont pas valides.")  # Print message if bounds are not valid, i.e., if there is no solution
        return None

    # Initialize surface_potential as the midpoint of a and b
    surface_potential = (a + b) / 2.0

    # Initialize iteration counter
    iterations = 0

    # Iterate until the solution is within the desired tolerance or max_iterations is reached
    while abs(b - a) > tolerance and iterations < max_iterations:
        
        # Update surface_potential as the midpoint of a and b
        surface_potential = (a + b) / 2.0

        # Determine in which half the root lies and update a or b accordingly
        if equation_to_solve(a) * equation_to_solve(surface_potential) < 0:
            b = surface_potential
        else:
            a = surface_potential

        # Increment the iteration counter
        iterations += 1

    # Print warning if max_iterations reached without finding a solution (you should increase the number of iterations)
    if iterations >= max_iterations:
        
        print("Nombre maximum d'itérations atteint. La solution peut ne pas être précise.")

    # Return the approximated solution
    return surface_potential




 # Call the bisection method to find the solution
surface_potential_solution = solve_with_bisection(a, b)                                             

# Print the solution if found
if surface_potential_solution is not None:
    print("The surface potential solution estimated with the dichotomy method :", surface_potential_solution, 'Volts')


''''RUNGE-KUTTA to estimate the surface potential and the charge in time'''


#Definition of parameters for integration
dt_RK4 = 0.01                                        #Time step
t_final_RK4 = 500                                    #Final time (seconds)
num_steps_RK4 = int(t_final_RK4 / dt_RK4)                    # Nombre de pas



#Definition of a function that realizes one step of the Runge-Kutta, later we will use this function and repeat it over time
def rk4_step(current_potential, dt):
    """
    One step of RK4 (Runge-Kutta 4th order) integration.

    Args:
    - current_potential (float): The potential at the current time step.
    - dt (float): The time step size.

    Returns:
    - delta_v (float): The change in potential over the time step dt using the RK4 method.
    """

    # Calculate the first slope (k1) based on the current potential
    k1 = total_current(current_potential) / capacitance

    # Calculate the second slope (k2) using the first slope
    k2 = total_current(current_potential + 0.5 * dt * k1) / capacitance

    # Calculate the third slope (k3) using the second slope
    k3 = total_current(current_potential + 0.5 * dt * k2) / capacitance

    # Calculate the fourth slope (k4) using the third slope
    k4 = total_current(current_potential + dt * k3) / capacitance

    # Combine the four slopes to get the weighted average slope (delta_v)
    delta_v = (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6.0

    return delta_v




def update_charge_and_potential(surface_potential, dt):
    """
    Updates the surface potential and the charge of the particle using RK4 (Runge-Kutta 4th order method).

    Args:
    - surface_potential (float): The current surface potential of the particle.
    - dt (float): The time step.

    Returns:
    - new_potential (float): The new surface potential.
    - new_charge (float): The new charge of the particle.
    """
    
    # Calculate the change in potential using the RK4 step method.
    delta_v = rk4_step(surface_potential, dt)
    
    # Update the surface potential by adding the change in potential.
    new_potential = surface_potential + delta_v
    
    # Calculate the new charge of the particle using the capacitance and the new potential.
    new_charge = capacitance * new_potential
    
    return new_potential, new_charge




'''BORIS ALGORITHM TO ESTIMATE THE POSITIONS AND VELOCITIES IN TIME'''

test_earth = []
test_drag = []

# Boris pusher function to that gives the new position and velocity for each time step
def boris_push(position, velocity, magnetic_field, electric_field, charge, mass, dt, environment, surface_potential):
    """
    Boris algorithm for the motion of charged particles in an electromagnetic field and other forces, including charge evolution.

    Args:
    position (array-like): Position vector of the particle.
    velocity (array-like): Velocity vector of the particle.
    magnetic_field (array-like): Magnetic field vector.
    electric_field (array-like): Electric field vector.
    charge (float): Charge of the grain.
    mass (float): Mass of the grain.
    dt (float): Time step.
    environment (string): Environment chosen (this affects only the expression of the gravitational force led by the comet or earth)
    surface_potential (float): Surface potential of the grain.

    Returns:
    tuple: Updated position vector, velocity vector, and surface potential.
    """
    
    # Update charge and surface potential
    surface_potential, charge = update_charge_and_potential(surface_potential, dt)
    
     
    
    # Select the gravitational force according to the environment
    if environment == "comet":
        
        comet_force = comet_gravity_force(comet_mass, grain_mass) * (position) / np.linalg.norm(position)**3
        radial_tidal_force = radial_tidal_force_magnitude(coeff_disable_tidal, solar_mass, grain_mass) * (position / np.linalg.norm(position)) / np.linalg.norm(solar_position - position)**3
        
        gravitational_force = comet_force + radial_tidal_force
        test_earth.append(gravitational_force)
    else:
        
        gravitational_force = earth_gravity_force(grain_mass, position)
        test_earth.append(gravitational_force)
    
    #Multiply the magnitude of the solar gravity force by the corresponding vectors to obtain the expression of the force
    solar_gravity_force = solar_gravity_magnitude(coeff_disable_solar, grain_mass, solar_mass) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3 - (solar_position / np.linalg.norm(solar_position)**3))
    
    #Multiply the magnitude of the solar radiation pressure force by the corresponding vectors to obtain the expression of the force
    solar_radiation_pressure_force = solar_radiation_pressure_force_magnitude(Q_pr, grain_radius) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3) 
    
    #Multiply the magnitude of the Poyting-Robertson force by the corresponding vectors to obtain the expression of the force
    poynting_robertson_drag_force = (poynting_robertson_drag_magnitude(coeff_disable_poynting, grain_radius)) * (velocity / np.linalg.norm(velocity)) * (1 / np.sqrt(np.linalg.norm(solar_position - position)**5)) 
    
    #Multiply the magnitude of the gas drag force by the corresponding vectors to obtain the expression of the force
    gas_drag_force = gas_drag_force_magnitude(drag_coeff, gas_mol_mass, gas_number_density, grain_radius) * np.linalg.norm(abs(velocity - gas_velocity)) * abs(velocity - gas_velocity)
    test_drag.append(gas_drag_force)
    # Total forces, divided by the charge to obtain something homogeneous with velocity
    total_force = (gravitational_force + solar_gravity_force + solar_radiation_pressure_force + poynting_robertson_drag_force + gas_drag_force) / charge
    

    # Update the velocity by half-step with electric impulsion
    v_minus = velocity + (charge / mass) * (electric_field + total_force) * (dt / 2.0)

    # Rotate the velocity vector by half-step with the magnetic field
    t = (charge / mass) * magnetic_field * (dt / 2.0)
    t_mag2 = np.dot(t, t)
    s = 2 * t / (1 + t_mag2)

    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)

    # Update the final velocity half-step
    velocity = v_plus + (charge / mass) * (electric_field + total_force) * (dt / 2.0)

    # Update the position
    position += velocity * dt

    return position, velocity, surface_potential, charge







        

