"""
This program is designed to simulate and validate the forces acting on charged particles
 (such as dust grains) in various space environments (such as around a comet or in Earth's orbit).

                                                      

                                                      
   MAIN OBJECTIVES :
       
1. Simulate the Movement of Charged Particles:
    
   - Use the Boris algorithm to estimate the positions and velocities of charged particles based 
   on electromagnetic and gravitational forces.

2. Calculate Charging Currents:
    
   - Calculate electron, ion, photoemission, secondary emission, thermionic emission, 
   and backscattered currents, which influence the charge and surface potential of the particles.

3. Validate Implemented Forces:
    
   - Compare simulation results (using the Boris algorithm) with analytical solutions to 
   validate the implementation of different forces:
     - Solar gravity
     - Earth's gravity
     - Comet's gravity
     - Radial tidal force
     - Solar radiation pressure
     - Gas drag force
     - Poynting-Robertson effect




   PROGRAM STRUCTURE :
       
1. Definition of Classic Physical Constants:
    
   - Elementary charge, electron mass, proton mass, Boltzmann's constant, speed of light, etc.

2. Definition of Input Parameters:
    
   - Electron and ion temperatures and densities, grain characteristics (radius, mass, surface potential), 
   and environmental characteristics (mass and radius of the comet or Earth).

3. Force Calculation Functions:
    
   - Functions to calculate various forces acting on the grains, such as gravity, radiation pressure, gas drag, etc.

4. Charging Current Calculation Functions:
    
   - Functions to calculate electron, ion, photoemission, secondary emission, thermionic emission, and backscattered currents.

5. Simulation Algorithms:
    
   - Use the Boris algorithm to simulate the movement of charged particles.
   - Bisection method to find the equilibrium surface potential.
   - Runge-Kutta method to estimate surface potential and charge over time.

6. Force Validation:
    
   - Compare simulation results with analytical solutions to validate the implementation of the forces.

7. Force Selection for Validation:
    
   - Select the force to validate by disabling other forces and calling the corresponding validation function.
"""




import matplotlib.pyplot as plt       # Importing the plotting module
import numpy as np                    # Importing the numerical computation module                       
import scipy.integrate as integrate   # Importing the numerical integration module
import sys                            # Importing the system module
import importlib

import inputs
import constants
import code_source

from constants import *
from inputs import *
from code_source import solar_position





first_grain_initial_conditions = initial_conditions[0]
initial_position = first_grain_initial_conditions['position']
initial_velocity = first_grain_initial_conditions['velocity']

altitude = np.linalg.norm(initial_position) - earth_radius
altitude_km = int(altitude * 10**-3)




''' USER CHOICE FOR THE ENVIRONMENT AND THE VERIFICATION TEST'''





# Ask the user to specify whether they want to work with forces or charge evolution and implementation in Boris
verification = input("Please specify the verification: choose between forces / charge :").strip().lower()

# Based on the selected environment, further specify the type of verification
if environment == "debris":
    
    # This section selects the test validation according to the answer given above
    if verification == "forces":
        
        print("Forces validation is chosen to proceed")
        
        # Ask which specific force to analyze
        validation = input("Which force would you like to analyze? (earth g / solar g / radiation pressure / poynting / drag): ").strip().lower()
    
    elif verification == "charge":
        
        print("Charge evolution and implementation is chosen to proceed")
        
    
    else:
        
        print("Error: unrecognized verification. Please enter 'force' or 'charge'.")
        sys.exit()
else:
    
    # This section selects the test validation according to the answer given above
    if verification == "forces":
        
        print("Forces validation is chosen to proceed")
        
        # Ask which specific force to analyze
        validation = input("Which force would you like to analyze? (comet g / solar g / radial tidal / radiation pressure / poynting / drag): ").strip().lower()
   
    elif verification == "charge":
        
        print("Charge evolution and implementation is chosen to proceed")
    
    
    else:
        
        print("Error: unrecognized verification. Please enter 'force' or 'charge'.")
        sys.exit()





''' VALIDATION OF THE FORCES IMPLEMENTATION '''

"Definition of test functions, each function corresponds to the validation of a force"

# Definition of a time interval for the plots
t = np.arange(num_steps_boris + 1) 

if verification == "forces":

    # Function for validating the solar gravity force implementation
    def solar_gravity_validation():
        counter = 0
        
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            # Calculate the acceleration due to solar gravity
            acceleration = (code_source.solar_gravity_magnitude(code_source.coeff_disable_solar, grain_mass, solar_mass) / grain_mass) * ((solar_position - position) / np.linalg.norm(solar_position - position)**3 - (solar_position / np.linalg.norm(solar_position)**3))
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
            
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            
            
            return new_position, new_velocity
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
        
        
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
        
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
        
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        if environment == "comet":
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Solar gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
            
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Solar gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
        
        else:
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Solar gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()
            
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Solar gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()
            
            

    # Function for validating the solar radiation pressure force implementation
    def radiation_pressure_validation():
        counter = 0
        
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            # Calculate the acceleration due to solar radiation pressure
            acceleration = code_source.solar_radiation_pressure_force_magnitude(Q_pr, grain_radius) / grain_mass * ((solar_position - position) / np.linalg.norm(solar_position - position)**3)
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
            
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            return new_position, new_velocity
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
        
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
        
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
        
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        if environment == "comet":
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Radiation Pressure test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
            
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Radiation pressure test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
        else:
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Radiation Pressure test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()
            
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Radiation pressure test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()

    # Function for validating the Earth's gravity force implementation
    def earth_gravity_validation():
        counter = 0
        
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            acceleration = code_source.earth_gravity_force(grain_mass, position) / grain_mass
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
            
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            return new_position, new_velocity
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
        
        # Function to perform a single step of the Runge-Kutta 4th order (RK4) method
        def rk4_step(f, y, t, dt):
            k1 = dt * f(y, t)
            k2 = dt * f(y + 0.5 * k1, t + 0.5 * dt)
            k3 = dt * f(y + 0.5 * k2, t + 0.5 * dt)
            k4 = dt * f(y + k3, t + dt)
            return y + (k1 + 2*k2 + 2*k3 + k4) / 6
    
        # Function to define the equations of motion for the RK4 method
        def equations_of_motion(state, t):
            position = state[:3]
            velocity = state[3:]
            acceleration = code_source.earth_gravity_force(grain_mass, position) / grain_mass
            return np.concatenate((velocity, acceleration))
    
        # Function to simulate the orbit using the RK4 method
        def simulate_orbit(initial_state, t, dt):
            state = np.zeros((num_steps_boris, 6))
            state[0] = initial_state
            
            t = 0
            counter = 0
            for i in range(1, num_steps_boris):
                counter += 1
                print(f"Step {counter} of {num_steps_boris}")
                state[i] = rk4_step(equations_of_motion, state[i-1], t, dt_boris)
                t += dt
            
            return state
    
        initial_state = np.concatenate((initial_position, initial_velocity))
    
        # Run simulation
        trajectory = simulate_orbit(initial_state, t, dt_boris)
        
        # Extract velocities from the RK4 trajectory
        velocities_rk4 = [initial_velocity]  # Initialize with the initial velocity
        velocities_rk4.extend(trajectory[:, 3:])  # Append velocities from the trajectory
    
        # Convert list to numpy array for calculations
        velocities_rk4 = np.array(velocities_rk4)
    
        # Compute the norm of the velocity at each time step
        norm_velocities_rk4 = np.linalg.norm(velocities_rk4, axis=1)
    
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
        
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
        
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        # Plot the speed evolution
        plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
        plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
        plt.plot(t, norm_velocities_rk4, label='Velocity RK4')
        plt.xlabel('Time steps (1 orbital period)')
        plt.ylabel('Velocity (m/s)')
        plt.title(f"Velocity norm in time, Earth gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
        plt.legend()
        plt.show()
        
        # Plot the relative error between analytical and Boris method velocities
        errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
        plt.plot(t, errors, label='Relative error')
        plt.xlabel('Time steps (1 orbital period)')
        plt.ylabel('Relative error')
        plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Earth gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
        plt.legend()
        plt.show()

    # Function for validating the Poynting-Robertson drag force implementation
    def poynting_robertson_validation():
        counter = 0
       
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            # Compute the acceleration due to Poynting-Robertson drag
            acceleration = (code_source.poynting_robertson_drag_magnitude(code_source.coeff_disable_poynting, grain_radius) / grain_mass) * (velocity / np.linalg.norm(velocity))  / np.sqrt(np.linalg.norm(solar_position - position)**5)
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
        
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            return new_position, new_velocity
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
           
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
       
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
               
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        if environment == "comet":
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Poynting-Robertson test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
            
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Poynting Robertson test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
        else:
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Poynting-Robertson test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()
            
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Poynting Robertson test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()

    # Function for validating the gas drag force implementation
    def gas_drag_validation():
        counter = 0
       
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            # Compute the acceleration due to gas drag force
            acceleration = (code_source.gas_drag_force_magnitude(code_source.drag_coeff, gas_mol_mass, gas_number_density, grain_radius) / grain_mass) * np.linalg.norm(abs(velocity - gas_velocity)) * (abs(velocity - gas_velocity))
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
        
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            return new_position, new_velocity, acceleration
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity, _ = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
           
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
       
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
               
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        if environment == "comet":
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic , label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Gas drag test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
        
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.plot(t, errors, label='Relative error')
            plt.xlabel('Time steps')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Gas drag test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
            plt.legend()
            plt.show()
        else:
            # Plot the speed evolution
            plt.plot(t, norm_velocities_analytic , label='Velocity (analytical)')
            plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Velocity (m/s)')
            plt.title(f"Velocity norm in time, Gas drag test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()
        
            # Plot the relative error between analytical and Boris method velocities
            errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
            plt.loglog(t, errors, label='Relative error')
            plt.xlabel('Time steps (1 orbital period)')
            plt.ylabel('Relative error')
            plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Gas drag test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg, altitude = {altitude_km}km")
            plt.legend()
            plt.show()
    
    # Function for validating the comet gravity force implementation
    def comet_gravity_validation():
        counter = 0
        
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            # Calculate the acceleration due to comet gravity
            acceleration = (code_source.comet_gravity_force(comet_mass, grain_mass) / grain_mass) * (position) / np.linalg.norm(position)**3
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
            
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            return new_position, new_velocity
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
            
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
        
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
        
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        # Plot the speed evolution
        plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
        plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
        plt.xlabel('Time steps')
        plt.ylabel('Velocity (m/s)')
        plt.title(f"Velocity norm in time, Comet gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
        plt.legend()
        plt.show()
        
        # Plot the relative error between analytical and Boris method velocities
        errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
        plt.plot(t, errors, label='Relative error')
        plt.xlabel('Time steps')
        plt.ylabel('Relative error')
        plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Comet gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
        plt.legend()
        plt.show()

    # Function for validating the radial tidal force implementation
    def radial_tidal_validation():
        counter = 0
        
        # Function to compute the analytical trajectory
        def analytical_trajectory(t, position, velocity):
            # Calculate the acceleration due to radial tidal force
            acceleration = (code_source.radial_tidal_force_magnitude(code_source.coeff_disable_tidal, solar_mass, grain_mass) / grain_mass) * (position / np.linalg.norm(position)) / np.linalg.norm(solar_position - position)**3
            
            # Calculate the new position
            new_position = position + velocity * t + 0.5 * acceleration * t**2
            
            # Calculate the new velocity
            new_velocity = velocity + acceleration * t
            
            return new_position, new_velocity
    
        # Initialize position and velocity arrays for the analytical calculation
        positions_analytic = [initial_position]
        velocities_analytic = [initial_velocity]
        
        current_position = initial_position.copy()
        current_velocity = initial_velocity.copy()
        
        # Initialize position and velocity arrays for the Boris algorithm
        positions = [initial_position]
        velocities = [initial_velocity]
        charges = [initial_charge]
        potentials = [initial_potential]
        
        position = initial_position.copy()
        velocity = initial_velocity.copy()
        charge = initial_charge
        potential = initial_potential

        # Perform the analytical calculation for each time step
        for _ in range(num_steps_boris):
            # Compute new position and velocity using the analytical method
            new_position, new_velocity = analytical_trajectory(dt_boris, current_position, current_velocity)
            current_position, current_velocity = new_position, new_velocity
            positions_analytic.append(current_position)
            velocities_analytic.append(current_velocity)
            
            # Counter to track the running of the simulation
            counter += 1
            print(f"Step {counter} of {num_steps_boris}")
            
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
            positions.append((position[0], position[1], position[2]))
            velocities.append((velocity[0], velocity[1], velocity[2]))
            charges.append(charge)
            potentials.append(potential)
            
        # Convert positions to a numpy array for easier plotting
        positions = np.array(positions)
        velocities = np.array(velocities)
        
        # Convert positions and velocities to numpy arrays for easier manipulation
        positions_analytic = np.array(positions_analytic)
        velocities_analytic = np.array(velocities_analytic)
        
        # Calculate the speed over time using the analytical method
        norm_velocities_analytic = np.linalg.norm(velocities_analytic, axis=1)
        
        # Calculate the speed over time using the Boris algorithm
        norm_velocities_boris = np.linalg.norm(velocities, axis=1)
        
        # Plot the speed evolution
        plt.plot(t, norm_velocities_analytic, label='Velocity (analytical)')
        plt.plot(t, norm_velocities_boris, label='Velocity Boris', linestyle=':', linewidth=4)
        plt.xlabel('Time steps')
        plt.ylabel('Velocity (m/s)')
        plt.title(f"Velocity norm in time, Radial tidal test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
        plt.legend()
        plt.show()
        
        # Plot the relative error between analytical and Boris method velocities
        errors = abs((norm_velocities_boris - norm_velocities_analytic)) / norm_velocities_analytic
        plt.plot(t, errors, label='Relative error')
        plt.xlabel('Time steps')
        plt.ylabel('Relative error')
        plt.title(f"Relative error - Velocity norm Analytic/Boris in time, Radial tidal gravity test \n dt = {dt_boris}, radius = {grain_radius}m, mass ={grain_mass}kg")
        plt.legend()
        plt.show()

elif verification == "charge" :
    
    # Calculate the total energy (kinetic + potential)
    def calculate_total_energy(position, velocity, mass, charge):
        kinetic_energy = 0.5 * mass * np.linalg.norm(velocity)**2
        potential_energy = charge * np.linalg.norm(position)  # Simplified for this example
        return kinetic_energy + potential_energy
    
    # Simulation constants
    comet_position = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    
    if environment == "comet":
        # Initialization of variables
        initial_position = np.array([1.0, 1.0, 1.0])
        initial_velocity = np.array([1.0, 1.0, 1.0])
    else:
        # Initialization of variables
        initial_position = np.array([altitude, 0, 0])
        initial_velocity = np.array([1, orbital_velocity, 0])
    
    # Initialization of variables for Full and Half Update
    initial_position_full = initial_position.copy()
    initial_velocity_full = initial_velocity.copy()
    initial_position_half = initial_position.copy()
    initial_velocity_half = initial_velocity.copy()
    surface_potential_full = surface_potential_half = surface_potential
    grain_charge_full = grain_charge_half = initial_charge
    
    time_values = np.arange(total_time_boris) * dt_boris
    
    positions_full = []
    velocities_full = []
    charges_full = []
    potentials_full = []
    
    positions_half = []
    velocities_half = []
    charges_half = []
    potentials_half = []
    
    # Simulation for Half Update
    for step in range(num_steps_boris):
        initial_position_half, initial_velocity_half, surface_potential_half, grain_charge_half = boris_push_half_update(
            initial_position_half, initial_velocity_half, magnetic_field, electric_field, grain_charge_half, grain_mass, dt_boris, environment, surface_potential_half)
        positions_half.append(initial_position_half.copy())
        velocities_half.append(initial_velocity_half.copy())
        charges_half.append(grain_charge_half)
        potentials_half.append(surface_potential_half)
        
        # Simulation for Full Update
        initial_position_full, initial_velocity_full, surface_potential_full, grain_charge_full = boris_push_full_update(
            initial_position_full, initial_velocity_full, magnetic_field, electric_field, grain_charge_full, grain_mass, dt_boris, environment, surface_potential_full)
        positions_full.append(initial_position_full.copy())
        velocities_full.append(initial_velocity_full.copy())
        charges_full.append(grain_charge_full)
        potentials_full.append(surface_potential_full)
    
    # Calculate relative errors on velocity
    relative_errors = np.abs((np.linalg.norm(velocities_full, axis=1) - np.linalg.norm(velocities_half, axis=1)) / np.linalg.norm(velocities_half, axis=1))
    
    # Calculate relative errors on charge
    relative_charge_errors = np.abs((np.array(charges_full) - np.array(charges_half)) / np.array(charges_half))
    
    # Calculate total energy
    total_energy_full = [calculate_total_energy(pos, vel, grain_mass, charge) for pos, vel, charge in zip(positions_full, velocities_full, charges_full)]
    total_energy_half = [calculate_total_energy(pos, vel, grain_mass, charge) for pos, vel, charge in zip(positions_half, velocities_half, charges_half)]
    
    # Calculate total momentum
    momentum_full = [grain_mass * np.linalg.norm(vel) for vel in velocities_full]
    momentum_half = [grain_mass * np.linalg.norm(vel) for vel in velocities_half]
    
    # Plot individual velocities with adjusted scale
    fig, axs = plt.subplots(3, 1, figsize=(15, 15))
    
    time_values = np.linspace(0, num_steps_boris * dt_boris, num_steps_boris)
    
    # Calculate the velocity over time using the Boris algorithm
    norm_velocities_half = [np.linalg.norm(vel) for vel in velocities_half]
    norm_velocities_full = [np.linalg.norm(vel) for vel in velocities_full]
    
    # Plot Half Update and Full Update velocities
    axs[0].plot(time_values, norm_velocities_half, label='Particle Velocity (Half Update)', linestyle='--')
    axs[0].plot(time_values, norm_velocities_full, label='Particle Velocity (Full Update)', linestyle='-')
    axs[0].set_xlabel('Time (s)', fontsize=18)
    axs[0].set_ylabel('Velocity (m/s)', fontsize=18)
    axs[0].set_title('Particle Velocity (Half Update and Full Update)', fontsize=20)
    axs[0].legend(fontsize=16)
    axs[0].tick_params(axis='both', which='major', labelsize=16)
    
    # Plot relative errors with thinner lines
    axs[1].semilogy(time_values, relative_errors, linewidth=1)
    axs[1].set_xlabel('Time (s)', fontsize=18)
    axs[1].set_ylabel('Relative Error', fontsize=18)
    axs[1].set_title('Relative Error between Full Update and Half Update Methods', fontsize=20)
    axs[1].tick_params(axis='both', which='major', labelsize=16)
    
    # Plot total energy and total momentum evolution
    axs[2].semilogy(time_values, total_energy_half, label='Total Energy (Half Update)')
    axs[2].semilogy(time_values, total_energy_full, label='Total Energy (Full Update)', linestyle='--')
    axs[2].semilogy(time_values, momentum_half, label='Total Momentum (Half Update)', linestyle='-.')
    axs[2].semilogy(time_values, momentum_full, label='Total Momentum (Full Update)', linestyle=':')
    axs[2].set_xlabel('Time (s)', fontsize=18)
    axs[2].set_ylabel('Total Energy / Momentum (J / kg*m/s)', fontsize=18)
    axs[2].set_title('Total Energy and Total Momentum over Time', fontsize=20)
    axs[2].legend(fontsize=16)
    axs[2].tick_params(axis='both', which='major', labelsize=16)
    axs[2].set_xlim(left=1)  # Start at 10^0
    
    # Adjust the spacing between subplots
    plt.subplots_adjust(hspace=0.5)  # Space of 0.5 inch (~1.27 cm)
    
    # Plot charge and charge errors
    fig_charge, axs_charge = plt.subplots(2, 1, figsize=(15, 10))
    
    # Plot Half Update and Full Update charges
    axs_charge[0].plot(time_values, charges_half, label='Charge (Half Update)', linestyle='--')
    axs_charge[0].plot(time_values, charges_full, label='Charge (Full Update)', linestyle='-')
    axs_charge[0].set_xlabel('Time (s)', fontsize=18)
    axs_charge[0].set_ylabel('Charge (C)', fontsize=18)
    axs_charge[0].set_title('Charge (Half Update and Full Update)', fontsize=20)
    axs_charge[0].legend(fontsize=16)
    axs_charge[0].tick_params(axis='both', which='major', labelsize=16)
    
    # Plot relative errors on charge
    axs_charge[1].semilogy(time_values, relative_charge_errors, linewidth=1)
    axs_charge[1].set_xlabel('Time (s)', fontsize=18)
    axs_charge[1].set_ylabel('Relative Charge Error', fontsize=18)
    axs_charge[1].set_title('Relative Charge Error between Full Update and Half Update Methods', fontsize=20)
    axs_charge[1].tick_params(axis='both', which='major', labelsize=16)
    
    # Adjust the spacing between subplots
    plt.subplots_adjust(hspace=0.5)  # Space of 0.5 inch (~1.27 cm)
    
    plt.tight_layout()
    plt.show()



    
'''CHOIX DE LA FORCE'''

if verification =="forces" :
    # This section selects the gravitational force according to the answer given above 
    if validation == "radial tidal":
        
        # Disable other forces
        code_source.comet_mass = 0                            # Disable comet's gravity
        code_source.earth_mass = 0                            # Disable Earth's gravity
        code_source.Q_pr = 0                                  # Disable radiation pressure
        code_source.drag_coeff = 0                            # Disable gas drag force
        code_source.coeff_disable_poynting = 0                # Disable Poynting-Robertson effect
        code_source.coeff_disable_solar = 0                   # Disable Solar gravity
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        print("Validation of the Radial tidal force implementation. All other forces are zero.")
        
        # Call the Earth gravity validation function
        radial_tidal_validation()
    
    # This section selects the gravitational force according to the answer given above 
    elif validation == "comet g":
        
        # Disable other forces
        code_source.coeff_disable_tidal = 0                   # Disable tidal effect
        code_source.earth_mass = 0                            # Disable Earth's gravity
        code_source.solar_mass = 0                            # Disable Solar gravity
        code_source.Q_pr = 0                                  # Disable radiation pressure
        code_source.drag_coeff = 0                            # Disable gas drag force
        code_source.coeff_disable_poynting = 0                # Disable Poynting-Robertson effect
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        print("Validation of the comet's gravity implementation. All other forces are zero.")
        
        # Call the Earth gravity validation function
        comet_gravity_validation()
    
    elif validation == "earth g":
        
        # Disable other forces
        code_source.solar_mass = 0                            # Disable Solar gravity
        code_source.Q_pr = 0                                  # Disable radiation pressure
        code_source.drag_coeff = 0                            # Disable gas drag force
        code_source.coeff_disable_poynting = 0                # Disable Poynting-Robertson effect
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        print("Validation of the earth's gravity implementation. All other forces are zero.")
        
        # Call the Earth gravity validation function
        earth_gravity_validation()
        
    elif validation == "solar g":
        print("Validation of the solar's gravity implementation. All other forces are zero.")
        
        # Disable other forces
        code_source.comet_mass = 0                            # Disable comet's gravity
        code_source.coeff_disable_tidal = 0                   # Disable tidal effect
        code_source.earth_mass = 0                            # Disable Earth's gravity
        code_source.Q_pr = 0                                  # Disable radiation pressure
        code_source.drag_coeff = 0                            # Disable gas drag force
        code_source.coeff_disable_poynting = 0                # Disable Poynting-Robertson effect
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        # Call the solar gravity validation function
        solar_gravity_validation()
    
    elif validation == "radiation pressure":
        print("Validation of the radiation pressure force implementation. All other forces are zero.")    
        
        # Disable other forces
        code_source.comet_mass = 0                            # Disable comet's gravity
        code_source.earth_mass = 0                            # Disable Earth's gravity
        code_source.solar_mass = 0                            # Disable Solar gravity
        code_source.drag_coeff = 0                            # Disable gas drag force
        code_source.coeff_disable_poynting = 0                # Disable Poynting-Robertson effect
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        # Call the radiation pressure validation function
        radiation_pressure_validation()
        
    elif validation == "poynting":
        print("Validation of the P-R force implementation. All other forces are zero.")
        
        # Disable other forces
        code_source.comet_mass = 0                            # Disable comet's gravity
        code_source.coeff_disable_tidal = 0                   # Disable tidal effect
        code_source.earth_mass = 0                            # Disable Earth's gravity
        code_source.drag_coeff = 0                            # Disable gas drag force
        code_source.coeff_disable_solar = 0                   # Disable Solar gravity
        code_source.Q_pr = 0                                  # Disable radiation pressure
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        
        # Call the Poynting-Robertson validation function
        poynting_robertson_validation()
        
    elif validation == "drag":
        print("Validation of the drag force implementation. All other forces are zero.")
        
        # Disable other forces
        code_source.comet_mass = 0                            # Disable comet's gravity
        code_source.coeff_disable_tidal = 0                   # Disable tidal effect
        code_source.earth_mass = 0                            # Disable Earth's gravity
        code_source.solar_mass = 0                            # Disable Solar gravity
        code_source.coeff_disable_poynting = 0                # Disable Poynting-Robertson effect
        code_source.Q_pr = 0                                  # Disable radiation pressure
        code_source.magnetic_field = np.array([0, 0, 0])      # Disable magnetic field (T)
        code_source.electric_field = np.array([0, 0, 0])      # Disable electric field (V/m)
        
        # Call the gas drag validation function
        gas_drag_validation()      
        
    else:
        print("Error: unrecognized environment. Please enter one of the options.")
        sys.exit()
else : 
    pass
