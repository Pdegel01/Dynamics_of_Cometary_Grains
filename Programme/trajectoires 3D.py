import matplotlib.pyplot as plt       # Importing the plotting module
import plotly.graph_objects as go
import plotly.offline as pyo
import numpy as np

import inputs
import code_source
from constants import *
from inputs import *

from code_source import *



'''REPRESENTATION 3D DES TRAJECTOIRES '''


if environment == "comet":
    
    '''Calcul des positions et vitesses avec boris pusher'''
    
    
    # Liste pour suivre les positions des différents grains
    all_positions = []
    

    # Conditions initiales pour les grains
    for conditions in initial_conditions:
        positions = []
        velocites = []
        position = conditions['position']
        velocity = conditions['velocity']
        
        positions.append(position)
        all_positions.append({'positions': positions,'velocity': velocity,'surface_potential': initial_potential})


    # Ajouter de nouveaux grains à chaque étape temporelle pour simuler la queue de la comète
    for step in range(num_steps_boris):
        
        # Ajuster la fréquence d'ajout des nouveaux grains si nécessaire
        if step % 100 == 0:  
            print(f"Adding new grains at step {step} of {num_steps_boris}")

            # Générer de nouvelles conditions initiales pour les nouveaux grains
            # Ajuster le nombre de nouveaux grains ajoutés à chaque étape
            for _ in range(10):  
                r = np.random.uniform(comet_radius, (comet_radius + 10))
                theta = np.random.uniform(0, np.pi)
                phi = np.random.uniform(0, 2 * np.pi)
                position = comet_position + np.array([
                    r * np.sin(theta) * np.cos(phi),
                    r * np.sin(theta) * np.sin(phi),
                    r * np.cos(theta)
                ])
                velocity = np.random.uniform(-1e3, 1e3, 3)
                
                all_positions.append({'positions': [position],'velocity': velocity,'surface_potential': initial_potential})

        # Mettre à jour la position de la comète
        comet_position += comet_velocity * dt_boris
        
        # Mettre à jour la position solaire par rapport à la comète en mouvement
        solar_position = comet_position + distance_sun * unit_vector_x

        # Mettre à jour toutes les positions et vitesses des grains
        for grain in all_positions:
            
            positions = grain['positions']
            position = positions[-1]
            velocity = grain['velocity']
            potential = initial_potential
            charge = initial_charge

            # Mettre à jour la position, la vitesse et le potentiel de surface en utilisant l'algorithme de Boris
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)

            # Vérifier la collision avec la comète
            if np.linalg.norm(position - comet_position) <= comet_radius:
                # Immobiliser le grain en cas de collision
                velocity = np.zeros(3)
            
            # Suivre la position relative de la particule par rapport à la comète
            relative_position = position - comet_position
            positions.append(relative_position.copy())



    '''Plot des trajectoires en 3D avec Plotly'''
    

    # Convertir les positions dans un format adapté à Plotly
    traces = []
    
    for grain in all_positions:
        
        positions = grain['positions']
        traces.append(go.Scatter3d(
            x=[pos[0] for pos in positions], 
            y=[pos[1] for pos in positions], 
            z=[pos[2] for pos in positions], 
            mode='lines', 
            line=dict(width=2)
        ))

    # Générer la position actuelle de la sphère de la comète
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = comet_radius * np.cos(u) * np.sin(v)
    y = comet_radius * np.sin(u) * np.sin(v)
    z = comet_radius * np.cos(v)

    # Définir la couleur de la comète en gris foncé
    colorscale = [[0, 'rgb(105, 105, 105)'], [1, 'rgb(105, 105, 105)']]
    comet_sphere = go.Surface(x=x, y=y, z=z, colorscale=colorscale, showscale=False)

    # Définir la mise en page pour le graphique avec un fond noir
    layout = go.Layout(
        title='Particle Trajectories Near a Comet',
        scene=dict(
            xaxis=dict(title='X (m)', showgrid=False, zeroline=False, showbackground=False, backgroundcolor='black'),
            yaxis=dict(title='Y (m)', showgrid=False, zeroline=False, showbackground=False, backgroundcolor='black'),
            zaxis=dict(title='Z (m)', showgrid=False, zeroline=False, showbackground=False, backgroundcolor='black'),
            aspectmode='data',
            bgcolor='black'
        ),
        paper_bgcolor='black',
        plot_bgcolor='black',
        font=dict(color='white'),
        showlegend=False
    )

    # Créer la figure avec les données et la mise en page définies
    fig = go.Figure(data=[comet_sphere] + traces, layout=layout)

    # Afficher la figure
    fig.show()

    # Enregistrer la figure sous forme de fichier HTML et l'ouvrir dans le navigateur
    pyo.plot(fig, filename='particle_trajectories.html')
    
    
        
else:
    
    # Initialization of the counter to track the simulation progress
    counter = 0

    # Lists to store positions and velocities for all debris
    all_positions = []
    all_velocities = []
    all_charges = []
    
    # Earth's position at the origin
    earth_position = np.array([0, 0, 0])
    
    # Loop over each set of initial conditions (all of the debris)
    for conditions in initial_conditions:
        
        # Initialize lists to store positions and velocities for the current debris
        positions = []
        velocities = []
        charges = [initial_charge]
        potentials = [initial_potential]
        
        # Set initial position and velocity
        position = conditions['position']
        velocity = conditions['velocity']
        charge = initial_charge
        potential = initial_potential
        # Loop over the number of steps for the Boris algorithm
        for _ in range(num_steps_boris):
            
            # Increment the counter and print progress every 100 steps
            counter += 1
            print(f"Step {counter} of {num_steps_boris * num_debris}")
                
            # Update position, velocity, and surface potential using the Boris algorithm
            position, velocity, potential, charge = code_source.boris_push(position, velocity, magnetic_field, electric_field, charge, grain_mass, dt_boris, environment, potential)
                
            # Check for collision with the Earth and immobilize the grain upon collision
            if np.linalg.norm(position - earth_position) <= earth_radius:
                velocity = np.zeros(3)
                
            # Record the position and velocity norm for the current step
            positions.append(position.copy())
            velocities.append(velocity.copy())
            charges.append(charge)
            
        # Store all velocities and positions for the current initial condition
        all_velocities.append(np.array(velocities))
        all_positions.append(np.array(positions))
        all_charges.append(charges)
    
    # Convert positions to a format suitable for Plotly
    traces = []
    
    # Define colors for each trajectory
    colors = ['rgb(31, 119, 180)', 'rgb(255, 127, 14)', 'rgb(44, 160, 44)', 'rgb(214, 39, 40)', 'rgb(148, 103, 189)']
    
    # Create a trace for each particle trajectory
    for i, positions in enumerate(all_positions):
        trace = go.Scatter3d(
            x=positions[:, 0], 
            y=positions[:, 1], 
            z=positions[:, 2], 
            mode='lines', 
            line=dict(width=2, color=colors[i % len(colors)])
        )
        
        # Add a sphere at the starting point of the trajectory
        trace['marker'] = dict(
            size=10,
            symbol='circle',
            color=colors[i % len(colors)],
            opacity=0.7,
            line=dict(color='black', width=2)
        )
        
        # Append the trace to the list of traces
        traces.append(trace)
    
    # Generate coordinates for the sphere representing the Earth
    phi, theta = np.mgrid[0.0:2.0*np.pi:50j, 0.0:np.pi:50j]
    x = earth_radius * np.sin(theta) * np.cos(phi)
    y = earth_radius * np.sin(theta) * np.sin(phi)
    z = earth_radius * np.cos(theta)
    
    # Center the Earth sphere at the origin
    earth_center = np.array([0, 0, 0])
    
    # Create the Plotly figure
    fig = go.Figure()
    
    # Add the particle trajectories to the Plotly figure
    for trace in traces:
        fig.add_trace(trace)
    
    # Add the sphere representing the Earth to the Plotly figure
    fig.add_trace(go.Surface(
        x=x + earth_center[0],  
        y=y + earth_center[1],  
        z=z + earth_center[2],  
        colorscale='Blues',  
        showscale=False, 
        opacity=0.5,  
        name="Earth"  
    ))
    
    # Parameters text for display, legend
    params_text = f"""
    Simulation Parameters:
        Grain Mass: {grain_mass} kg
        Electric Field: {electric_field} V/m
        Magnetic Field: {magnetic_field} T
        Time Step: {dt_boris} s
        Number of Steps: {num_steps_boris}
    """
    
    # Define layout for the plot
    layout = go.Layout(
        title='Particle Trajectories Near the Earth',
        scene=dict(
            xaxis=dict(title='X (m)', showgrid=False, zeroline=False, showbackground=False),
            yaxis=dict(title='Y (m)', showgrid=False, zeroline=False, showbackground=False, showticklabels=False),
            zaxis=dict(title='Z (m)', showgrid=False, zeroline=False, showbackground=False, showticklabels=False),
            aspectmode='data'
        ),
        showlegend=False,
        annotations=[
            dict(
                text=params_text,
                x=1.05,
                y=1.05,
                xref='paper',
                yref='paper',
                showarrow=False,
                font=dict(size=12)
            )
        ]
    )
    
    # Update the layout of the figure
    fig.update_layout(layout)
    
    # Show the figure
    fig.show()
    
    # Save the figure as an HTML file and open it in the browser
    pyo.plot(fig, filename='particle_trajectories.html', auto_open=True)
    
