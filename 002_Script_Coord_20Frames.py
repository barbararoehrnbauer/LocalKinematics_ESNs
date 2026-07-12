# -*- coding: mbcs -*-
'''This code opens an Abaqus output database (ODB) file, extracts nodal coordinates from 21 frames
and writes these coordinates to a text file.'''
from abaqus import *
from abaqusConstants import *
from odbAccess import *

files = ['2024-05-26_R_200x6xlsm5_plast_viso0_91_UniAx_1.odb',
'2024-05-27_R_200x6xlsm5_plast_viso0_92_UniAx_2.odb']

for file in files:
    # Open the ODB file
    odb = openOdb(path=file)

    # Define the frames of interest
    frame_numbers = range(21)

    # Open the text file to write the node positions
    with open(str(file[:-4]) + '_COORD.txt', 'w') as file:
        # Write the header
        header = 'Part_Instance, Node_ID, x_ref, y_ref'
        for i in frame_numbers:
            header += ', Frame_{}_X, Frame_{}_Y'.format(i, i)
        file.write(header + '\n')
        
        # Access the step (assuming the first step)
        step = odb.steps[odb.steps.keys()[0]]

        # Access the desired frames
        coord_fields = [step.frames[i].fieldOutputs['COORD'] for i in frame_numbers]
        
        # Iterate through each part instance in the model
        for instanceName, instance in odb.rootAssembly.instances.items():
            # Get the COORD values for the specific instance for each frame
            instance_coords = [coord_fields[i].getSubset(region=instance, position=NODAL) for i in frame_numbers]

            # Create a dictionary to store coordinates for each frame
            coord_dicts = [{value.nodeLabel: value.data for value in instance_coords[i].values} for i in frame_numbers]
            
            # Iterate through each node in the instance
            for node in instance.nodes:
                # Get the node label
                label = node.label

                # Get the coordinates for reference frame (frame 0)
                x_ref, y_ref = coord_dicts[0][label][:2]
                
                # Initialize the line with the reference frame coordinates
                line = '{0}, {1}, {2}, {3}'.format(instanceName, label, x_ref, y_ref)

                # Add the coordinates for each frame to the line
                for i in frame_numbers:
                    x, y = coord_dicts[i][label][:2]
                    line += ', {0}, {1}'.format(x, y)

                # Write the node information to the file
                file.write(line + '\n')

    # Close the ODB file
    odb.close()
