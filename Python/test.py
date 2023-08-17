import numpy as np
import csv
import bim_helper as bim
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits import mplot3d
import os

freq = 0.9e9

def main():
    dirname = os.path.dirname(__file__)
    files_by_source_csv = os.path.join(dirname, 'filenames_by_source.csv')

    # Define arrays
    base_incident_field_data = [] 
    base_scatter_field_data = []
    bleed_incident_field_data = [] 
    bleed_scatter_field_data = []
    x_vector = []
    y_vector = []
    
    # _______________________________
    # Read simulation data from files
    # _______________________________
    # Parse the csv for the source files containing field data
    with open(files_by_source_csv,'r') as csv_file:
        reader = csv.reader(csv_file)    
        row_count = 0 

        # Source count
        m = 0  
        M = 0 # Number of antennas
        N = 0 # Number of pixels

        for row in reader:
            row_count = row_count + 1
            if row_count > 1:

                source_id = row[0]  

                # Row format:
                # Col 0         Col 1           Col 2                       
                # Source ID     Alternate ID    Folder name     
                #
                # Col 3                         Col 4                     
                # Scatter file for bleed        Incident file for bleed   
                # 
                # Col 5                         Col 6
                # Scatter file for baseline     Incident file for baseline

                # baseline domain
                scatter_field_fname = os.path.join(dirname, '..', 'COMSOL', row[2], row[5])
                incident_field_fname = os.path.join(dirname, '..', 'COMSOL', row[2], row[6])        
                               
                # with bleed
                bleed_scatter_field_fname = os.path.join(dirname, '..', 'COMSOL', row[2], row[3])
                bleed_incident_field_fname =  os.path.join(dirname, '..', 'COMSOL', row[2], row[4])

                ## BASELINE DATA
                # Collect point coordinates as we ...
                # Parse the incident field data file line by line and determine N
                data_list = []
                with open(incident_field_fname,'r') as base_inc_csv_file:
                    base_inc_reader = csv.reader(base_inc_csv_file)
                    i = 0
                    count = 0       
                    for row in base_inc_reader:
                        i = i + 1
                        ez = 0
                        if i >= 10 : # Data row
                            if row[2] != 'NaN': # if NaN, this point is outside the antenna limits
                                val = row[2].replace('i', 'j')
                                ez = complex(val)
                                new_node = bim.Node(count, row[0], row[1], ez)
                                if m==0:
                                    x_vector.append(row[0])
                                    y_vector.append(row[1])
                                data_list.append(new_node)
                                count = count + 1 # not to exceed N
                    
                    base_incident_field_data.append(data_list)

                    base_inc_csv_file.close()
                    if m == 0: N = count # Check

                data_list = []
                with open(scatter_field_fname,'r') as base_scat_csv_file:
                    base_scat_reader = csv.reader(base_scat_csv_file)
                    i = 0
                    count = 0       
                    for row in base_scat_reader:
                        i = i + 1
                        if i >= 10 : # Data row
                            if row[2] != 'NaN':
                                val = row[2].replace('i', 'j')
                                e_field = complex(val)
                                new_node = bim.Node(count, row[0], row[1], e_field)
                                data_list.append(new_node)
                                count = count + 1 
                    
                    base_scatter_field_data.append(data_list)

                    base_scat_csv_file.close()
                    if m == 0: M = count

                ## BLEED DATA
                # Parse the incident field data file line by line and determine N
                data_list = []
                with open(bleed_incident_field_fname,'r') as bleed_inc_csv_file:
                    bleed_inc_reader = csv.reader(bleed_inc_csv_file)
                    i = 0
                    count = 0       
                    for row in bleed_inc_reader:
                        i = i + 1
                        if i == 5:
                            N_total = int(row[1]) # for debug purposes; will likely exceed N
                        ez = 0
                        if i >= 10 : # Data row
                            if row[2] != 'NaN': # if NaN, this point is outside the antenna limits
                                val = row[2].replace('i', 'j')
                                ez = complex(val)
                                new_node = bim.Node(count, row[0], row[1], ez)
                                data_list.append(new_node)
                                count = count + 1 # not to exceed N
                    
                    bleed_incident_field_data.append(data_list)

                    bleed_inc_csv_file.close()

                # Parse the scatter data file line by line
                data_list = []
                with open(bleed_scatter_field_fname,'r') as bleed_scat_csv_file:
                    bleed_scat_reader = csv.reader(bleed_scat_csv_file)
                    i = 0
                    count = 0       
                    for row in bleed_scat_reader:
                        i = i + 1
                        ez = 0
                        if i >= 10 : # Data row
                            if row[2] != 'NaN':
                                val = row[2].replace('i', 'j')
                                e_field = complex(val)
                                new_node = bim.Node(count, row[0], row[1], e_field)
                                data_list.append(new_node)
                                count = count + 1 
                    
                    bleed_scatter_field_data.append(data_list)
        
                    bleed_scat_csv_file.close()
                                                
                m = m + 1

        csv_file.close()
   
    # _____________________________________
    # Solve using the Born Iterative Method 
    # _____________________________________
    # Define parameters
    parameters = bim.Params(freq,N,M) 
    
    # Solve base field
    baseline_solution = bim.run(parameters, base_incident_field_data, base_scatter_field_data, [])
    
    # Solve field of interest
    solution = bim.run(parameters, bleed_incident_field_data, bleed_scatter_field_data, baseline_solution)

    # _______________________
    # Perform post-processing
    # _______________________
    X = np.array(x_vector)
    Y = np.array(y_vector)

    # Plot the absolute value
    fig = plt.figure("Solution")
    ax = plt.axes(projection = '3d')
    plt.axis([X[0], X[N-1], Y[0], Y[N-1]])
    ax.scatter(X, Y, abs(np.array(solution)), color='red')   

    # Save data to a Matlab data file
    name = 'epsilon_src' + str(source_id) + '.mat'
    scipy.io.savemat(name, dict(x=X, y=Y, solution=solution))

    plt.show()  
    

if __name__ == '__main__':
    main()  