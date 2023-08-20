import numpy as np
import csv
import bim_helper as bim
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits import mplot3d
import os

# Define the frequency used in the COMSOL siimulation
freq = 0.9e9

def solve_baseline(incident_field_fname, scatter_field_fname):
    incident_field_data = []
    scatter_field_data = []

    # Parse the incident field data file line by line and determine N
    with open(incident_field_fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        count = 0       
        for row in reader:
            i = i + 1
            if i == 5:
                N_total = int(row[1]) # for debug purposes; will likely exceed N
            ez = 0
            if i >= 10 : # Data row
                if row[2] != 'NaN': # if NaN, this point is outside the antenna limits
                    val = row[2].replace('i', 'j')
                    ez = complex(val)
                    new_node = bim.Node(count, row[0], row[1], ez)
                    incident_field_data.append(new_node)
                    count = count + 1 # not to exceed N
        csv_file.close()
        N = count # Check

    # Parse the scatter data file line by line and determine M
    M = 0
    with open(scatter_field_fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        count = 0       
        for row in reader:
            i = i + 1
            ez = 0
            if i >= 10 : # Data row
                if row[2] != 'NaN':
                    val = row[2].replace('i', 'j')
                    e_field = complex(val)
                    new_node = bim.Node(count, row[0], row[1], e_field)
                    scatter_field_data.append(new_node)
                    count = count + 1 
        csv_file.close()
        M = count
        print("M = ", str(M))

    #for multiple sources/freqs, iterate here? Could grab data from different pages in excel
    # ^No. Each independent measurement improves accuracy of the solution
    parameters = bim.Params(freq,N,M) 
    solution = bim.run(parameters, incident_field_data, scatter_field_data, [])
    
    return solution

def solve_bleed(bleed_incident_field_fname, bleed_scatter_field_fname, baseline_solution):
    incident_field_data = []
    scatter_field_data = []

    # Parse the incident field data file line by line and determine N
    with open(bleed_incident_field_fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        count = 0       
        for row in reader:
            i = i + 1
            if i == 5:
                N_total = int(row[1]) # for debug purposes; will likely exceed N
            ez = 0
            if i >= 10 : # Data row
                if row[2] != 'NaN': # if NaN, this point is outside the antenna limits
                    val = row[2].replace('i', 'j')
                    ez = complex(val)
                    new_node = bim.Node(count, row[0], row[1], ez)
                    incident_field_data.append(new_node)
                    count = count + 1 # not to exceed N
        csv_file.close()
        N = count # Check

    # Parse the scatter data file line by line and determine M
    M = 0
    with open(bleed_scatter_field_fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        count = 0       
        for row in reader:
            i = i + 1
            ez = 0
            if i >= 10 : # Data row
                if row[2] != 'NaN':
                    val = row[2].replace('i', 'j')
                    e_field = complex(val)
                    new_node = bim.Node(count, row[0], row[1], e_field)
                    scatter_field_data.append(new_node)
                    count = count + 1         
        csv_file.close()
        M = count
        print("M = ", str(M))

    #for multiple sources/freqs, iterate here? Could grab data from different pages in excel
    # ^No. Each independent measurement improves accuracy of the solution
    parameters = bim.Params(freq,N,M) 
    solution = bim.run(parameters, incident_field_data, scatter_field_data, baseline_solution)

    return solution

def main():
    dirname = os.path.dirname(__file__)
    files_by_source_csv = os.path.join(dirname, 'filenames_by_source.csv')

    print(dirname)

    fig = plt.figure("Solution")
    ax = plt.axes(projection = '3d')
    colors = ['blue', 'green']

    # Parse the csv for the source files containing field data
    with open(files_by_source_csv,'r') as csv_file:
        reader = csv.reader(csv_file)    
        row_count = 0   

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

                # with bleed
                bleed_scatter_field_fname = os.path.join(dirname, '..', 'COMSOL', row[2], row[3])
                bleed_incident_field_fname =  os.path.join(dirname, '..', 'COMSOL', row[2], row[4])
                                                           
                # baseline domain
                scatter_field_fname = os.path.join(dirname, '..', 'COMSOL', row[2], row[5])
                incident_field_fname = os.path.join(dirname, '..', 'COMSOL', row[2], row[6])

                # Solve base field
                baseline_solution = solve_baseline(incident_field_fname, scatter_field_fname)


                # Solve field of interest
                solution = solve_bleed(bleed_incident_field_fname, bleed_scatter_field_fname, baseline_solution)

                X = bim.get_x_vector()
                Y = bim.get_y_vector()

                name = 'epsilon_src' + str(source_id) + '.mat'

                scipy.io.savemat(name, dict(x=X, y=Y, solution=solution))

                # Plot solution
                legend = 'Source ' + str(source_id)

                # Plot absolute value; currently dropping imag component
                ax.scatter(X, Y, abs(solution), label=legend, color=colors[row_count-2], lw=0.5)    
                
        csv_file.close()
    
    ax.legend()  
    plt.show()  
    

if __name__ == '__main__':
    main()  