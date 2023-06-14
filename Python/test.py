import numpy as np
import csv
import bim_helper as bim
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import os

freq = 0.9e9

def solve_baseline(incident_field_fname, scatter_field_fname, tx_id):
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
                    count = count + 1 # not to exceed N
                    val = row[2].replace('i', 'j')
                    ez = complex(val)
                    new_node = bim.Node(count, row[0], row[1], ez)
                    incident_field_data.append(new_node)
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
                    count = count + 1 
                    val = row[2].replace('i', 'j')
                    e_field = complex(val)
                    new_node = bim.Node(count, row[0], row[1], e_field)
                    scatter_field_data.append(new_node)
        csv_file.close()
        M = count

    #for multiple sources/freqs, iterate here? Could grab data from different pages in excel
    # ^No. Each independent measurement improves accuracy of the solution
    parameters = bim.Params(freq,N,M) 
    solution = bim.run(parameters, incident_field_data, scatter_field_data, [], tx_id)
    
    return solution

def solve_bleed(bleed_incident_field_fname, bleed_scatter_field_fname, baseline_solution, tx_id):
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
                    count = count + 1 # not to exceed N
                    val = row[2].replace('i', 'j')
                    ez = complex(val)
                    new_node = bim.Node(count, row[0], row[1], ez)
                    incident_field_data.append(new_node)
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
                    count = count + 1 
                    val = row[2].replace('i', 'j')
                    e_field = complex(val)
                    new_node = bim.Node(count, row[0], row[1], e_field)
                    scatter_field_data.append(new_node)
        csv_file.close()
        M = count

    #for multiple sources/freqs, iterate here? Could grab data from different pages in excel
    # ^No. Each independent measurement improves accuracy of the solution
    parameters = bim.Params(freq,N,M) 
    solution = bim.run(parameters, incident_field_data, scatter_field_data, baseline_solution, tx_id)

    return solution

# Define main function
def main():
    dirname = os.path.dirname(__file__)
    files_by_source_csv = os.path.join(dirname, 'filenames_by_source.csv')

    print(dirname)

    fig = plt.figure("Solution")
    ax = plt.axes(projection = '3d')
    colors = ['blue', 'green']
    cmap_temp = ['hot', 'winter']

    # Parse the csv for the source files containing field data
    with open(files_by_source_csv,'r') as csv_file:
        reader = csv.reader(csv_file)    
        row_count = 0   

        
        for row in reader:
            row_count = row_count + 1
            if row_count > 1:

                source_id = row[0]  

                # Row format:
                # Col 1     Col 2         Col 3                   Col 4                     Col 5                      Col 6
                # Source #  Folder name   Scatter file for bleed  Incident file for bleed   Scatter file for baseline  Incident file for baseline

                # with bleed
                bleed_scatter_field_fname = os.path.join(dirname, '..', 'COMSOL', row[1], row[2])
                bleed_incident_field_fname =  os.path.join(dirname, '..', 'COMSOL', row[1], row[3])
                                                           
                # baseline domain
                scatter_field_fname = os.path.join(dirname, '..', 'COMSOL', row[1], row[4])
                incident_field_fname = os.path.join(dirname, '..', 'COMSOL', row[1], row[5])

                # Solve base field
                baseline_solution = solve_baseline(incident_field_fname, scatter_field_fname, int(source_id))
                # Plot absolute value; currently dropping imag component
                #ax.scatter(X, Y, abs(baseline_solution), label=legend, color=colors[row_count-2], lw=0.5)    
                
                # Solve field of interest
                solution = solve_bleed(bleed_incident_field_fname, bleed_scatter_field_fname, baseline_solution, int(source_id))

                # Plot solution
                legend = 'Source ' + str(source_id)

                X = bim.get_x_vector()
                Y = bim.get_y_vector()
                # Plot absolute value; currently dropping imag component
                ax.scatter(X, Y, abs(solution), label=legend, color=colors[row_count-2], lw=0.5)    
                
                #my_cmap = plt.get_cmap(cmap_temp[row_count-2])
                #surf = ax.plot_surface(X,Y,solution, cmap=my_cmap, edgecolor='none', label=legend)
                #surf._edgecolors2d = surf._edgecolor3d
                #surf._facecolors2d = surf._facecolor3d
                #fig.colorbar(surf, ax=ax)
                
                break
                
        csv_file.close()
    
    ax.legend()  
    plt.show()  
    

if __name__ == '__main__':
    main()  