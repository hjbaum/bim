import numpy as np
import csv
import bim_helper as bim
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

freq = 1e9
scatter_field_fname = "C:/Users/hanzf/OneDrive/Documents/GitHub/microwave_imaging/bim/COMSOL/Homogeneous/pt2_scatter.csv"
incident_field_fname = "C:/Users/hanzf/OneDrive/Documents/GitHub/microwave_imaging/bim/COMSOL/Homogeneous/pt2_field_regular.csv"
antenna_coords_fname = "C:/Users/hanzf/OneDrive/Documents/GitHub/microwave_imaging/bim/COMSOL/Homogeneous/AntennaLocations.csv"

# Define main function
def main():

    N = 0
    x_num = 20
    y_num = 20

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
    solution = bim.run(parameters, incident_field_data, scatter_field_data)

    # plot solution
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    X = bim.get_x_vector()
    Y = bim.get_y_vector()
    ax.scatter(X, Y, solution)
    #ax.plot_surface(X,Y, solution)
    plt.show()


if __name__ == '__main__':
    main()  