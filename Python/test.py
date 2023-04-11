import numpy as np
import csv
import bim_helper as bim

freq = 1e9
fname = "C:/Users/hanzf/OneDrive/Documents/GitHub/microwave_imaging/bim/COMSOL/pt_2.csv"

antenna_coords_fname = "C:/Users/hanzf/OneDrive/Documents/GitHub/microwave_imaging/bim/COMSOL/pt_2.csv"

# Define main function
def main():

    I = 10
    J = 10

    N = 0
    field_data = []

    # Parse the Ez data file line by line
    with open(fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        count = 0       
        for row in reader:
            i = i + 1
            if i == 5:
                N = int(row[1])
            ez = 0
            if i >= 10 : # Data row
                count = count + 1 # not to exceed N
                if row[2] != 'NaN':
                    val = row[2].replace('i', 'j')
                    ez = complex(val)
                else: 
                    ez = 0
                new_node = bim.Node(count, row[0], row[1], ez)
                field_data.append(new_node)
                print(new_node.Ez)              

        csv_file.close()

    antennas = []

    # These are the grid coordinates used to define the E field at each point
    # Parse the antenna coordinates file line by line
    with open(antenna_coords_fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        count = -1      
        for row in reader:
            if count > -1:
                antennas.append(bim.Antenna(count, row[0], row[1]))
            count = count + 1 # not to exceed N
         
        csv_file.close()

    print(N)

    #if N != I*J:
    #    print("N mismatch")

    [first_x, first_y] = field_data[0].coord
    [last_x, last_y] = field_data[N-1].coord 
    
    #dx = (float(last_x)-float(first_x)) / (I-1)
    #dy = (float(last_y)-float(first_y)) / (J-1)

    #print([dx,dy])

    #for multiple sources/freqs, iterate here? Could grab data from different pages in excel
    # ^No. Each independent measurement improves accuracy of the solution
    parameters = bim.Params(freq,0,0,0,0) 
    solution = bim.run(parameters, field_data, antennas)

    # plot solution
    print(solution)


if __name__ == '__main__':
    main()  