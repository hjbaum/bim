import numpy as np
import csv
import bim_helper as bim

freq = 1e9
fname = "C:/Users/hanzf/OneDrive/Documents/GitHub/microwave_imaging/bim/Python/pt_2_sig_2.csv"

# Define main function
def main():

    I = 10
    J = 10

    N = 0
    field_data = []

    # Parse the data file line by line
    with open(fname,'r') as csv_file:
        reader = csv.reader(csv_file)
        i = 0
        for row in reader:
            i = i + 1
            count = 0
            if i == 5:
                N = int(row[1])
            ez = 0
            if i >= 10 : # Data row
                count = count + 1 # not to exceed N
                if row[2] != 'NaN':
                    ez = row[2]
                else: 
                    ez = 0
                new_node = bim.Node(count, row[0], row[1], ez)
                field_data.append(new_node)
                print(new_node.Ez)              

        csv_file.close()

    print(N)

    if N != I*J:
        print("N mismatch")

    [first_x, first_y] = field_data[0].coord
    [last_x, last_y] = field_data[N-1].coord 
    
    dx = (float(last_x)-float(first_x)) / (I-1)
    dy = (float(last_y)-float(first_y)) / (J-1)

    print([dx,dy])

    #for multiple sources/freqs, iterate here? Could grab data from different pages in excel
    # ^No. Each independent measurement improves accuracy of the solution
    parameters = bim.Params(freq,dx,dy,I,J) 
    solution = bim.run(parameters, field_data)

    # plot solution
    print(solution)


if __name__ == '__main__':
    main()  