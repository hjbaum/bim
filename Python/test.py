import numpy as np
import csv
#from bim_helper import bim
#from bim_helper import Params

class Node:
     
    id = 0
    Ez = 0.0
    coord = np.array([0.0,0.0])
    def __init__(self, id:int, x : float, y:float,  Ez:float):
        self.id = id
        self.Ez = Ez
        self.coord = np.array([x,y])


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
                new_node = Node(count, row[0], row[1], ez)
                field_data.append(new_node)
                print(new_node.Ez)              

        csv_file.close()

    print(N)

    # Perform parsing


    if N != I*J:
        print("N mismatch")
    # dx = (first_x + last_x) / I
    # dy = (first_y + last_y) / J

    # Until file end
    #if es == 'NaN':
    #    es = 0

    #Ez = 
    #bim(Params(), fname)   

if __name__ == '__main__':
    main()  