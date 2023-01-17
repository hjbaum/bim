import numpy as np

from bim_helper import bim
from bim_helper import Params


global k, omega
omega = 2 * np.pi
k = omega/3e8



# Define main function
def main():
    I = 10
    J = 10
    h = 0.25 / I
    bim(Params(1e9, I, J, h))
    

if __name__ == '__main__':
    main()  