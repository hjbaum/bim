import numpy as np

class Params:

    # Constants
    eps_naught = 8.854187817E-12 #const

    def __init__(self, freq, I, J, h, M):
        self.freq = freq
        self.omega = 2 * np.pi * freq
        self.k = self.omega/3e8
        self.I = I # number of points along x axis
        self.J = J # number of points along y axis
        self.h = h # space between grid points
        self.M = M

# Calculate Green's function
# r_obj, r_src: coordinate array (x,y,z) for observation and source points
def greens(r_obs, r_src, k):
    # Scalar distance value
    R = np.sqrt((r_obs(1)-r_src(1))^2 + (r_obs(2)-r_src(2))^2 + (r_obs(3)-r_src(3))^2)
    gs = np.exp(-1j * k * R) / (4* np.pi * R)
    return gs


# BIM Algorithm:
# 1. Solve linear inverse problem using Born
    # Calculate Green's function

# 2. Solve scattering problem at object and observation points

# 3. a) Substitute field in 2 into integrand in integral eq and 
#    b) Solve inverse problem to recover next order distribution ftn

# 4. Repeat 2 and compare fields with measured data. If difference < 5% 
#    of scattered field, terminate successfully. Otherwise repeat until soln converges


# params is a struct containing all required parameters (eg. I, J, Epsilon, etc)
def bim(params: Params):

    # Loop over I and J
    I = np.array(0,1,params.I)
    J = np.array(0,1,params.J)

    N = I * J # Define N
    # Define M (# of sampled measurements - 1 freq * # of RX)

    for i in I: 
        for j in J:

            # Get the coordinate of this point
            obs_coord = [0,0,0]
            src_coord = [0,0,0]
            gs = greens(obs_coord, src_coord, params.k)

            # Balculate the A matrix using Greens
            Aij = []


            # Define basis functions for permittivity profile - this is what we will solve for
            # pulse: d = sum(a*delta) 

            # Calculate the local scattered electric field
            # Es_ij = 1j * params.k * params.eta * Aij - 1j * params.eta / params.k * Bij

            # Calculate RHS


            # Perform regularization
            H = np.matlib.identity(N)

            # Define arbitrary regularization param
            gamma = 1E-12 # Adjust if needed