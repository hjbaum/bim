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

# Questions
# Does get_incident_field alter coarsemodel?
# What's the diff btn coarse and finemodels
# why have finemodel, coarsemodel, and parameters?
# What portion is the Born approximation? 
#   Ez ^(r) is the solution of the forward problem,
#   the forward scattering solution of order r, when r = 0. It is 
#   the incident field in the object (the Born approximation).


# params is a struct containing all required parameters (eg. I, J, Epsilon, etc)
def bim(params: Params):

    # Loop over I and J
    I = np.array(0,1,params.I)
    J = np.array(0,1,params.J)

    # Define pulse functions: the same pulse functions are used in both the 
    # forward and inverse procedures,

    N = I * J # Define N
    # Define M (# of sampled measurements - 1 freq * # of RX)
    #   the number of the independent measurement data
    #   the product of the number of receivers and the number of incident waves or transmitters. 

    for i in I: 
        for j in J:

            
            #rho_i is the coordinate of the center of patch i.

            # Get the coordinate of this point
            obs_coord = [0,0,0]
            src_coord = [0,0,0]
            gs = greens(obs_coord, src_coord, params.k)

            # Calculate the K matrix using Greens
            Kji = []
            # Kji = k^2 * Ez_r(r_i) * double_integral(greens(pj - p') * dx' * dy' )


            # Define basis functions for permittivity profile - this is what we will solve for
            # pulse: d = sum(a*delta) 

            # Calculate LHS (Total electric field)
            # Column vector
            b = []

            # Solve matrix equation for permittivity profile
            # Column vector
            # a = A \ B



            # Perform regularization
            H = np.matlib.identity(N)

            # Define arbitrary regularization param
            gamma = 1E-12 # Adjust if needed