import numpy as np

class Params:

    # Constants
    eps_naught = 8.854187817E-12 #const

    def __init__(self, freq, dx, dy, I, J, L, M, eps_rb, eps_min, eps_max, sig_rb, sig_min, sig_max):
        self.freq = freq
        self.omega = 2 * np.pi * freq
        self.k = self.omega/3e8
        self.dx = dx
        self.dy = dy
        self.I = I # number of points along x axis
        self.J = J # number of points along y axis
        self.N = I*J # number of pixels
        self.M = M # number of independent measurements
        self.L = L # number of sources

        self.eps_rb = eps_rb
        self.eps_min = eps_min
        self.eps_max = eps_max

        self.sig_rb = sig_rb
        self.sig_min = sig_min
        self.sig_max = sig_max

        self.misf = 5 # Number of iterations for bim sets; Can be initialized dynamically as well

# Calculate Green's function (if not given)
# r_obj, r_src: coordinate array (x,y,z) for observation and source points
def greens(r_obs, r_src, k):
    # Scalar distance value
    R = np.sqrt((r_obs(1)-r_src(1))^2 + (r_obs(2)-r_src(2))^2 + (r_obs(3)-r_src(3))^2)
    gs = np.exp(-1j * k * R) / (4* np.pi * R)
    return gs

def greens_d():
    # TODO: Check what this is
    stub = 1

def solve_forward(params: Params):
    # Solve the forward problem
    stub = 1

def solve_inverse(params: Params):
    # Etc.
    stub = 1

global K, b

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
def bim(ei, es, gs, params: Params):

    # Loop over I and J
    I = np.array(0,1,params.I)
    J = np.array(0,1,params.J)
    L = np.array(0,1,params.L)
    N = np.array(0,1,params.N)
    M = np.array(0,1,params.M)

    # Define pulse functions: the same pulse functions are used in both the 
    # forward and inverse procedures,


    # Calculate initial matrix. Can just fill with 


    
    for num_iter in np.array(0,1,params.misf):


 

            gs = greens(obs_coord, src_coord, params.k)


            # Calculate forward solution
            for l in L:
                solve_forward(params)

            # Calculate the K (aka A) matrix using Greens
            # Pseudocode: 
            #   K = [ Re{gs.*es} ...                  |  Im{gs.*es}                      ]
            #       [ Im{gs.*es} / (omega*eps0*eps_rb)| -Re{gs.*es} / (omega*eps0*eps_rb)]


            # Define basis functions for permittivity profile - this is what we will solve for
            # pulse: d = sum(a*delta) 

            # Calculate LHS (Total electric field)
            # Column vector
            # b = [b(index)+ current]

            # Solve matrix equation for permittivity profile
            # Column vector
            # a = A \ B



            # Perform regularization
            H = np.matlib.identity(N)

            # Define arbitrary regularization param
            gamma = 1E-12 # Adjust if needed