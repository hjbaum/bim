import numpy as np
import scipy as sci

global K, Es, epsilon
global parameters, node_data

# BIM Algorithm:
# 1. Solve linear inverse problem for first order using Born approximation
# 2. Solve forward scattering problem at object and observation points ( USING SOLVER )
# 3. Substitute field (SOLVED FORWARD) in 2 into integrand in integral eq 
#       and solve inverse problem (WITH_REGULARIZATION) to recover next order distribution ftn
# 4. Repeat 2 and compare fields with measured data. If difference < 5% 
#    of scattered field, terminate successfully. Otherwise repeat until soln converges

class Node:
    id = 0
    Ez = 0.0
    coord = np.array([0.0,0.0])
    def __init__(self, id:int, x : float, y:float,  Ez:float):
        self.id = id
        self.Ez = Ez
        self.coord = np.array([float(x),float(y)])
    
class Params:
    # Constants
    eps_naught = 8.854187817E-12 #const
    freq = 0.0
    omega = 0.0
    k = 0
    dx = 0.0 # distance between points along x axis
    dy = 0.0 # distance between points along y axis
    I = 0 # number of points along x axis
    J = 0 # number of points along y axis
    N = 0 # number of pixels
    M = 1 # number of independent measurements
    L = 1 # number of sources
    eps_rb = 0
    eps_min = 0
    eps_max = 0
    sig_rb = 0
    sig_min = 0
    sig_max = 0
    misf = 0 # Number of iterations for bim sets; Can be initialized dynamically as well
    
    def __init__(self, freq: int, dx : float, dy : float, I : int, J : int):
        self.freq = freq
        self.omega = 2 * np.pi * freq
        self.k = self.omega/3e8
        self.dx = dx # distance between points along x axis
        self.dy = dy # distance between points along y axis
        self.I = I  # number of points along x axis
        self.J = J # number of points along y axis
        self.N = int(I*J) # number of pixels

    #def __init__(self, freq, dx, dy, I, J, L, M, eps_rb, eps_min, eps_max, sig_rb, sig_min, sig_max):
    #    self.freq = freq
    #    self.omega = 2 * np.pi * freq
    #    self.k = self.omega/3e8
    #    self.dx = dx # distance between points along x axis
    #    self.dy = dy # distance between points along y axis
    #    self.I = I # number of points along x axis
    #    self.J = J # number of points along y axis
    #    self.N = I*J # number of pixels
    #    self.M = M # number of independent measurements
    #    self.L = L # number of sources
    #    self.eps_rb = eps_rb
    #    self.eps_min = eps_min
    #    self.eps_max = eps_max
    #    self.sig_rb = sig_rb
    #    self.sig_min = sig_min
    #    self.sig_max = sig_max
    #    self.misf = 5 # Number of iterations for bim sets; Can be initialized dynamically as well
    
# Globals
K = []
Es = []
epsilon = []
node_data = []


# Calculate Green's function (if not given)
# r_obj, r_src: coordinate array (x,y,z) for observation and source points
def greens(r_obs, r_src):
    global parameters

    # Scalar distance value
    R = np.sqrt((r_obs(1)-r_src(1))^2 + (r_obs(2)-r_src(2))^2 + (r_obs(3)-r_src(3))^2)
    gs = 1j/4 * sci.hankel1(0,parameters.k * R)
    return gs

# Discretized version of the dyadic Green function
def greens_disc(r_m, r_n):
    global parameters
    
    an = np.sqrt(parameters.dx*parameters.dy/np.pi)
    J1 = sci.jv(1,parameters.k*an) # first-order bessel function
    R = np.sqrt((r_m[0]-r_n[0])^2 + (r_m[1]-r_n[1])^2) # distance vector
    H0_2 = sci.hankel2(0, parameters.k * R) # 0th-order, 2nd degree hankel function
    green_d = -1j/2 * np.pi * parameters.k * an * J1 * H0_2
    return green_d

def get_field_data(node):
    global node_data

    # Read data point from csv using the source # and coordinate point
    Ez = node_data[node].Ez
    return Ez


def get_complete_field_data():
    global node_data

    field_data = []
    for node in node_data:
        field_data.append(node.Ez)

    return field_data


def get_grid_coord(node): 
    global node_data

    [x,y] = node_data[node].coord
    return [float(x),float(y)]

# TODO: collect antenna locations
def get_source_coord(id): 
    # These should be pre-mapped
    x = 1
    y = 1
    return np.array[x,y]

# TODO: what is L vs M? M: # independent measurements, L: # sources?
# TODO: Fix this
# TODO: Unit test this function this week; use simple domain
# This is the Born approximation of the first order, 
# as it solves the integral equations assuming the incident
# field is an approximation to the total electric field
def initial_guess(): # "Lobel et al. (1996) - multifrequency" - Batista
    global node_data
    global parameters
    global K
    global epsilon
    global Es

    # Estimate permittivity using the min value provided
    K = np.zeros([parameters.M * parameters.L, parameters.N])
    epsilon = np.zeros([parameters.N, 1]) # Permittivity values
    Es = np.zeros([parameters.M*parameters.L,1])

    G_sum = 0
    for l in range(1,parameters.L):
        for m in range(1,parameters.M):
            for n in range(1,parameters.N): 
                # Get coordinates of pixel and antenna
                src_pt = get_grid_coord(n)
                obs_pt = get_source_coord(m)

                # Get field for this pixel 
                Ez_r = get_field_data(l, n)

                # integrate over x and y range
                gs = greens(obs_pt, src_pt)

                # Add Kji element to K matrix 
                K[m,n] = Ez_r * parameters.k^2 * gs

    # Estimate initial distribution function using zeros vector
    # TODO: check this
    epsilon = np.transpose(K) * Es

# params is a struct containing all required parameters (eg. I, J, Epsilon, etc)
def run(params: Params, data: np.array(Node)):
    # Assign globals

    global K
    global Es
    global epsilon    
    global node_data
    global parameters

    node_data = data
    parameters = params 

    # Perform first-order born approximation to 
    # find initial permittivity solution
    epsilon = initial_guess()

    # Loop until solution approximates measured data
    awaiting_solution = True
    while awaiting_solution:    
        # M: # of transmitters x # of receivers
        # Assume all transmitters are capable of receiving
        for l in range(1,parameters.L): # Loop over each source (tx)
            for m in range(1,parameters.M): # Loop over each independent measurement (rx)

                # Observation coordinate
                obs_pt = get_source_coord(m) # We will have one tx at a time and 11 rx -> M = 11

                # This loop goes over each element in N = I * J pixels
                for n in range(1,parameters.N): 

                    # This may not be necessary. Depends on how data is arranged
                    # Source coordinate
                    # Find coordinates of pixel
                    src_pt = get_grid_coord(n)

                    # _________________________________________________________
                    #
                    # Step 2: Solve forward scattering problem (using COMSOL)
                    # _________________________________________________________
                    # Get data for this pixel 
                    Ez = get_field_data(n)

                    # TODO: check this
                    # Perform discretized surface integral over x and y range
                    # K_sum = 0                    
                    # for x in range(parameters.dx/2, parameters.I - parameters.dx/2, parameters.dx): # Use centers of each pixel
                    #     for y in range(parameters.dy/2, parameters.J - parameters.dy/2, parameters.dy):
                    #         coord = [x,y]
                    #         K_sum = K_sum + greens(obs_pt, coord)
                
                    # TODO: integrate over x and y range
                    gs = greens(obs_pt, src_pt)

                    # Add Kji element to K matrix 
                    K[m,n] = Ez * parameters.k^2 * gs


        # __________________________________________________________
        #
        # STEP 3: Solve matrix eq for permittivity coefficients
        # ___________________________________________________________
        # Calculate LHS (Total electric field) using previous permittivity solution to solve forward problem
        # b = K * a
        Es = K * epsilon

        # ----------------------
        # Perform regularization
        # ----------------------
        H = np.matrix(np.matlib.identity(parameters.N))
        H_t = np.transpose(H)
        K_t = np.transpose(K)   

        # Define arbitrary regularization param
        gamma = 1E-12 # Should be btn 1e-10 and 1e-15

        # Regularization matrix
        R = K_t.dot(K) + gamma * H_t.dot(H)

        # ------------------------------
        # Solve for permittivity profile
        # ------------------------------
        epsilon = R^-1 * K_t * Es

        # ____________________________
        #
        # STEP 4: Solve forward problem to recalculate Ez_s
        # ____________________________
        # 
        Es = K * epsilon
                
        # ____________________________
        #
        # STEP 5: Compare to measured data
        # ____________________________
        # Calculate error: 

        # TODO: fix this. need to compare all data at once, or piece by piece
        err = abs((Es - get_complete_field_data())/Es) * 100 # Error (%)

        if err < 5.0:
            awaiting_solution = False
            break
    
    return epsilon # Permittivity distribution

