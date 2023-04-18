import numpy as np
import numpy.matlib
import scipy as sci
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

global K, Es, epsilon, parameters, incident_data, scatter_data

# BIM Algorithm:
# 1. Solve linear inverse problem for first order using Born approximation
# 2. Solve forward scattering problem at object and observation points ( USING SOLVER )
# 3. Substitute field (SOLVED FORWARD) in 2 into integrand in integral eq 
#       and solve inverse problem (WITH_REGULARIZATION) to recover next order distribution ftn
# 4. Repeat 2 and compare fields with measured data. If difference < 5% 
#    of scattered field, terminate successfully. Otherwise repeat until soln converges

class Node:
    id = 0
    E = 0.0
    coord = np.array([0.0,0.0])
    def __init__(self, id:int, x : float, y:float,  E:float):
        self.id = id
        self.E = E
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
    M = 12 # number of independent measurements
    eps_rb = 0
    eps_min = 40    # eps_min is 40 in the homogeneous domain
    eps_max = 0
    sig_rb = 0
    sig_min = 0
    sig_max = 0
    misf = 0 # Number of iterations for bim sets; Can be initialized dynamically as well
    
    def __init__(self, freq: int, N : int, M : int):
        self.freq = freq
        self.omega = 2 * np.pi * int(freq)
        self.k = float(self.omega)/3e8
        self.N = N # number of pixels
        self.M = M # number of independent measurements
    
# Globals
K = []
Es = []
epsilon = []
incident_data = []
scatter_data = []

# Remember to subtract background medium (E matrix of constant value)
# Use circle example from Boyu - coordinate.

# Calculate Green's function (if not given)
# r_obj, r_src: coordinate array (x,y,z) for observation and source points
def greens(r_obs, r_src):
    global parameters

    # Scalar distance value
    R = np.sqrt(float((r_obs[0]-r_src[0])**2 + (r_obs[1]-r_src[1])**2))

    gs = 1j/4 * sci.special.hankel1(0,parameters.k * R)
    return gs

# Discretized version of the dyadic Green function
def greens_disc(r_m, r_n):
    global parameters
    
    an = np.sqrt(parameters.dx*parameters.dy/np.pi)
    J1 = sci.special.jv(1,parameters.k*an) # first-order bessel function
    R = np.sqrt((r_m[0]-r_n[0])**2 + (r_m[1]-r_n[1])**2) # distance vector
    H0_2 = sci.special.hankel2(0, parameters.k * R) # 0th-order, 2nd degree hankel function
    green_d = -1j/2 * np.pi * parameters.k * an * J1 * H0_2
    return green_d

def get_field_data(node):
    global incident_data

    # Read data point from csv using the source # and coordinate point
    Ez = incident_data[node].E
    return Ez

def get_scatter_measurements():
    global scatter_data

    E_field = []
    for node in scatter_data:
        E_field.append(node.E)
    return E_field

def get_grid_coord(node): 
    global incident_data

    [x,y] = incident_data[node].coord
    return [float(x),float(y)]

def get_x_vector(): 
    global incident_data

    X = []
    #x_vec = set(X)
    
    # Not correct
    #for point in incident_data:
    #    x_val = float(point.coord[0])
    #    if x_val not in x_vec:
    #        X.append(x_val)
    #        x_vec = set(X)
    #return X

    for point in incident_data:
        x_val = float(point.coord[0])
        X.append(x_val)
    return X

def get_y_vector(): 
    global incident_data

    Y = []
    y_vec = set(Y)
    
    #for point in incident_data:
    #    y_val = float(point.coord[1])
    #    if y_val not in y_vec:
    #        Y.append(y_val)
    #        y_vec = set(Y)
    for point in incident_data:
        x_val = float(point.coord[1])
        Y.append(x_val)
    return Y

# collect antenna locations
def get_receiver_coord(id): 
    global scatter_data

    [x,y] = scatter_data[id].coord
    return [float(x),float(y)]

# TODO: Unit test this function this week; use simple domain
# This is the Born approximation of the first order, 
# as it solves the integral equations assuming the incident
# field is an approximation to the total electric field
def initial_guess(): # "Lobel et al. (1996) - multifrequency" - Batista
    global incident_data
    global parameters
    global K 
    global epsilon
    global Es

    # 2, roughly
    K = np.zeros([parameters.M, parameters.N])
    epsilon = np.zeros([parameters.N, 1]) # Permittivity values
    Es = np.zeros([parameters.M,1])

    # G_sum = 0

    #     for m in range(1,parameters.M):
    #         for n in range(1,parameters.N): 
    #             # Get coordinates of pixel and antenna
    #             src_pt = get_grid_coord(n)
    #             obs_pt = get_source_coord(m)

    #             # Get field for this pixel 
    #             Ez_r = get_field_data(l, n)

    #             # integrate over x and y range
    #             gs = greens(obs_pt, src_pt)

    #             # Add Kji element to K matrix 
    #             K[m,n] = Ez_r * parameters.k**2 * gs

    # # Estimate initial distribution function using zeros vector
    # # TODO: check this
    # guess = np.transpose(K) * Es

    # 4
    epsilon = np.ones([parameters.N,1]) * parameters.eps_min

    ax = plt.axes(projection = '3d')
    X = get_x_vector()
    Y = get_y_vector()
    count = 0

    #for point in incident_data:
    #    ax.scatter(point.coord[0], point.coord[1], epsilon[count])
    #    count = count + 1
    ax.scatter(X, Y, epsilon)

    #ax.plot_surface(X,Y,epsilon)
    plt.show()

    # 5
    #guess = parameters.eps_max

    return epsilon

# params is a struct containing all required parameters (eg. I, J, Epsilon, etc)
def run(params: Params, inc_data: np.array(Node), scat_data: np.array(Node)):
    # Assign globals
    global K
    global Es
    global epsilon    
    global incident_data
    global parameters
    global scatter_data

    incident_data = inc_data
    scatter_data = scat_data
    parameters = params 

    # Perform first-order born approximation to 
    # find initial permittivity solution
    epsilon = initial_guess()
    iteration = 0

    MAX_ITER = 2

    # Loop until solution approximates measured data
    awaiting_solution = True
    while awaiting_solution:    
        iteration = iteration + 1
        # M: # of transmitters x # of receivers
        # Assume all transmitters are capable of receiving
        for m in range(0,parameters.M): # Loop over each independent measurement (rx)
            # Observation coordinate
            obs_pt = get_receiver_coord(m) 

            # This loop goes over each element in N = I * J pixels
            for n in range(0,parameters.N): 

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
                K[m,n] = Ez * parameters.k**2 * gs


        # __________________________________________________________
        #
        # STEP 3: Solve matrix eq for permittivity coefficients
        # ___________________________________________________________
        # Calculate LHS (Total electric field) using previous permittivity solution to solve forward problem
        # b = K * a
        Es = K.dot(epsilon) # Todo: wrong shape

        # eps: column of length N
        # Es: column of length M
        # K : M x N matrix

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
        epsilon = R**-1 * K_t * Es  # TODO: is this correct?

        # ____________________________
        #
        # STEP 4: Solve forward problem to recalculate Ez
        # ____________________________
        # b = K * a
        Es = K.dot(epsilon) # TODO: check

        # ____________________________
        #
        # STEP 5: Compare to measured data
        # ____________________________
        # Calculate error: 
        scatter_measured = get_scatter_measurements()
        err_array = Es.transpose()/scatter_measured
        abs_error = abs(err_array)
        err = np.max(abs_error) * 100 # Error (%)

        # Plot error vs iteration
        plt.scatter(iteration, err)
        if iteration > MAX_ITER - 1:
            break

        if err < 5.0:
            awaiting_solution = False
            break

    plt.show()
    
    return epsilon # Permittivity distribution

