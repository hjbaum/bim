import numpy as np
import scipy as sci
import numpy.matlib
import matplotlib
matplotlib.use("TkAgg")
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
    dx = 0.0072 # distance between points along x axis
    dy = 0.0072 # distance between points along y axis
    I = 50 # number of points along x axis
    J = 50 # number of points along y axis
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

# Calculates the mean squared error of two arrays
def calculate_mse(actual, predicted):
    actual = np.array(actual, dtype=np.complex_)
    predicted = np.array(predicted, dtype=np.complex_)
    diff = np.subtract(actual, predicted)
    diff_squared = np.square(diff)
    return abs(diff_squared.mean())

# Calculates the relative mean squared error of two arrays
#  as defined by Chew et al
def calculate_relative_mse(actual, predicted):
    actual = np.array(actual, dtype=np.complex_)
    predicted = np.array(predicted, dtype=np.complex_)
    diff = np.subtract(predicted, actual)
    fraction = np.sum(np.square(diff)) / np.sum(np.square(actual))

    return abs(np.sqrt(fraction))

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

    for point in incident_data:
        x_val = float(point.coord[0])
        X.append(x_val)
    return X

def get_y_vector(): 
    global incident_data

    Y = []

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
    global parameters

    # 4
    epsilon = np.ones([parameters.N,1]) * parameters.eps_min

    return epsilon

# params is a struct containing all required parameters (eg. I, J, Epsilon, etc)
def run(params: Params, inc_data: np.array(Node), scat_data: np.array(Node), base_solution: np.array(float)):
    # Assign globals
    global K
    global Es
    global epsilon    
    global incident_data
    global parameters
    global scatter_data

    # Initialize data structures
    incident_data = inc_data
    scatter_data = scat_data
    parameters = params
    
    # Initialize matrices and arrays
    # M x N matrix
    K = np.zeros([parameters.M, parameters.N],dtype=np.complex_)
    # length N column 
    epsilon = np.zeros([parameters.N, 1]) # Permittivity values
    # length M row 
    Es = np.zeros([parameters.M,1])

    # Indicates if we are calculating the solution for 
    #   an empty domain vs the bleed domain
    running_baseline = False

    # Perform first-order born approximation to 
    #   find initial permittivity solution
    if np.any(base_solution):
        epsilon = base_solution
        running_baseline = False
    else:
        epsilon = initial_guess() 
        running_baseline = True 

    base_epsilon = np.zeros([parameters.N, 1],dtype=np.complex_) 

    if not running_baseline:
        base_epsilon = base_solution

    # Define max number of iterations to compute. Otherwise,
    #   the program will run until the error is less than 5%
    MAX_ITER = 5

    # Define arbitrary regularization param
    gamma = 0.05 #E-10 # Recommendation between 1e-10 and 1e-15

    #fig = plt.figure("Error")
    # Perform pre-calculations by storing Greens summations
    greens_by_source = np.zeros([parameters.M,1], dtype=np.complex_)
    for m in range(0,parameters.M): # Loop over each independent measurement (rx)
        # Observation coordinate
        obs_pt = get_receiver_coord(m) 

        # TODO: Check. Calculate and store Greens integral for this M
        sum = 0
        for p in range(0,parameters.N): # Sum over the domain
            # Find coordinates of pixel
            src_pt = get_grid_coord(p)
            gs = greens(obs_pt, src_pt) # Save for reference
            sum = sum + gs * parameters.dx * parameters.dy
        greens_by_source[m] = sum
    
    # Loop until solution approximates measured data
    iteration = 0    
    awaiting_solution = True
    while awaiting_solution:    
        iteration = iteration + 1
        
        if iteration > MAX_ITER:
            break

        print("Iteration " + str(iteration))
        
        # M: # of transmitters x # of receivers
        # Assume all transmitters are capable of receiving
        for m in range(0,parameters.M): # Loop over each independent measurement (rx)

            # This loop goes over each element in N = I * J pixels
            for n in range(0,parameters.N): 
                # _________________________________________________________
                # Step 2: Solve forward scattering problem (using COMSOL)
                # _________________________________________________________
                # Get data for this pixel 
                Ez_n = get_field_data(n) # n,m)

                # Add Kji element to K matrix 
                K[m,n] = Ez_n * parameters.k**2 * greens_by_source[m] # Todo: Check casting error
       
        # _________________________________________________________
        # STEP 3: Solve matrix eq for permittivity coefficients
        # ___________________________________________________________
        # Calculate LHS (Total electric field) using previous permittivity solution to solve forward problem
        # b = K * a
        # Es = get_scatter_data(m)
        Es = K.dot(epsilon) # Todo: Check casting error

        # ----------------------
        # Perform regularization
        # ----------------------
        H = np.matrix(np.matlib.identity(parameters.N))
        H_t = np.transpose(H)
        K_t = np.transpose(K)   

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
        Es_check = K.dot(epsilon-base_epsilon) # TODO: check

        # ____________________________
        #
        # STEP 5: Compare to measured data
        # ____________________________
        # Calculate error: 

        if not running_baseline:

            # Calculate the average percent error [ (experimental - actual) / actual ]
            scatter_measured = get_scatter_measurements()
            err_array = (Es_check.transpose() - scatter_measured) / scatter_measured
            ave_error = abs(np.average(err_array.transpose()))
            percent_err = ave_error * 100 # Error (%)
            print("Error: ", str(percent_err), "%")
        
            # Calculate relative Mean Squared Error per Chew et al
            mse = calculate_relative_mse(scatter_measured, Es_check.transpose())
            print("MSE: ", str(mse*100), "%")

            if percent_err < 5.0:
                print("Error threshold met!")
                awaiting_solution = False
                break

    final = epsilon - base_epsilon
    return final # Permittivity distribution

