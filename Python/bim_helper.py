import numpy as np
import scipy as sci

global K, b, a

global p; Xn

class Params:

    # Constants
    eps_naught = 8.854187817E-12 #const

    def __init__(self, freq, dx, dy, I, J, L, M, eps_rb, eps_min, eps_max, sig_rb, sig_min, sig_max):
        self.freq = freq
        self.omega = 2 * np.pi * freq
        self.k = self.omega/3e8

        self.dx = dx # distance between points along x axis
        self.dy = dy # distance between points along y axis
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
def greens(r_obs, r_src):
    # Scalar distance value
    R = np.sqrt((r_obs(1)-r_src(1))^2 + (r_obs(2)-r_src(2))^2 + (r_obs(3)-r_src(3))^2)
    gs = 1j/4 * sci.hankel1(0,p.k * R)
    return gs

# Discretized version of the dyadic Green function
def greens_disc(r_m, r_n):
    an = np.sqrt(p.dx*p.dy/np.pi)
    J1 = sci.jv(1,p.k*an) # first-order bessel function
    R = np.sqrt((r_m[0]-r_n[0])^2 + (r_m[1]-r_n[1])^2) # distance vector
    H0_2 = sci.hankel2(0, p.k * R) # 0th-order, 2nd degree hankel function
    green_d = -1j/2 * np.pi * p.k * an * J1 * H0_2
    return green_d


def get_field_data(source_num, element):
    # Read data point from csv using the source # and coordinate point
    Ei = 1
    return Ei

def get_coord(N): 
    row = N / p.J  
    col = N % p.J  
    return np.array[row,col]

def get_antenna_coord(id): 
    # These should be pre-mapped
    x = 1
    y = 1

    return np.array[x,y]

# TODO: Fix this

# TODO: Unit test this function this week; use simple domain

# This is the Born approximation of the first order, 
# as it solves the integral equations assuming the incident
# field is an approximation to the total electric field
def solve_first_order(tx):
    # Estimate permittivity using the min value provided
    a_1 = p.eps_min * np.ones(p.N, 1) # Permittivity values
    b_1 = np.zeros(p.M,1)

    m = 0
    G_sum = 0
    for tx in range(1,p.L):
        for rx in range(1,p.L):
            m = m + 1
            if rx == tx:
                continue

            obs_pt = get_antenna_coord(rx) # We will have one tx at a time and 11 rx -> M = 11
            # This loop goes over each element in N = I * J pixels
            for n in range(1,p.N): 

                # This may not be necessary. Depends on how data is arranged
                # Source coordinate
                # Find coordinates of pixel
                src_pt = get_coord(n)

                # Get data for this pixel 
                Ez_r = get_field_data(rx, n)

                #Perform discretized surface integral over x and y range
                K_sum = 0                    
                for x in range(p.dx/2, p.I - p.dx/2, p.dx): # Use centers of each pixel
                    for y in range(p.dy/2, p.J - p.dy/2, p.dy):
                        coord = [x,y]
                        K_sum = K_sum + greens(obs_pt, coord)

                # Add Kji element to K matrix 
                K[m,n] = Ez_r * p.k^2 * K_sum

    # Estimate initial distribution function using zeros vector
    # TODO: check this
    a_1 = K / b

    return a_1

# BIM Algorithm:
# 1. Solve linear inverse problem for first order using Born approximation
# 2. Solve forward scattering problem at object and observation points ( USING SOLVER )
# 3. Substitute field (SOLVED FORWARD) in 2 into integrand in integral eq 
#       and solve inverse problem (WITH_REGULARIZATION) to recover next order distribution ftn
# 4. Repeat 2 and compare fields with measured data. If difference < 5% 
#    of scattered field, terminate successfully. Otherwise repeat until soln converges


# params is a struct containing all required parameters (eg. I, J, Epsilon, etc)
def bim(ei, es, params: Params):

    p = params
    k = p.k

    # Perform first-order born approximation to find initial permittivity solution
    a = solve_first_order()

    K = np.zeros(p.M, p.N) # Initialize the MxN matrix
    b = np.zeros(p.M,1) # Mx1 vector

    # Correct location?
    # Loop until solution approximates measured data
    while True:    
        # Measurement point that tracks location in array 
        m = 0

        # M: # of transmitters x # of receivers
        # Assume all transmitters are capable of receiving
        for tx in range(1,p.L): # Available Transmitters
            for rx in range(1,p.L): # Available receivers
                # increment measurement point
                m = m + 1

                if rx == tx:
                    b[m] = 0
                    continue # Antenna cannot receive and transmit simultaneously

                # Observation coordinate
                obs_pt = get_antenna_coord(rx) # We will have one tx at a time and 11 rx -> M = 11

                Ez_s = 0


                # This loop goes over each element in N = I * J pixels
                for n in range(1,p.N): 

                    # This may not be necessary. Depends on how data is arranged
                    # Source coordinate
                    # Find coordinates of pixel
                    src_pt = get_coord(n)

                    # __________________________________________
                    #
                    # Step 2: Solve forward scattering problem
                    # __________________________________________
                    # Get data for this pixel 
                    Ez_r = get_field_data(rx, n)


                    # Perform discretized surface integral over x and y range
                    K_sum = 0                    
                    for x in range(p.dx/2, p.I - p.dx/2, p.dx): # Use centers of each pixel
                        for y in range(p.dy/2, p.J - p.dy/2, p.dy):
                            coord = [x,y]
                            K_sum = K_sum + greens(obs_pt, coord)

                    # Add Kji element to K matrix 
                    K[m,n] = Ez_r * k^2 * K_sum


        # ____________________________
        #
        # STEP 3: Solve matrix eq for permittivity coefficients
        # ____________________________
        # Calculate LHS (Total electric field) using previous permittivity solution to solve forward problem
        b = K * a

        # ----------------------
        # Perform regularization
        # ----------------------
        H = np.matrix(np.matlib.identity(N))
        H_t = np.transpose(H)
        K_t = np.transpose(K)   

        # Define arbitrary regularization param
        gamma = 1E-12 # Should be btn 1e-10 and 1e-15

        # Regularization matrix
        R = K_t.dot(K) + gamma * H_t.dot(H)

        # ------------------------------
        # Solve for permittivity profile
        # ------------------------------
        a = R^-1 * K_t * b

        # ____________________________
        #
        # STEP 4: Solve forward problem to recalculate Ez_s
        # ____________________________
        # 
        b = K * a
                
        # ____________________________
        #
        # STEP 5: Compare to measured data
        # ____________________________
        # Calculate error: 

        # TODO: need to compare all data at once, or piece by piece
        err = (Ez_s - get_field_data())/Ez_s

        if err < 0.05:
            break

