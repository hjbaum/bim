import numpy as np
import scipy as sci

global K, b

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
def greens(r_obs, r_src, k):
    # Scalar distance value
    R = np.sqrt((r_obs(1)-r_src(1))^2 + (r_obs(2)-r_src(2))^2 + (r_obs(3)-r_src(3))^2)
    gs = np.exp(-1j * k * R) / (4* np.pi * R)
    return gs

# Discretized version of the dyadic Green function
def greens_disc(r_m, r_n, dx, dy, k):
    an = np.sqrt(dx*dy/np.pi)
    J1 = sci.jv(1,k*an) # first-order bessel function
    
    R = np.sqrt((r_m[0]-r_n[0])^2 + (r_m[1]-r_n[1])^2) # distance vector

    H0_2 = sci.hankel2(0, k * R) # 0th-order, 2nd degree hankel function
    
    green_d = -1j/2 * np.pi * k * an * J1 * H0_2

    return green_d

# This is the Born approximation of the first order, 
# as it solves the integral equations assuming the incident
# field is an approximation to the total electric field
def estimate_initial_contrast(I, J, L, M, dx, dy, k):
    Xn = np.array([])
    n = 0

    for i in range(0,I):
        for j in range(0,J):
            gmn = greens_disc(dx, dy, k)
            Ei = solve_forward()
            Es = solve_inverse(gmn, Xn, Ei)

            # May need to reconfigure this since Es depends on Xn
            Xn[n] = sum(sum(gmn * Ei * Es, 1, L), 1, M)
            n = n + 1

    return Xn


def solve_forward():
    # Solve the forward problem
    Ei = 1
    return Ei

def solve_inverse(gmn, Xn, En, N):
    # Solves the inverse scattering problem
    Es = sum(gmn*Xn*En, 1, N)


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
    I = np.array(0,1,params.I) # use inrange
    J = np.array(0,1,params.J)
    L = np.array(0,1,params.L)
    N = np.array(0,1,params.N)
    M = np.array(0,1,params.M)
    k = params.k

    # Determine contrast function
    estimate_initial_contrast(I, J, L, M, dx, dy, k)

    
    for num_iter in range(1,params.misf):
        for l in range(1,L):
            for m in range(1,M):
                for n in range(1,N):
                    # Calculate forward solution
                    En_tot = 1

                    # Calculate inverse solution
                    # _r denotes relative value at point r
                    Xn = eps_r / eps_rb - 1 - 1j * (sig_r - sigma_b) / (omega * eps_b)
                    Em_s = sum(greens_disc(dx, dy, k) * Xn * En_tot, 1, N)


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