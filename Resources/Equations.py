"""
Equations
- simple equations (exp, sin) that are rewritten to be used by the Crocodile.Resources.Mathematics.fit() function
- complicated equations

"""


from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import math
import scipy.special as SS


# SIMPLE EQUATIONS

def quadratic(A, t):
    return A[0] + A[1] * t + A[2] * t**2  

def rb_cos(A, t):
    """
    4 parameters
    function: A[0] + A[1] * numpy.cos(2 * numpy.pi * A[2] * t + A[3])
    A[0]: offset
    A[1]: amplitude
    A[2]: frequency
    A[3]: phase
    """
    return A[0] + A[1] * numpy.cos(2 * numpy.pi * A[2] * t + numpy.pi*A[3])


def linear(A, t):
    return A[0] + A[1] * t  


def rb_exp(A,t):
    return A[0] * numpy.exp(-t / A[1]) 


def double_exp(A,t):
    return A[0] * numpy.exp(-t / A[1]) + A[2] * numpy.exp(-t / A[3])


def single_exp_offset(A,t):
    return A[0] * numpy.exp(-t / A[1]) + A[2]


# DISTRIBUTIONS
def rb_gaussian(A, t):
    """
    A[0]: sigma (sigma^2 = variance)
    A[1]: mu (mean)
    A[2]: offset 
    A[3]: scale, before offset

    """
    y = ( 
        A[3] / (A[0] * numpy.sqrt(2*numpy.pi))
    ) * numpy.exp( 
        -(t - A[1])**2 / (2 * A[0]**2)
    ) + A[2]
    return y



def rb_lorentzian(A, t):
    """
    A[0]: gamma
    A[1]: mean
    A[2]: offset
    A[3]: scale
    """

    return A[3]/(numpy.pi * A[0] * (1 + ((t - A[1])/A[0])**2)) + A[2]


def rb_two_lorentzians(A, t,):

    return rb_lorentzian(A[:4], t) + rb_lorentzian(A[4:], t)





# RESPONSE FUNCTIONS

# to calculate the non-rephasing and rephasing diagrams
def g(t, delta, t_corr):
    """
    20101209/RB: started/continued
    These are the g(t) functions used to make the non-rephasing and rephasing diagrams.

    INPUT:
    t: time array

    OUTPUT:
    An array with the results
    """
    # delta= 1/1000
    # t_corr = 1000
    return delta**2 * t_corr**2 * ( numpy.exp(-t/t_corr) + t/t_corr - 1)

# calculate the non-rephasing diagram
def non_rephasing(t1, t2, t3, w, anh, delta, t_corr):
    """
    20101209/RB: started/continued
    Calculates the non-rephasing diagram

    INPUT:
    t1, t3 (mesh): coherence times in fs
    t1_axis = numpy.arange(start, stop, step)
    t3_axis = numpy.arange(start, stop, step)
    t1, t3 = numpy.meshgrid(t1_axis, t3_axis)

    t2 (integer/float?): population time in fs
    w (integer): frequency (in cm-1)
    anh (integer): anharmonicity (in cm-1)
    delta (float): ?
    t_corr (float): correlation time

    OUTPUT:
    A 2D array of the response in for the coherence times
    """

    return (numpy.exp(-1j * (t3 + t1) * w) - numpy.exp(-1j * (t3 * (w - anh) + t1 * w))) * (numpy.exp(-g(t1, delta, t_corr) - g(t2, delta, t_corr) - g(t3, delta, t_corr) + g(t1+t2, delta, t_corr) + g(t2+t3, delta, t_corr) - g(t1+t2+t3, delta, t_corr)))

#calculate the rephasing diagram
def rephasing(t1, t2, t3, w, anh, delta, t_corr):
    """
    20101209/RB: started/continued
    Calculates the non-rephasing diagram

    see for details the non_rephasing function
    """
    return (numpy.exp(-1j * (t3 - t1) * w) - numpy.exp(-1j * (t3 * (w - anh) - t1 * w))) * (numpy.exp(-g(t1, delta, t_corr) + g(t2, delta, t_corr) - g(t3, delta, t_corr) - g(t1+t2, delta, t_corr) - g(t2+t3, delta, t_corr) + g(t1+t2+t3, delta, t_corr)))


# REFRACTIVE INDEX

def gvd_formula_1(x, s): 
    """
    Formula 1 and 2 are the same, except that some coefficients are squared. 
    """
    if len(s) < 9:
        _s = numpy.zeros(9)
        _s[:len(s)] = s
        s = _s

    y = (
        (6 * s[1] * s[2]**2) / (x**4 * (1 - s[2]**2 / x**2)**2) 
        + (8 * s[1] * s[2]**4) / (x**6 * (1 - s[2]**2 / x**2)**3) 
        + (6 * s[3] * s[4]**2) / (x**4 * (1 - s[4]**2 / x**2)**2) 
        + (8 * s[3] * s[4]**4) / (x**6 * (1 - s[4]**2 / x**2)**3) 
        + (6 * s[5] * s[6]**2) / (x**4 * (1 - s[6]**2 / x**2)**2) 
        + (8 * s[5] * s[6]**4) / (x**6 * (1 - s[6]**2 / x**2)**3) 
        + (6 * s[7] * s[8]**2) / (x**4 * (1 - s[8]**2 / x**2)**2) 
        + (8 * s[7] * s[8]**4) / (x**6 * (1 - s[8]**2 / x**2)**3)
    ) / (
        2 * numpy.sqrt(
            s[1] / (1 - s[2]**2 / x**2) 
            + s[3] / (1 - s[4]**2 / x**2) 
            + s[5] / (1 - s[6]**2 / x**2) 
            + s[7] / (1 - s[8]**2 / x**2) 
            + s[0] 
            + 1
        )
    ) - (
        -(2 * s[1] * s[2]**2)/(x**3 * (1 - s[2]**2 / x**2)**2) 
        - (2 * s[3] * s[4]**2)/(x**3 * (1 - s[4]**2 / x**2)**2) 
        - (2 * s[5] * s[6]**2)/(x**3 * (1 - s[6]**2 / x**2)**2) 
        - (2 * s[7] * s[8]**2)/(x**3 * (1 - s[8]**2 / x**2)**2)
    )**2 / (
        4 * (
            s[1]/(1 - s[2]**2 / x**2) 
            + s[3] / (1 - s[4]**2 / x**2) 
            + s[5] / (1 - s[6]**2 / x**2) 
            + s[7] / (1 - s[8]**2 / x**2) 
            + s[0] 
            + 1
        )**(3/2)
    )
    
    return y


def gvd_formula_3(x, s):

    """
    http://www.wolframalpha.com/input/?i=second+derivative+of+sqrt(a+%2B+b*x%5Ec+%2B+d*x%5Ee+%2B+f*x%5Eg+%2B+h*x%5Ei+%2B+j*x%5Ek+%2B+l*x%5Em+%2B+n*x%5Eo+%2Bp*x%5Eq)+with+respect+to+x
    
    """

    if len(s) < 17:
        _s = numpy.zeros(17, dtype = "float")
        _s[:len(s)] = s
        s = _s

    y = (
            s[1] * (s[2] - 1) * s[2] * x**(s[2] - 2) 
            + s[3] * (s[4] - 1) * s[4] * x**(s[4] - 2) 
            + s[5] * (s[6] - 1) * s[6] * x**(s[6] - 2) 
            + s[7] * (s[8] - 1) * s[8] * x**(s[8] - 2) 
            + s[9] * (s[10] - 1) * s[10] * x**(s[10] - 2) 
            + s[11] * (s[12] - 1) * s[12] * x**(s[12] - 2) 
            + s[13] * (s[14] - 1) * s[14] * x**(s[14] - 2) 
            + s[15] * (s[16] - 1) * s[16] * x**(s[16] - 2)
        ) / (
            2 * numpy.sqrt(
                s[0] 
                + s[1] * x**s[2] 
                + s[3] * x**s[4] 
                + s[5] * x**s[6] 
                + s[7] * x**s[8] 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )
        ) - (
            s[1] * s[2] * x**(s[2] - 1) 
            + s[3] * s[4] * x**(s[4] - 1) 
            + s[5] * s[6] * x**(s[6] - 1) 
            + s[7] * s[8] * x**(s[8] - 1) 
            + s[9] * s[10] * x**(s[10] - 1) 
            + s[11] * s[12] * x**(s[12] - 1) 
            + s[13] * s[14] * x**(s[14] - 1) 
            + s[15] * s[16] * x**(s[16] - 1)
        )**2 / (
            4 * (
                s[0] 
                + s[1] * x**s[2] 
                + s[3] * x**s[4] 
                + s[5] * x**s[6] 
                + s[7] * x**s[8] 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )**(3/2)
        )
    return y


def gvd_formula_4(x, s):
    """
    http://www.wolframalpha.com/input/?i=second+derivative+of+sqrt(a%2B+(b*x%5Ec)%2F(x%5E2-d%5Ee)+%2B+(f*x%5Eg)%2F(x%5E2-h%5Ei)+%2B+j*x%5Ek+%2B+l*x%5Em+%2B+n*x%5Eo+%2B+p*x%5Eq+)+with+respect+to+x
    
    
    """
    if len(s) < 17:
        _s = numpy.zeros(17, dtype = "float")
        _s[:len(s)] = s
        s = _s

    y = (
            (s[1] * (s[2] - 1) * s[2] * x**(s[2] - 2))/(x**2 - s[3]**s[4]) 
            + (8 * s[1] * x**(s[2] + 2))/(x**2 - s[3]**s[4])**3 
            - (2 * s[1] * s[2] * x**s[2])/(x**2 - s[3]**s[4])**2 
            - (2 * s[1] * (s[2] + 1) * x**s[2])/(x**2 - s[3]**s[4])**2 
            + (s[5] * (s[6] - 1) * s[6] * x**(s[6] - 2))/(x**2 - s[7]**s[8]) 
            + (8 * s[5] * x**(s[6] + 2))/(x**2 - s[7]**s[8])**3 
            - (2 * s[5] * s[6] * x**s[6])/(x**2 - s[7]**s[8])**2 
            - (2 * s[5] * (s[6] + 1) * x**s[6])/(x**2 - s[7]**s[8])**2 
            + s[9] * (s[10] - 1) * s[10] * x**(s[10] - 2) 
            + s[11] * (s[12] - 1) * s[12] * x**(s[12] - 2) 
            + s[13] * (s[14] - 1) * s[14] * x**(s[14] - 2) 
            + s[15] * (s[16] - 1) * s[16] * x**(s[16] - 2)
        ) / (
            2 * numpy.sqrt(
                s[0] 
                + (s[1] * x**s[2])/(x**2 - s[3]**s[4]) 
                + (s[5] * x**s[6])/(x**2 - s[7]**s[8]) 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )
        ) - (
            (s[1] * s[2] * x**(s[2] - 1))/(x**2 - s[3]**s[4]) 
            - (2 * s[1] * x**(s[2] + 1))/(x**2 - s[3]**s[4])**2 
            + (s[5] * s[6] * x**(s[6] - 1))/(x**2 - s[7]**s[8]) 
            - (2 * s[5] * x**(s[6] + 1))/(x**2 - s[7]**s[8])**2 
            + s[9] * s[10] * x**(s[10] - 1) 
            + s[11] * s[12] * x**(s[12] - 1) 
            + s[13] * s[14] * x**(s[14] - 1) 
            + s[15] * s[16] * x**(s[16] - 1)
        )**2 / (
            4 * (
                s[0] 
                + (s[1] * x**s[2])/(x**2 - s[3]**s[4]) 
                + (s[5] * x**s[6])/(x**2 - s[7]**s[8]) 
                + s[9] * x**s[10] 
                + s[11] * x**s[12] 
                + s[13] * x**s[14] 
                + s[15] * x**s[16]
            )**(3/2)
        )
    return y


    
    
def gvd_formula_5(x, s):

    """
    Formula 5 
    http://www.wolframalpha.com/input/?i=second+derivative+of+a%2B+b*x%5Ec+%2B+d*x%5Ee+%2B+f*x%5Eg+%2Bh*x%5Ei+%2B+j*x%5Ek+with+respect+to+x
    """
    
    if len(s) < 11:
        _s = numpy.zeros(11, dtype = "float")
        _s[:len(s)] = s
        s = _s  

    y = ( 
        s[1] * (s[2] - 1) * s[2] * x**(s[2] - 2) 
        + s[3] * (s[4] - 1) * s[4] * x**(s[4] - 2) 
        + s[5] * (s[6] - 1) * s[6] * x**(s[6] - 2) 
        + s[7] * (s[8] - 1) * s[8] * x**(s[8] - 2) 
        + s[9] * (s[10] - 1) * s[10] * x**(s[10] - 2)
    )
    
    return y
    
    
def gvd_formula_6(x, s):  
    """
    http://www.wolframalpha.com/input/?i=second+derivative+of+1%2Ba%2B+b%2F(c-x%5E-2)+%2B+d%2F(e-x%5E-2)+%2B+f%2F(g-x%5E-2)+%2B+h%2F(i-x%5E-2)+%2B+j%2F(k-x%5E-2)+with+respect+to+x
    """
    
    if len(s) < 11:
        _s = numpy.zeros(11, dtype = "float")
        _s[:len(s)] = s
        s = _s
        
        
    y = (
        s[1] * (
            8/((x**6) * (s[2] - 1/x**2)**3) 
            + 6/((x**4) * (s[2] - 1/x**2)**2)
        ) + s[3] * (
            8/((x**6) * (s[4] - 1/x**2)**3) 
            + 6/((x**4) * (s[4] - 1/x**2)**2)
        ) + s[5] * (
            8/((x**6) * (s[6] - 1/x**2)**3) 
            + 6/((x**4) * (s[6] - 1/x**2)**2)
        ) + s[7] * (
            8/((x**6) * (s[8] - 1/x**2)**3) 
            + 6/((x**4) * (s[8] - 1/x**2)**2)
        ) + s[9] * (
            8/((x**6) * (s[10] - 1/x**2)**3) 
            + 6/((x**4) * (s[10] - 1/x**2)**2)
        )   
    ) 
    
    return y
 
def gvd_formula_7(x, s):  
    """
    http://www.wolframalpha.com/input/?i=second+derivative+of+a+%2B+b%2F(x%5E2-0.028)%2Bc*(1%2F(x%5E2-0.028))%5E2+%2B+d*x%5E2+%2Be*x%5E4+%2Bf*x%5E6+with+respect+to+x
    """

    if len(s) < 6:
        _s = numpy.zeros(6, dtype = "float")
        _s[:len(s)] = s
        s = _s

    y = (
        s[1] * (
            (8 * x**2)/(x**2 - 0.028)**3 
            - 2/(x**2 - 0.028)**2
        ) 
        + s[2] * (
            (24 * x**2)/(x**2 - 0.028)**4 
            - 4/(x**2 - 0.028)**3
        ) 
        + 2 * s[3] 
        + 12 * s[4] * x**2 
        + 30 * s[5] * x**4   
    )
    
    return y
    
    
    
    
    


