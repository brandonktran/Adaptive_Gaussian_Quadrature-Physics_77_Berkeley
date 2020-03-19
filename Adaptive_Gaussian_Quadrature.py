import numpy as np
from math import sqrt
from numpy import *

#regular three-point gaussian quadrature with no regard for tolerance
def gaussian(f,a,b):
    #f: user-defined function that needs to be integrated
    #a,b: the intervals over which to integrate, with a<b
    #this function is just the normal gaussian quadrature calculation with known weights
    #also uses an interval transformation so instead of just the standard interval of -1 to 1, you can use any interval
    u=(b-a)/2.*(5./9*f((b-a)/2.*-1.*sqrt(3./5)+(b+a)/2.)+8./9*f((b+a)/2.)+5./9*f((b-a)/2.*sqrt(3./5)+(b+a)/2.))
    return u

#this is the adaptive part of the gaussian quadrature and allows a tolerance to be met
def sub_division(f,a,b,tol,entire):
    #tol: user-defined tolerance
    #total: the integral on the whole interval a,b using above gaussian function
    #this function splits an integral into the right and left half and compares it to the integral on the whole function
    #if tolerance is not satisfied, the function splits the left and right into smaller intervals and repeats the function
    #when integral interval values satisfy the tolerance, they are added to the sum
    a_z=a+(b-a)/2.  #sub-divides intervals
    b_k=a+(b-a)/2.
    entire=gaussian(f,a,b)
    right=gaussian(f,a_z,b)
    left=gaussian(f,a,b_k)
    if abs(entire-(left+right))<tol*max(abs(entire), (abs(left)+abs(right))):
        return entire
    return sub_division(f,a_z,b,tol,right)+sub_division(f,a,b_k,tol,left)

def adaptive_gaussian_quadrature(f,a,b,tol):
    #returns the approximate integral of f from a to b with an upper bound on the error given by the tolerance
    return sub_division(f,a,b,tol,gaussian(f,a,b))

    
f=lambda x: x**(-x)   #define a function here 
print'Integral value is', adaptive_gaussian_quadrature(f,0,1,10**(-15))  #input intervals and tolerance here 