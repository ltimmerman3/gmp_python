import numpy as np
from numba import jit
import sympy as sym
from math import comb, factorial
from scipy.special import factorial2

class MCSHGaussFxnIntegrals():
    """Gaussian Function Integrals Summary

    Written by Lucas R Timmerman, Medford Group, Februrary 2023

    This class is used to calculate the integrals of the product of gaussian and
    Maxwell Cartesian Spherical Harmonics functions from -inf to inf as implemented in 
    the GMP descriptors. The code is broken up into several functions to make it easier 
    to read and understand. These integrals form the basis of the GMP descriptors. The 
    integrals are calculated using general properties of the gaussian function and the 
    binomial theorem. The goal of this class is not numerical computation, but rather to 
    provide a systematic and reproducible way to validate the most fundamental poriton of 
    the GMP descriptors. This class may also be used with minimal alteration to calculate 
    the integrals of other gaussian based functions."""

    def __init__(self):
        self.L = sym.Symbol('L')
        self.G = sym.Symbol('G')

    def gaussian_function_integral(self,n_omega: int, k: int, cart_component: str):
        """Gaussian Function Integral Summary
        
        Written by Lucas R Timmerman, Medford Group, Februrary 2023
        
        This function calculates the integral of the gaussian function from -inf to inf.
        Note that the factor sqrt(pi / constant) is not included in the integral. It is a
        constant that can be factored out of the remainder of the integral calculation.
        
        Arguments:
            n: Order of the component of MCSH function (x,y,z) for which the integral is being calculated
            k: Order of lamda*r in the gaussian function"""    

        # Integral of odd function * gaussian == 0
        if (n_omega-k)%2 != 0:
            return 0
        else:
            # L = lambda
            L = self.L
            # r = generic cartesian coordinate (x,y,z)
            r = sym.Symbol(cart_component)
            # G = gamma
            G = self.G
            return comb(n_omega,k) * (L*r) ** k * factorial2((n_omega - k - 1),exact=True) / (2 ** ((n_omega-k) / 2) * G ** ((n_omega - k) / 2))

    def sum_of_gaussian_function_integrals(self, n_omega: int, cart_component: str):
        """Sum of Gaussian Function Integrals Summary
        
        Written by Lucas R Timmerman, Medford Group, Februrary 2023
        
        This function calculates the sum of the integrals of the gaussian function
        for a given order of a single cartesian component.
        
        Arguemnts:
            n_omega: Order of the component of MCSH function (x,y,z) for which the integral is being calculated"""

        return sum([self.gaussian_function_integral(n_omega,k,cart_component) for k in range(n_omega+1)])
    
    def get_raw_integrals(self,n_omega: int, cart_component: str):
        """Get Raw Integrals Summary
        
        Written by Lucas R Timmerman, Medford Group, Februrary 2023
        
        This function calls all others to return the integrals of the MCSH * gaussian function
        of a given order for a given cartesian axis
        
        Arguments:
            n_omega: Order of the component of MCSH function (x,y,z) for which the integral is being calculated
            cart_component: Cartesian component of the MCSH function (x,y,z) for which the integral is being calculated

        Returns:
            Expression for the integrals of a single cartesian component of the MCSH * gaussian function
            as a function of:
                L: lambda
                r: generic cartesian coordinate (x,y,z)
                G: gamma"""
        
        return self.sum_of_gaussian_function_integrals(n_omega, cart_component)
    
    def get_n_cartesian_integrals(self,n_omega: int, cart_component: str):
        """Get n Cartesian Integrals Summary
        
        Written by Lucas R Timmerman, Medford Group, Februrary 2023
        
        This function calls all others to return the integrals of the MCSH * gaussian function
        of a given order for all cartesian components
        
        Arguments:
            n_omega: Order of the component of MCSH function (x,y,z) for which the integral is being calculated

        Returns:
            Expression for the integrals of all cartesian components of the MCSH * gaussian function
            as a function of:
                L: lambda
                r: generic cartesian coordinate (x,y,z)
                G: gamma"""
        
        return self.get_raw_integrals(n_omega,cart_component)
    
    def get_all_integrals(self,n_x: int, n_y: int, n_z: int):
        """Get Integrals Summary
        
        Written by Lucas R Timmerman, Medford Group, Februrary 2023
        
        This function calls gives the complete form of the integrals of the MCSH * gaussian function
        for a given group of MCSH order n = n_x + n_y + n_z
        
        Arguments:
            n_x: Order of the x component of MCSH function for which the integral is being calculated
            n_y: Order of the y component of MCSH function for which the integral is being calculated
            n_z: Order of the z component of MCSH function for which the integral is being calculated
            
        Returns:
            Expression for the integrals of all cartesian components of the MCSH * gaussian function for a given group of MCSH order n = n_x + n_y + n_z"""
        return [self.get_n_cartesian_integrals(n_omega, cart_component) for n_omega,cart_component in zip([n_x,n_y,n_z],['x','y','z'])]
    
    def maxwell_cartesian_spherical_harmonics(self,permutation_group, order): 
        """Maxwell Cartesian Spherical Harmonics Summary
        Calculates the spherical harmonics for a given permutation group and overall order.
        Goal is to get polynomial coefficients for the spherical harmonics.
        
        Arguments:
            permutation_group: List of integers that represent the permutation group
            order: Order of the spherical harmonics (not strictly necessary, but useful for readability)
            
        Returns:
            _function: Full expression for MSCH
            _coefficients: Polynomial coefficients for MSCH
            poly_degrees: Polynomial degrees for each cartesian component of MSCH for each term"""
        
        # Variables corresponding to cartesian coordinates
        x = sym.Symbol('x')
        y = sym.Symbol('y')
        z = sym.Symbol('z')

        # radius of sphere on which MCSH is defined
        r = sym.Symbol('r')

        # Overall MCSH order
        n = order

        # Assign degree of each cartesian component based on permutation group
        n_1 = permutation_group[0]
        n_2 = permutation_group[1]
        n_3 = permutation_group[2]

        # Collect all unique terms separately
        terms = []

        # Compuational loop
        for m_1 in range(int(n_1/2+1)):
            for m_2 in range(int(n_2/2+1)):
                for m_3 in range(int(n_3/2+1)):
                    m = m_1 + m_2 + m_3
                    terms += [(-1)**m * factorial2(2*n - 2*m - 1) \
                            * factorial(n_1)/(2**m_1*factorial(m_1)*factorial(n_1-2*m_1)) \
                            * factorial(n_2)/(2**m_2*factorial(m_2)*factorial(n_2-2*m_2)) \
                            * factorial(n_3)/(2**m_3*factorial(m_3)*factorial(n_3-2*m_3)) \
                            * r**(2*m) * x**(n_1 - 2*m_1) * y**(n_2 - 2*m_2) * z**(n_3 - 2*m_3)]
                    
        # Collect all terms into a single expression
        _function = sym.expand(sum(terms)/r**n)

        # Get polynomial coefficients
        coefficients = sym.poly(_function).coeffs()

        # Get polynomial degrees by Cartesian component
        poly_degrees = [{} for i in range(len(coefficients))]
        for i,term in enumerate(terms):
            _poly = sym.poly(term)
            for cart,var in zip([x,y,z],['x','y','z']):
                try:
                    poly_degrees[i][var] = sym.degree(_poly,gen=cart) 
                except:
                    poly_degrees[i][var] = 0
                    
        return _function, coefficients, poly_degrees

if __name__ == "__main__":

    test = MCSHGaussFxnIntegrals()
    
    print(test.maxwell_cartesian_spherical_harmonics([2,2,0], 4))

# Write a python script that calculates the coefficients for the maxwell cartesian spherical harmonics given a permuattion group and an order