import time
import numpy as np

'''
This script contains the code for Exercise 2.
'''

## 2a)

# Define the equilibrium function
def equilibrium_1(T):
    term1 = psi * Tc * k 
    term2 = (0.684 - 0.0416 * np.log(T / (1e4 * Z * Z))) * T * k 
    return term1 - term2

# Measure the time taken and find the root of the equilibrium function using the bisection method
def bisection_method(f, a, b, tol=1e-10, max_iter=100):
    '''
    Find the root of the function f(x) = 0 using the bisection method
    Takes the function f, the endpoints of the initial interval [a, b],
    and optional arguments tol and max_iter
    '''
    if f(a) * f(b) >= 0:
        print("Bisection method fails.")
        return None, None
    
    # Initialize values
    num_steps = 0
    start_time = time.time()
    while (b - a) / 2 > tol and num_steps < max_iter:
        c = (a + b) / 2
        if f(c) == 0:
            end_time = time.time()
            return c, num_steps + 1, end_time - start_time
        
        # Update the interval
        if f(c) * f(a) < 0:
            b = c
        else:
            a = c
        num_steps += 1

    end_time = time.time()
    return (a + b) / 2, num_steps, end_time - start_time

# Define the constants
k = 1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s
Z = 0.015 # metallicity 
Tc = 1e4 * Z**2 # stellar temperature in K 
psi = 0.929 

# Measure the time taken and find the root of the equilibrium function using the bisection method
a = 1
b = 1e7
T_eq, num_steps, time_taken = bisection_method(equilibrium_1, a, b)

# Print the result
if T_eq is not None:
    print(f"The equilibrium temperature is {T_eq:.2f} K.")
    print(f"The bisection method found the root in {num_steps} steps.")
    print(f"The time taken was {time_taken:.10f} seconds.")
else:
    print("The bisection method did not converge.")
    
## 2b)
        
# Define the equilibrium function
def equilibrium_2(T, n_e):
    T4 = T / 1e4
    term1 = 0.54 * T4**0.37 * aB * n_e * n_e * k * T
    term2 = A * n_e * epsilon_CR
    term3 = 8.9e-26 * n_e * T4
    return term1 - term2 - term3

# Bisection method
def bisection_method_ne(f, a, b, n_e, tol=1e-10, max_iter=100):
    # Check if the function values at the interval endpoints have the same sign
    if f(a, n_e) * f(b, n_e) >= 0:
        print("Bisection method fails.") 
        return None, None, None  

    num_steps = 0 
    start_time = time.time()  
    # Iterate until convergence or maximum iterations
    while (b - a) / 2 > tol and num_steps < max_iter: 
        # Midpoint of the interval
        c = (a + b) / 2  
         # Check if the midpoint is already the root
        if f(c, n_e) == 0: 
            end_time = time.time() 
            return c, num_steps + 1, end_time - start_time 

        # Oheck if the root lies between a and c
        if f(c, n_e) * f(a, n_e) < 0:  
             # Update the upper bound 
            b = c 
        else:
            # Update the lower bound
            a = c 
        num_steps += 1  

    end_time = time.time()  
    return (a + b) / 2, num_steps, end_time - start_time 

# Secant method
def secant_method(f, x0, x1, n_e, tol=1e-10, max_iter=100):
    num_steps = 0
    start_time = time.time()
    while num_steps < max_iter:
        x_new = x1 - f(x1, n_e) * (x1 - x0) / (f(x1, n_e) - f(x0, n_e))
        if abs(x_new - x1) < tol:
            end_time = time.time()
            return x_new, num_steps + 1, end_time - start_time
        
        x0, x1 = x1, x_new
        num_steps += 1

    end_time = time.time()
    return x_new, num_steps, end_time - start_time

# Densities to consider
densities = [1e-4, 1, 1e4]

with open('2b.txt', 'w') as f:
    for n_e in densities:
        # Initial interval for bisection method
        a_bisection = 1
        b_bisection = 1e15
        T_eq_bisection, num_steps_bisection, time_taken_bisection = bisection_method_ne(equilibrium_2, a_bisection, b_bisection, n_e)

            # Write the result to the file
            if T_eq_bisection is not None:
                f.write(f"For n_e = {n_e} cm^-3, the equilibrium temperature is {T_eq_bisection:.2f} K.\n")
                f.write(f"The bisection method found the root in {num_steps_bisection} steps.\n")
                f.write(f"The time taken was {time_taken_bisection:.6f} seconds.\n\n")
            else:
                f.write(f"For n_e = {n_e} cm^-3, the bisection method failed to converge.\n\n")
    
