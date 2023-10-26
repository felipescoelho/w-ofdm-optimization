"""nonconvex_gen.py

Script with algorithms for general nonconvex problems.

luizfelipe.coelho@smt.ufrj.br
Out 25, 2023
"""


import numpy as np


def scp(f, g, h, x0, rho, epsilon):
    """Sequential convex programming for nonconvex problems.
    
    From Algorithm 15.1 in Antoniou, A. and Lu W. "Practical
    Optimizaiton: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.
    
    Solves the problem:
        min f(x)
        s.t. a_i(x) = 0, for i = 1, 2, ..., p
             c_j(x) <= 0 for j = 1, 2, ..., q
    
    Parameters
    ----------
    f : function
        A function that calculates f(x), a_i(x), and c_j(x).
    g : function
        A function containing the gradients for f(x) and a(x).
    h : function
        A function containing the Hessian.
    x0 : np.ndarray
        Feasible initial point.
    rho : np.ndarray
        Parameters for trust region.
    epsilon : float
        Threshold for convergence criterion. Algorithm stops if
        ||x_{k+1} - x_k||_2 < epsilon.
    """


    xk = x0.copy()
    n_iter = 0
    dist = 1
    while dist < epsilon:
        n_iter += 1
        f_val, a_val, c_val = f(xk)
        gf_val, ga_val = g(xk)

