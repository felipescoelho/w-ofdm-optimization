"""quadratic_programming.py

Script with wrapper for quadratic programming optimization functions.

luizfelipe.coelho@smt.ufrj.br
Oct 6, 2023
"""


import numpy as np
from optimization_core import (eliminate_equality, convert_inequality,
                               recover_variable, nfi_pd_pf_cqp)


def quadratic_solver(H:np.ndarray, p:np.ndarray, A:np.ndarray, b:np.ndarray,
                     C:np.ndarray, d:np.ndarray, rho=None, epsilon=1e-6):
    """Method to solve quadratic programming problems, such as:
    
    min (1/2) x.T H x + x.T p
    s.t. A x = b
         C x <= d
    
    By using variable transformations and a Primal-Dual path following
    algorithm for convex QP problems.
    
    Parameters
    ----------
    H : np.ndarray
        (Positive semi-definite) Hessian matrix of the optimization
        problem.
    p : np.ndarray
        Linear term of the cost function.
    A : np.ndarray
        Linear equality constraints matrix.
    b : np.ndarray
        Values for the linear equality constrains.
    C : np.ndarray
        Linear inequality constraints matrix.
    d : np.ndarray
        Values for the linear inequality constraints.
    rho : float
        A parameter to calculate tau's value, where
        rho >= sqrt(len(x))
    epsilon : float
        Tolerance for duality gap, 1e-9 <= epsilon <= 1e-6 is good
        practice.
        
    Returns
    -------
    x_sol : np.ndarray
        Solution vector.
    n_iter : int
        Number of iterations.
    """

    H_hat0, p_hat0 = eliminate_equality(H, p, A, b)
    H_hat, p_hat, C_hat = convert_inequality(H_hat0, p_hat0, C)
    x_sol_hat, _, n_iter = nfi_pd_pf_cqp(H_hat, p_hat, C_hat, d, None, rho,
                                         epsilon)
    x_sol = recover_variable(x_sol_hat, H.shape[0])

    return x_sol, n_iter


# EoF