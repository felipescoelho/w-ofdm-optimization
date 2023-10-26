"""seq_quadratic_programmin.py

Script with wrapper for sequential quadratic optimization function.

luizfelipe.coelho@smt.ufrj.br
Oct 23, 2023
"""

import numpy as np
from .optimization_core import p_sqp


def sqp_solver(H:np.ndarray, p:np.ndarray, A:np.ndarray, a: np.ndarray,
               b:np.ndarray, C:np.ndarray, d:np.ndarray, eta:float, tau:float,
               rho=None, epsilon=1e-6):
    """Method to solve the nonconvex problem.

    min (1/2) x.T H x + x.T p
    s.t. x.T A x + x.T a = b
         Cx <= d
    
    Parameters
    ----------
    H : np.ndarray
        Hessian matrix.
    p : np.ndarray
        Linear term of the cos function.
    A : np.ndarray
        Quadratic term of the equality constraints.
    a : np.ndarray
        Linear term of the equality constriants.
    b : np.ndarray
        Values for the equality constriants.
    C : np.ndarray
        Linear inequality constraints.
    d : np.ndarray
        Values for the inequality constraints.
    eta"""


    def f(x):
        """
        Method to evaluate the function's value and its constraints at a
        given point x.
        """
        a_k = x.T@A@x + x.T@a - b
        c_k = C@x - d
        f_k = .5*x.T@H@x + x.T@p
        
        return f_k, a_k, c_k
    
    def g(x):
        """
        Method to evaluate the gradient of the function and its
        constraints for a given point x.
        """
        A_ek = a + A @ x
        A_ik = C
        g_k = p + H @ x

        return g_k, A_ek, A_ik
    
    x_sol, _, n_iter = p_sqp(f, g, eta, tau, rho, (H.shape[0], ))