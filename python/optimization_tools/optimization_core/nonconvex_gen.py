"""nonconvex_gen.py

Script with algorithms for general nonconvex problems.

luizfelipe.coelho@smt.ufrj.br
Out 25, 2023
"""


import numpy as np
from .utils import (var_trans_equality, recover_var_eq_trans,
                    slack_var_inequality, recover_var_inequality)


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


def __qp_solver(H, p, A, b, C, d, epsilon):
    """A solver for the QP."""
    # Converting problem to correct format:
    H_aux, p_aux, C_aux, d_hat = var_trans_equality(H, p, A, b, C, d)
    var_len = H_aux.shape[0]
    H_hat, p_hat, C_hat = slack_var_inequality(H_aux, p_aux, C_aux)
    # Algorithm settings:
    n_rows, n_cols = C_hat.shape
    xk = np.ones((n_cols, 1), dtype=np.float64)
    lambk = np.ones((n_rows, 1), dtype=np.float64)
    muk = np.ones((n_cols, 1), dtype=np.float64)
    rho = n_cols + 20*np.sqrt(n_cols)
    den = n_cols + rho
    gap = xk.T @ muk
    n_iter = 0  # Start iteration count
    while gap > epsilon:
        n_iter += 1
        tau = gap/den
        rd = H_hat @ xk + p_hat - C_hat.T@lambk - muk
        rp = d_hat - C_hat@xk
        M = np.diagflat(muk)
        X = np.diagflat(xk)
        Gamma = np.linalg.inv(M + X@H_hat)
        Gamma_aux = Gamma@X@C_hat.T
        Y0 = np.linalg.inv(C_hat@Gamma_aux)
        yd = Gamma @ (xk*(muk+rd) - tau)
        d_lamb = Y0 @ (C_hat@yd + rp)
        d_x = Gamma_aux@d_lamb - yd
        d_mu = H_hat@d_x - C_hat.T@d_lamb + rd
        alpha_k = (1 - 1e-6) * np.min((
            np.min(-xk.flatten()[d_x.flatten() < 0]
                    / d_x.flatten()[d_x.flatten() < 0], initial=1e12),
            np.min(-muk.flatten()[d_mu.flatten() < 0]
                    / d_mu.flatten()[d_mu.flatten() < 0], initial=1e12)
        ))
        xk += alpha_k*d_x
        lambk += alpha_k*d_lamb
        muk += alpha_k*d_mu
        gap = xk.T @ muk
    
    x_aux = recover_var_inequality(xk, var_len)
    mu_sol = recover_var_inequality(muk, var_len)
    x_sol = recover_var_eq_trans(x_aux, A, b)

    return x_sol, lambk, mu_sol


def __search_alpha(f, xk, dx, tau):
    """Search alpha in our merit function (L1)."""
    L = 100
    alpha_axis = np.linspace(0, 1, L)
    psi = np.zeros((L, 1), dtype=np.float64)
    for idx, alpha in enumerate(alpha_axis):
        x = xk + alpha*dx
        fval, ak, ck = f(x)
        psi[idx] = fval + tau*np.sum(np.abs(ak)) + tau*np.sum(ck[ck > 0])
    alpha_k = alpha_axis[np.argmin(psi)]

    return alpha_k


def sqp_basic(f, g, h, x0:np.ndarray, lamb0:np.ndarray, mu0:np.ndarray,
              tau, epsilon=1e-6):
    """Basic SQP algorithm for nonconvex problems.
    
    From Algorithm 15.5 in Antoniou, A. and Lu W. "Practical
    Optimizaiton: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.

    Parameters
    ----------
    f : function
        A function that returns the value of the function at a given
        point, with the value of each constraint.
    g : function
        A function to returns the gradient at given point for the cost
        function and its constraints.
    x0 : np.ndarray
        A feasible initial point.
    epsilon : float
        Threshold for convergence criterion. Algorithm stops if
        ||x_{k+1} - x_k||_2 < epsilon.
    
    Returns
    -------
    x_sol : np.ndarray
        Solution point.
    n_iter : int
        Number of iterations.
    """

    xk = x0.copy()
    _, ak, ck = f(x0)
    p = len(ak)
    q = len(ck)
    lamb0 = np.ones((p, 1), dtype=np.float64)
    mu0 = np.ones((q, 1), dtype=np.float64)
    lambk = lamb0.copy()
    muk = mu0.copy()
    n_iter = 0
    distance = 1
    while distance >= epsilon:
        n_iter += 1
        _, ak, ck = f(xk)
        gk, Aek, Aik = g(xk)
        H, A, C = h(xk)
        Zk = H + lambk * A + muk * C
        dx, lamb_qp, mu_qp = __qp_solver(Zk, gk, Aek, -ak, Aik, -ck, epsilon)
        dlamb = lamb_qp - lambk
        dmu = mu_qp - muk
        alpha_k = __search_alpha(f, xk, dx, tau)
        xk += alpha_k*dx
        lambk += alpha_k*dlamb
        muk += alpha_k*dmu
        distance = np.linalg.norm(alpha_k*np.vstack((
            dx, dlamb, dmu
        )))
    
    return xk, n_iter


def sqp_practical(f, g, x0, eta, tau, rho, epsilon):
    """Practical SQP algorithm for nonconvex problems.
    
    From Algorithm 15.5 in Antoniou, A. and Lu W. "Practical
    Optimizaiton: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.

    Parameters
    ----------
    f : function
        A function that returns the value of the function at a given
        point, with the value of each constraint.
    g : function
        A function to returns the gradient at given point for the cost
        function and its constraints.
    x0 : np.ndarray
    eta : np.ndarray
    tau : np.ndarray
    rho : np.ndarray
    epsilon : float

    Returns
    -------
    xk : np.ndarray
    n_iter : np.ndarray
    """

    xk = x0.copy()
    _, ak, ck = f(x0)
    p = len(ak)
    q = len(ck)
    lamb0 = np.ones((p, 1), dtype=np.float64)
    mu0 = np.ones((q, 1), dtype=np.float64)
    lambk = lamb0.copy()
    muk = mu0.copy()
    n_iter = 0
    distance = 1
    while distance >= epsilon:
        n_iter += 1
        _, ak, ck = f(xk)
        gk, Aek, Aik = g(xk)
        dx, lamb_qp, mu_qp = __qp_solver(Zk, gk, Aek, -ak, Aik, -ck, epsilon)
        dlamb = lamb_qp - lambk
        dmu = mu_qp - muk
        alpha_k = __search_alpha(f, xk, dx, tau)
        dk = alpha_k*dx
        xk += dk
        lambk += alpha_k*dlamb
        muk += alpha_k*dmu
        # Modified BFGS
        gk1, Aek1, Aik1 = g(xk)
        gammak = (gk1 - gk) + (Aek1 - Aek).T@lambk + (Aik1 - Aik).T@muk
        aux0 = dk.T@gammak
        aux1 = dk.T@Zk@dk
        theta = 1 if aux0 >= .2*aux1 else .8*aux1/(aux1-aux0)
        etak = theta*gammak + (1-theta)*Zk@dk
        Zk = Zk + etak@etak.T/(dk.T@etak) - Zk@dk@dk.T@Zk/aux1
        # Distance from last point:
        distance = np.linalg.norm(alpha_k*np.vstack((
            dx, dlamb, dmu
        )))
    
    return xk, n_iter


# EoF