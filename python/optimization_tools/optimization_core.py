"""optimization_algorithms.py

Script with optimization algorithms used in this work.

luizfelipe.coelho@smt.ufrj.br
Jul 31, 2023
"""


import numpy as np
from numba import jit


def eliminate_equality(H:np.ndarray, p:np.ndarray, A:np.ndarray,
                       b:np.ndarray, C=None, d=None):
    """Method to eliminate the equality constraints from the convex
    optmization problem. Following Eq. (13.6a) from Antoniou, A. and
    Lu W. "Practical Optimization: Algorithms and Engineering
    Applications", 2nd Ed., Springer, 2021.

    Converts the problem:
        min 1/2 x.T H x + x.T p
        s.t. A x = b
    into
        min 1/2 x.T hat_H x + x.T hat_p
    eliminating, more likely implying, the equality constraints.
    
    Parameters
    ----------
    H : np.ndarray
        (Positive semidefinite) Hessian matrix of the optimization
        problem.
    p : np.ndarray
        Input vector.
    A : np.ndarray
        Linear equality constraint matrix (full row rank).
    b : np.ndarray
        Input vector.
    C : np.ndarray (default=None)
        Linear inequality constraint matrix.
    d : np.ndarray (default=None)
        Values for linear inequality constraints.
    
    Returns
    -------
    hat_H : np.ndarray
        Converted Hessian matrix.
    hat_p : np.ndarray
        Converted input vector.
    hat_C : np.ndarray
        Converted linear inequality constraints.
    hat_D : np.ndarray
        Converted values for linear inequality contraints.
    """

    n_rows, n_cols = A.shape
    Q, hat_R = np.linalg.qr(A.T, mode='complete')
    R = hat_R[:n_rows, 0:n_rows]

    Q_1 = Q[0:n_cols, 0:n_rows]
    Q_2 = Q[0:n_cols, n_rows:n_cols]

    hat_H = np.matmul(Q_2.T, np.matmul(H, Q_2))
    hat_p = np.matmul(Q_2.T, (np.matmul(H, np.matmul(Q_1, np.matmul(
        np.linalg.inv(R).T, b
    ))) + p))

    if C is None and d is None:
        
        return hat_H, hat_p
    
    hat_C = C @ Q_2
    hat_d = d - C@Q_1@np.linalg.inv(R).T @ b

    return hat_H, hat_p, hat_C, hat_d


def recover_variable_eq_trans(x_sol, A, b):
    n_rows, n_cols = A.shape
    Q, hat_R = np.linalg.qr(A.T, mode='complete')
    R = hat_R[:n_rows, 0:n_rows]
    Q_1 = Q[0:n_cols, 0:n_rows]
    Q_2 = Q[0:n_cols, n_rows:n_cols]

    return Q_2@x_sol + Q_1@np.linalg.inv(R).T @ b


def convert_inequality(H:np.ndarray, p:np.ndarray, A:np.ndarray):
    """Method to convert inequality constraints into a set of equality
    constraints by adding the slack vector eta. Based on the solution
    for Example 13.4 from Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.

    In practice, converts H, p, and A so that the constraint set is now:
        
        hat_A hat_x = b, with hat_x > 0

    instead of Ax <= b.

    Parameters
    ----------
    H : np.ndarray
        (Positive semidefinte) Hessian matrix of the optimization
        problem.
    p : np.ndarray
        Linear term of cost function.
    A : np.ndarray
        Linear inequality constraint matrix.
    
    Returns
    -------
    hat_H : np.ndarray
        Converted Heassian matrix.
    hat_p : np.ndarray
        Converted linear term of cost function.
    hat_A : np.ndarray
        Converted linear constrint matrix.
    """

    n_rows, n_cols = A.shape
    hat_A = np.hstack((A, -A, np.eye(n_rows)))
    hat_H = np.vstack((
        np.hstack((H, -H, np.zeros((n_cols, n_rows)))),
        np.hstack((-H, H, np.zeros((n_cols, n_rows)))),
        np.zeros((n_rows, 2*n_cols+n_rows))
    ))

    hat_p = np.vstack((p, -p, np.zeros(p.shape)))

    return hat_H, hat_p, hat_A


def recover_variable_inq_trans(hat_x:np.ndarray, var_len:int):
    """Method to recover solution point from solution using slack
    variables.
    
    Based on Section 11.2.1 from Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.

    x = x^{+} - x^{-}
    with hat_x = np.vstack((x^{+}, x^{-}, \eta))

    Parameters
    ----------
    hat_x : np.ndarray
        Solution point.
    var_len : np.ndarray
        Correct length or solution point.

    Returns
    -------
    x : np.ndarray
        Recovered solution point.
    """

    x = hat_x[0:var_len, :] - hat_x[var_len:2*var_len, :]

    return x


def pd_pf_cqp(H:np.ndarray, p:np.ndarray, A:np.ndarray, b:np.ndarray,
              init_point:tuple, rho=None, epsilon=None):
    """Primal-dual path-following algorithm for convex QP problems.
    
    From Algorithm 13.2 in Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.
    
    Solves the problem:
        min (1/2) x.T H x + x.T p
        s.t. A x = b, x >= 0

    Parameters
    ----------
    H : np.ndarray
        (Positive semidefinite) Hessian matrix of the optimization
        problem.
    p : np.ndarray
        Linear term of cost function.
    A : np.ndarray
        Linear equality constraint matrix.
    b : np.ndarray
        Values for the linear equality.
    init_point : tuple
        Strictly feasible initial point 3 element tuple, with:
            x0 : np.ndarray
                Initial point.
            lambda0 : np.ndarray
                Lagrange multiplier (eq. constraints).
            mu0 : np.ndarray
                Lagranfe multiplier (active ineq. constraints).
    rho : float
        A parameter to calculate of tau's value, with
        rho >= sqrt(len(x)) (Eq. 13.36).
    epsilon : float
        Tolerance for duality gap, 1e-9 <= epsilon <= 1e-6 is good
        practice.

    Returns
    -------
    x_sol : np.ndarray
        Solution vector.
    y_sol : float
        Cost function value at the solution vector.
    n_iter : int
        Number of iterations before stop.
    """

    # Initialization:
    x0, lamb0, mu0 = init_point
    if rho is None:
        rho = np.ceil(np.sqrt(len(x0)))
    if epsilon is None:
        epsilon = 1e-6
    den = len(x0) + rho  # Saves computation
    gap = np.matmul(x0.T, mu0)
    mu = mu0.copy()
    x_sol = x0.copy()
    lamb = lamb0.copy()
    n_iter = 0
    # Computing:
    while gap > epsilon:
        n_iter += 1
        tau = gap/den
        M = np.diagflat(mu)
        X = np.diagflat(x_sol)
        Gamma = np.linalg.inv(M + np.matmul(X, H))
        Gamma_aux = np.matmul(Gamma, np.matmul(X, A.T))  # Saves computation
        Y = np.matmul(np.linalg.inv(np.matmul(A, Gamma_aux)), A)
        y = np.matmul(Gamma, (x_sol*mu - tau))
        delta_lamb = np.matmul(Y, y)
        delta_x = np.matmul(Gamma_aux, delta_lamb) - y
        delta_mu = np.matmul(H, delta_x) - np.matmul(A.T, delta_lamb)
        alpha_k = (1 - 1e-6) * np.min(np.vstack((
            np.min(-x_sol[delta_x < 0] / delta_x[delta_x < 0]),
            np.min(-mu[delta_mu < 0] / delta_mu[delta_mu < 0])
        )))
        x_sol = x_sol + alpha_k*delta_x
        lamb = lamb + alpha_k*delta_lamb
        mu = mu + alpha_k*delta_mu
        gap = np.matmul(x_sol.T, mu)
    y_sol = .5* np.matmul(x_sol.T, np.matmul(H, x_sol)) + np.matmul(x_sol.T, p)

    return x_sol, y_sol, n_iter


# @jit(nopython=True)
def nfi_pd_pf_cqp(H:np.ndarray, p:np.ndarray, A:np.ndarray, b:np.ndarray,
                  init_point=None, rho=None, epsilon=None):
    """Nonfeasible-initialization primal-dual path following algorithm
    for convex QP problems.
    
    From Algorithm 13.3 in Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021
    
    Solves the problem:
        min (1/2) x.T H x + x.T p
        s.t. A x = b, x >= 0

    Parameters
    ----------
    H : np.ndarray
        (Positive semidefinite) Hessian matrix of the optimization
        problem.
    p : np.ndarray
        Linear term of cost function.
    A : np.ndarray
        Linear equality constraint matrix.
    b : np.ndarray
        Values for the linear equality.
    init_point : tuple
        Initial point 3 element tuple, with:
            x0 : np.ndarray
                Initial point > 0.
            lambda0 : np.ndarray
                Lagrange multiplier (eq. constraints).
            mu0 : np.ndarray
                Lagranfe multiplier > 0 (active ineq. constraints).
    rho : float
        A parameter to calculate of tau's value, with
        rho >= sqrt(len(x)) (Eq. 13.36).
    epsilon : float
        Tolerance for duality gap, 1e-9 <= epsilon <= 1e-6 is good
        practice.

    Returns
    -------
    x_sol : np.ndarray
        Solution vector.
    y_sol : float
        Cost function value at the solution vector.
    n_iter : int
        Number of iterations before stop.
    """

    # Initialization:
    n_rows, n_cols = A.shape
    if init_point is None:
        x0 = .5*np.ones((n_cols, 1), dtype=np.float64)
        lamb0 = .5*np.ones((n_rows, 1), dtype=np.float64)
        mu0 = .5*np.ones((n_cols, 1), dtype=np.float64)
    else:
        x0, lamb0, mu0 = init_point
    if not p.any():  # Allocate 0's if input is empty vector.
        p = np.zeros((n_cols, 1), dtype=np.float64)
    if rho is None:
        rho = n_cols + 20*np.sqrt(n_cols)
    if epsilon is None:
        epsilon = 1e-5
    den = n_cols + rho
    x_sol = x0.copy()
    lamb = lamb0.copy()
    mu = mu0.copy()
    gap = x_sol.T@mu
    n_iter = 0
    # Computing:
    while gap > epsilon:
        n_iter += 1
        tau = gap/den
        r_d = H@x_sol + p - A.T@lamb - mu
        r_p = b - A@x_sol
        M = np.diag(mu.flatten())
        X = np.diag(x_sol.flatten())
        Gamma = np.linalg.inv(M + X@H)
        Gamma_aux = Gamma@X@A.T  # Saves computation
        Y_0 = np.linalg.inv(A@Gamma_aux)
        y_d = Gamma@ (x_sol*(mu+r_d) - tau)
        delta_lamb = Y_0 @ (A@y_d + r_p)
        delta_x = Gamma_aux@delta_lamb - y_d
        delta_mu = H@delta_x - A.T@delta_lamb + r_d
        alpha_k = (1 - 1e-6) * np.min(np.array((
            np.min(-x_sol.flatten()[delta_x.flatten() < 0]
                   / delta_x.flatten()[delta_x.flatten() < 0], initial=1e12),
            np.min(-mu.flatten()[delta_mu.flatten() < 0]
                   / delta_mu.flatten()[delta_mu.flatten() < 0], initial=1e12)
        )))
        x_sol = x_sol + alpha_k*delta_x
        lamb = lamb + alpha_k*delta_lamb
        mu = mu + alpha_k*delta_mu
        gap = x_sol.T @ mu
        print(gap)
        if n_iter == 40:
            return x_sol, 10, n_iter
        
    y_sol = (.5*x_sol.T@H@x_sol) + (x_sol.T@p)

    return x_sol, y_sol, n_iter


def nfi_ip_mlc_cqp(H:np.ndarray, p:np.ndarray, A:np.ndarray, b:np.ndarray,
                   init_point=None, rho=None, epsilon=None):
    """Nonfeasible-initialization interior-point algorithm for mixed LC
    problems.

    From Algorithm 13.4 in Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021
    
    Parameters
    ----------
    H : np.ndarray
        (Positive semidefinite) Hessian matrix of the optimization
        problem.
    p : np.ndarray
        Linear term of cost function.
    A : np.ndarray
        Linear equality constraint matrix.
    b : np.ndarray
        Values for the linear equality.
    init_point : tuple
        Initial point 3 element tuple, with:
            x0 : np.ndarray
                Initial point > 0.
            lambda0 : np.ndarray
                Lagrange multiplier (eq. constraints).
            mu0 : np.ndarray
                Lagranfe multiplier > 0 (active ineq. constraints).
    rho : float
        A parameter to calculate of tau's value, with
        rho >= sqrt(len(x)) (Eq. 13.36).
    epsilon : float
        Tolerance for duality gap, 1e-9 <= epsilon <= 1e-6 is good
        practice.

    Returns
    -------
    x_sol : np.ndarray
        Solution vector.
    y_sol : float
        Cost function value at the solution vector.
    n_iter : int
        Number of iterations before stop.
    """

    # Initialization:
    n_rows, n_cols = A.shape
    if init_point is None:
        x0 = .5*np.ones((n_cols, 1), dtype=np.float64)
        lamb0 = .5*np.ones((n_rows, 1), dtype=np.float64)
        mu0 = .5*np.ones((n_cols, 1), dtype=np.float64)
    else:
        x0, lamb0, mu0 = init_point
    if not p.any():  # Allocate 0's if input is empty vector.
        p = np.zeros((n_cols, 1), dtype=np.float64)
    if rho is None:
        rho = n_cols + 20*np.sqrt(n_cols)
    if epsilon is None:
        epsilon = 1e-5
    den = n_cols + rho
    x_sol = x0.copy()
    lamb = lamb0.copy()
    mu = mu0.copy()
    gap = x_sol.T@mu
    n_iter = 0
    I = np.eye(n_cols)
    Z = np.zeros((n_rows, n_cols))
    Z2 = np.zeros((n_rows, n_rows))
    e = np.ones((n_cols, 1))
    while gap > epsilon:
        n_iter += 1
        tau = gap/den
        r1 = H @ x_sol + A.T @ lamb - mu + p
        r2 = A@x_sol + Z2@lamb - b
        C = np.vstack((np.hstack((-H, -A.T, I)),
                       np.hstack((-A, Z2, Z)),
                       np.hstack((np.diagflat(mu), Z.T, np.diagflat(x_sol)))))
        delta_w = np.linalg.inv(C)@np.vstack((r1, r2, tau*e-x_sol*mu))
        delta_x = delta_w[0:n_cols]
        delta_lamb = delta_w[n_cols:n_cols+n_rows]
        delta_mu = delta_w[n_cols+n_rows:2*n_cols+n_rows]
        alpha_k = (1 - 1e-6) * np.min(np.array((
            np.min(-x_sol.flatten()[delta_x.flatten() < 0]
                   / delta_x.flatten()[delta_x.flatten() < 0], initial=1e12),
            np.min(-mu.flatten()[delta_mu.flatten() < 0]
                   / delta_mu.flatten()[delta_mu.flatten() < 0], initial=1e12)
        )))
        x_sol = x_sol + alpha_k*delta_x
        lamb = lamb + alpha_k*delta_lamb
        mu = mu + alpha_k*delta_mu
        gap = x_sol.T @ mu
        print(gap)
        if n_iter == 40:
            return x_sol, 10, n_iter
    y_sol = (.5*x_sol.T@H@x_sol) + (x_sol.T@p)

    return x_sol, y_sol, n_iter


# EoF
