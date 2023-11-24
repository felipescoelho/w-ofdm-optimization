"""qp_sdp_socp.py

Script with algorithms for convex quadratic problems.

luizfelipe.coelho@smt.ufrj.br
Out 25, 2023
"""


import numpy as np
from .utils import __svec, __mat, __symmetric_kron


def pd_pf_cqp(H:np.ndarray, p:np.ndarray, A:np.ndarray, b:np.ndarray,
              init_point:tuple, rho=None, epsilon=1e-6):
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
        rho = 20*np.sqrt(len(x0))
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

    return x_sol, n_iter


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
        x0 = np.ones((n_cols, 1), dtype=np.float64)
        lamb0 = np.ones((n_rows, 1), dtype=np.float64)
        mu0 = np.ones((n_cols, 1), dtype=np.float64)
    else:
        x0, lamb0, mu0 = init_point
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

    return x_sol, n_iter


def nfi_ip_mlc(K_11:np.ndarray, K_12:np.ndarray, K_21:np.ndarray,
               K_22:np.ndarray, q_1:np.ndarray, q_2:np.ndarray, init_point=None,
               rho=None, epsilon=1e-6):
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
    n_rows, n_cols = K_21.shape
    if init_point is None:
        x0 = np.ones((n_cols, 1), dtype=np.float64)
        lamb0 = np.ones((n_rows, 1), dtype=np.float64)
        mu0 = np.ones((n_cols, 1), dtype=np.float64)
    else:
        x0, lamb0, mu0 = init_point
    if rho is None:
        rho = n_cols + 15*np.sqrt(n_cols)
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
    e = np.ones((n_cols, 1))
    while gap > epsilon:
        n_iter += 1
        tau = gap/den
        r1 = K_11 @ x_sol + K_12 @ lamb - mu + q_1
        r2 = -K_21@x_sol - K_22@lamb - q_2
        C = np.vstack((np.hstack((-K_11, -K_12, I)),
                       np.hstack((K_21, K_22, Z)),
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
        # print(gap)

    return x_sol, n_iter


def pd_pf_sdp(A_i, b, C, init_set, gamma, sigma=None, epsilon=1e-6):
    """Primal-dual path-following algorithm for SDP problem.

    From Algorithm 13.5 in Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021

    Solves:
    min C X
    s.t. A X = b

    Parameters
    ----------
    A : list
    b : np.ndarray
    C : np.ndarray
    init_set : tuple
    gamma : float
    sigma : float (default=None)
    epsilon : float (default=1e-6)

    Returns
    -------
    X : np.ndarray
    y : np.ndarray
    S : np.ndarray
    n_iter : int
    """
    
    # Initialization:
    n = C.shape[0]
    if sigma is None:
        sigma = n/(15*np.sqrt(n) + n)
    p = len(b)
    X0, y0, S0 = init_set
    X = X0.copy()
    y = y0.copy()
    S = S0.copy()
    A = np.zeros((p, int(n*(n+1)/2)), dtype=np.float64)
    for idx, Ap in enumerate(A_i):
        A[idx, :] = __svec(Ap).T  # 13.81a
    I = np.eye(n)
    gap = np.trace(X@S)/n
    n_iter = 0
    while gap > epsilon:
        n_iter += 1
        tau = sigma*gap  # 13.88
        E = __symmetric_kron(S, I)  # 13.84c
        F = __symmetric_kron(X, I)  # 13.84d
        x = __svec(X)  # 13.85b
        rd = __svec(C - S + __mat(A.T@y))  # 13.85f
        rp = b - A@x  # 13.85c
        rc = __svec(tau*I - .5*(X@S + S@X))  # 13.84e
        E_inv = np.linalg.inv(E)
        M = A@E_inv@F@A.T
        Dy = -np.linalg.inv(M)@(rp + A@E_inv@(F@rd - rc))  # 13.87c
        Dx = -E_inv@(F@(rd + A.T@Dy) - rc)  # 13.87a
        Ds = rd + A.T@Dy
        DX = __mat(Dx)
        DS = __mat(Ds)
        Xk_hat_inv = np.linalg.inv(np.linalg.cholesky(X))
        Sk_hat_inv = np.linalg.inv(np.linalg.cholesky(S))
        min_eig_x = np.min(np.linalg.eigvals(Xk_hat_inv@DX@Xk_hat_inv.T))
        min_eig_s = np.min(np.linalg.eigvals(Sk_hat_inv@DS@Sk_hat_inv.T))
        alpha = 1 if min_eig_x >= 0 else np.min((1, -gamma/min_eig_x))  # 13.91
        beta = 1 if min_eig_s >= 0 else np.min((1, -gamma/min_eig_s))  # 13.92
        X += alpha*DX
        y += beta*Dy
        S += beta*DS
        gap = np.trace(X@S)/n

    return X, y, S, n_iter

def pc_sdp(A_list, b, C, init_set, gamma, epsilon=1e-6):
    """Predictor-corrector algorithm for SDP problems.
    
    From Algorithm 13.6 in Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021

    Parameters
    ----------
    A : list
    b : np.ndarray
    C : np.ndarray
    init_set : tuple
    gamma : float (.9 <= gamma <= .99 is good practice)
    sigma : float (default=None)
    epsilon : float (default=1e-6)

    Returns
    -------
    X : np.ndarray
    y : np.ndarray
    S : np.ndarray
    n_iter : int
    """

    n = C.shape[0]
    p = len(b)
    X0, y0, S0 = init_set
    X = X0.copy()
    y = y0.copy()
    S = S0.copy()
    A = np.zeros((p, int(n*(n+1)/2)), dtype=np.float64)
    for idx, Ap in enumerate(A_list):
        A[idx, :] = __svec(Ap).T  # 13.81a
    I = np.eye(n)
    gap = np.trace(X@S)/n
    print(gap)
    n_iter = 0
    while gap > epsilon:
        # Predictor
        print(gap)
        n_iter += 1
        E = __symmetric_kron(S, I)  # 13.84c
        F = __symmetric_kron(X, I)  # 13.84d
        x = __svec(X)  # 13.85b
        rd = __svec(C - S + __mat(A.T@y))  # 13.85f
        rp = b - A@x  # 13.85c
        rc = __svec(- .5*(X@S + S@X))  # 13.94
        E_inv = np.linalg.inv(E)
        M = A@E_inv@F@A.T
        Dy = -np.linalg.inv(M)@(rp + A@E_inv@(F@rd - rc))  # 13.87c
        Dx = -E_inv@(F@(rd + A.T@Dy) - rc)  # 13.87a
        Ds = rd + A.T@Dy
        DX = __mat(Dx)
        DS = __mat(Ds)
        Xk_hat_inv = np.linalg.inv(np.linalg.cholesky(X))
        Sk_hat_inv = np.linalg.inv(np.linalg.cholesky(S))
        min_eig_x = np.min(np.linalg.eigvals(Xk_hat_inv@DX@Xk_hat_inv.T))
        min_eig_s = np.min(np.linalg.eigvals(Sk_hat_inv@DS@Sk_hat_inv.T))
        alpha_p = 1 if min_eig_x >= 0 else np.min((1, -gamma/min_eig_x))  # 13.95a
        beta_p = 1 if min_eig_s >= 0 else np.min((1, -gamma/min_eig_s))  # 13.95b
        sigma = (np.trace((X+alpha_p*DX)@(S + beta_p*DS))/np.trace(X@S))**3
        tau = sigma*gap  # 13.88
        # Corrector
        rc = __svec(tau*I - .5*(X@S + S@X + DX@DS + DS@DX))
        Dy = -np.linalg.inv(M)@(rp + A@E_inv@(F@rd - rc))  # 13.87c
        Dx = -E_inv@(F@(rd + A.T@Dy) - rc)  # 13.87a
        Ds = rd + A.T@Dy
        DX = __mat(Dx)
        DS = __mat(Ds)
        Xk_hat_inv = np.linalg.inv(np.linalg.cholesky(X))
        Sk_hat_inv = np.linalg.inv(np.linalg.cholesky(S))
        min_eig_x = np.min(np.linalg.eigvals(Xk_hat_inv@DX@Xk_hat_inv.T))
        min_eig_s = np.min(np.linalg.eigvals(Sk_hat_inv@DS@Sk_hat_inv.T))
        alpha_c = 1 if min_eig_x >= 0 else np.min((1, -gamma/min_eig_x))  # 13.100a
        beta_c = 1 if min_eig_s >= 0 else np.min((1, -gamma/min_eig_s))  # 13.100b
        X += alpha_c*DX  # 13.99a
        y += beta_c*Dy  # 13.99b
        S += beta_c*DS  # 13.99c
        gap = np.trace(X@S)/n
        
        
    return X, y, S, n_iter


# EoF
