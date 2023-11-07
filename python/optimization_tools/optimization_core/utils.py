"""utils.py

Utilitary functions to the optimization core.

luizfelipe.coelho@smt.ufrj.br
Oct 25, 2023
"""


import numpy as np


def var_trans_equality(H:np.ndarray, p:np.ndarray, A:np.ndarray,
                       b:np.ndarray, C=None, d=None):
    """
    Method to transform the problem with equality and (with or without)
    inequality, into a problem without the equality constraints.
    Following Eq. (13.6a) from Antoniou, A. and Lu W. "Practical
    Optimization: Algorithms and Engineering Applications", 2nd Ed.,
    Springer, 2021.

    Converts the problem:
        min 1/2 x.T @ H @ x + x.T @ p
        s.t. A @ x = b
             C @ x <= d
    into
        min 1/2 x.T @ hat_H @ x + x.T @ hat_p
        s.t. C_hat @ x <= d_hat
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
    hat_d : np.ndarray
        Converted values for linear inequality contraints.
    """

    n_rows, n_cols = A.shape
    Q, hat_R = np.linalg.qr(A.T, mode='complete')
    R = hat_R[:n_rows, 0:n_rows]

    Q_1 = Q[0:n_cols, 0:n_rows]
    Q_2 = Q[0:n_cols, n_rows:n_cols]

    hat_H = Q_2.T@H@Q_2
    hat_p = Q_2.T@(H@Q_1@np.linalg.inv(R).T@b + p)

    if C is None and d is None:
        
        return hat_H, hat_p

    hat_C = C @ Q_2
    hat_d = d - C@Q_1@np.linalg.inv(R).T @ b

    return hat_H, hat_p, hat_C, hat_d


def recover_var_eq_trans(x_sol:np.ndarray, A:np.ndarray, b:np.ndarray):
    """
    Method to recover original solution variable from equality
    constraint transformation.

    Parameters
    ----------
    x_sol : np.ndarray
        Solution vector.
    A : np.ndarray
        Equality constraints matrix
    b : np.ndarray
        Values for the inequality constraints.
    """

    n_rows, n_cols = A.shape
    Q, hat_R = np.linalg.qr(A.T, mode='complete')
    R = hat_R[:n_rows, 0:n_rows]
    Q_1 = Q[0:n_cols, 0:n_rows]
    Q_2 = Q[0:n_cols, n_rows:n_cols]

    return Q_2@x_sol + Q_1@np.linalg.inv(R).T @ b


def slack_var_inequality(H:np.ndarray, p:np.ndarray, A:np.ndarray):
    """
    Method to convert inequality constraints into a set of equality
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
    hat_p = np.vstack((p, -p, np.zeros((n_rows, 1))))

    return hat_H, hat_p, hat_A


def recover_var_inequality(hat_x:np.ndarray, var_len:int):
    """
    Method to recover solution point from solution using slack variable.
    
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


def __svec(K:np.ndarray):
    """
    Symmetric Kronecker vectorization. (simplifies computation)
    
    Converts a symmetric matrix into a vector, eliminating redundant
    elements from the symmetric matrix.

    Parameters
    ----------
    K : np.ndarray
        n x n symmetric matrix.

    Returns
    -------
    k : np.ndarray
        Vectorized symmetric matrix with length n*(n+1)/2.
    """

    n = K.shape[0]
    nn = int(n*(n+1)/2)
    aux = np.sqrt(2)
    k = np.zeros((nn, 1), dtype=np.float64)
    for idx in range(n-1):
        start = int(n*idx)
        end = int((n-1)*idx + n)
        if start < nn:
            k[start] = K[idx, idx]
            k[start+1:end, 0] = aux*K[idx+1:, idx]
    k[-1] = K[-1, -1]

    return k


def __mat(k:np.ndarray):
    """
    Symmetric Kronecker Matrix recover.

    Parameters
    ----------
    k : np.ndarray
        Vectorized symmetric matrix

    K : np.ndarray
        Recovered symmetric matrix.
    """

    nn = len(k)
    n = int((-1+np.sqrt(1+8*nn))/2)
    aux = 1/np.sqrt(2)
    K = np.zeros((n, n), dtype=np.float64)
    for idx in range(n-1):
        start = int(n*idx)
        end = int((n-1)*idx + n)
        if start < nn:
            K[idx, idx] = k[start, 0]
            K[idx+1:, idx] = aux*k[start+1:end, 0]
    K[-1, -1] = k[-1]
    K += K.T - np.diag(np.diag(K))

    return K


def __symmetric_kron(M:np.ndarray, N:np.ndarray):
    """
    Kronecker product of symmetric matrices.
    
    Parameters
    ----------
    M : np.ndarray
    N : np.ndarray

    Returns
    -------
    Y : np.ndarray
    """

    n = M.shape[0]
    nn = int(n*(n+1)/2)
    aux = 1/np.sqrt(2)
    Y = np.zeros((nn, nn), dtype=np.float64)
    idx = -1
    for i in range(n):
        for j in range(i, n):
            K = np.zeros((n, n), dtype=np.float64)
            idx += 1
            K[i, j] = 1 if i == j else aux
            K[j, i] = 1 if i == j else aux
            Y[:, idx] = __svec(.5*(N@K@M.T + M@K@N.T)).flatten()
    
    return Y


# EoF
