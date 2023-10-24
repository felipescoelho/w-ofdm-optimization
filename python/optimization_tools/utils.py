"""utils.py

Utilitary functions for the optimization tools.

luizfelipe.coelho@smt.ufrj.br
Oct 7, 2023
"""


import numpy as np


def reduce_variable_tx(dft_len:int, cp_len:int, cs_len:int, tail_len:int):
    """
    Method to remove repetitive samples from our optimization variable
    for the transmitter window optimization.
    
    Parameters
    ----------
    dft_len : int
        Length of the DFT.
    cp_len : int
        Length of the cyclic prefix.
    cs_len : int
        Lenght of the cyclic suffix.
    tail_len : int
        Length of the transmitter window tail.
    
    Returns
    -------
    reduce_var_mat : np.ndarray
        Matrix responsible for reducing the optimization variable.
    """

    return np.vstack((
        np.hstack((np.zeros((tail_len, 1)), np.fliplr(np.eye(tail_len)))),
        np.hstack((np.ones((dft_len+cp_len+cs_len-2*tail_len, 1)),
                   np.zeros((dft_len+cp_len+cs_len-2*tail_len, tail_len)))),
        np.hstack((np.zeros((tail_len, 1)), np.eye(tail_len)))))


def reduce_variable_rx(dft_len:int,  tail_len:int):
    """
    Method to remove repetitive samples from our optimization variable
    for the receiver window optimization.

    Parameters
    ----------
    dft_len : int
        Length of the DFT.
    tail_len : int
        Lenght of the receiver window tail.

    Returns
    -------
    reduce_var_mat : np.ndarray
        Matrix responsible for reducaing the optimization variable.
    """

    return np.vstack((
        np.hstack((np.ones((int(tail_len/2), 1)), -np.eye(int(tail_len/2)))),
        np.hstack((np.zeros((int(tail_len/2), 1)),
                   np.fliplr(np.eye(int(tail_len/2))))),
        np.hstack((np.ones((dft_len-tail_len, 1)),
                   np.zeros((dft_len-tail_len, int(tail_len/2))))),
        np.hstack((np.zeros((int(tail_len/2), 1)), np.eye(int(tail_len/2)))),
        np.hstack((np.ones((int(tail_len/2), 1)),
                   -np.fliplr(np.eye(int(tail_len/2)))))
    ))


def gen_constraints_tx(tail_len:int):
    """Method to generate constraints for the optimization of Tx.
    
    Parameters
    ----------
    tail_len : int
        Length of the window tail.

    Returns
    -------
    A : np.nadarray
        Linear equality constraints.
    b : np.ndarray
        Values for linear equality constraints.
    C : np.ndarray
        Linear inequality constraints.
    d : np.ndarray
        Values for linear inquality constraints.
    """

    A = np.hstack((np.array([1.], ndmin=2),
                   np.zeros((1, tail_len), dtype=np.float64)))
    b = np.array([1.], ndmin=2)
    C = np.hstack((np.zeros((tail_len, 1),dtype=np.float64),
                   np.eye(tail_len, dtype=np.float64)))
    d = np.ones((tail_len, 1), dtype=np.float64)

    return A, b, C, d


def gen_constraints_rx(tail_len:int):
    """Method to generate constriants for the optimization of Rx.
    
    Parameters
    ----------
    tail_len : int
        Length of the window tail.
    
    Returns
    -------
    A : np.ndarray
        Linear equality constraints.
    b : np.ndarray
        Values for linear equality constraints.
    C : np.ndarray
        Linear inequality constraints.
    d : np.ndarray
        Values for linear inequality constraints.
    """

    A = np.hstack((np.array([1.], ndmin=2), np.zeros((1, int(tail_len/2)))))
    b = np.array([1.], ndmin=2)
    C = np.vstack((
        np.hstack((np.zeros((int(tail_len/2), 1)), np.eye(int(tail_len/2)))),
        np.hstack((np.ones((int(tail_len/2), 1)), -np.eye(int(tail_len/2))))
    ))
    d = np.vstack((np.ones((int(tail_len/2), 1)),
                   .5*np.ones((int(tail_len/2), 1))))

    return A, b, C, d


def gen_constraints_tx_rx(tail_tx_len:int, tail_rx_len:int):
    """Method to reduce variable before for both ends optimization.
    
    Parameters
    ----------
    tail_tx_len : int
        Number of samples at the transmitter's tail.
    tail_rx_len : int
        Number of samples at the receiver's tail.
    
    Returns
    -------
    Q : np.ndarray
        Matrix with quadratic constrains.
    """

    n_samples = int((tail_rx_len/2 + 1)*(tail_tx_len+1))
    Q = np.zeros((n_samples, n_samples), dtype=np.float64)
    for n in range(1, int(1+tail_rx_len/2)):
        i = n*(tail_tx_len+1)
        for j in range(1, tail_tx_len+1):
            Q[0, i+j] = 1
            Q[i+j, 0] = 1
            Q[i, j] = -1
            Q[j, i] = -1
            # m = int(tail_tx_len*np.floor(i/tail_tx_len))
            # n = i - m
            # Q[0, i] = 1
            # Q[m, n] = -1
            # Q[n, m] = -1
            # Q[i, 0] = 1

    return Q


# EoF
