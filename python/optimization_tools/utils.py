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
        Length of the window tail.
    
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


# EoF
