"""transmitter.py

Script with functions for the transmitter.

luizfelipe.coelho@smt.ufrj.br
Jul 22, 2023
"""


import numpy as np
import sympy as sym


def gen_add_redundancy_sym(idft_len:int, cp_len:int, cs_len:int):
    """Mathod to generate matrix that adds redundancy with sympy.
    
    Parameters
    ----------
    idft_len : int
        Length of IDFT.
    cp_len : int
        Length of cyclic prefix (CP).
    cs_len : int
        Length of cyclic suffix (CS).
        
    Returns
    -------
    add_redundancy_mat : sym.Matrix
        Matrix responsible for adding redundancy.
    """

    add_redundancy_mat = sym.Matrix((
        [sym.zeros(cp_len, idft_len-cp_len), sym.eye(cp_len)],
        [sym.eye(idft_len)],
        [sym.eye(cs_len), sym.zeros(cs_len, idft_len-cs_len)]
    ))

    return add_redundancy_mat


def gen_idft_sym(idft_len:int):
    """Method to generate IDFT matrix using sympy.

    Parameters
    ----------
    idft_len : int
        Length of IDFT.

    Returns
    -------
    idft_mat : sym.Matrix
        Matrix responsible for IDFT.
    """

    idft_mat = sym.Matrix(
        idft_len, idft_len, lambda i,j: sym.exp(sym.I*2*sym.pi*i*j/idft_len)
    )

    return idft_mat/idft_len


def gen_rc_window_tx_sym(idft_len:int, cp_len:int, cs_len:int, tail_len:int):
    """Method to generate the Raised Cosine window for transmitter using sympy.
    
    Parameters
    ----------
    idft_len : int
        Length of IDFT.
    cp_len : int
        Length of cyclic prefix.
    cs_len : int
        Length of cyclic suffix.
    tail_len : int
        Length of window tail.
    
    Returns
    -------
    rc_window : sym.Matrix
        Raised cosine window for the transmitter.
    """

    rc_fun = lambda x: sym.sin((sym.pi*x) / (2*tail_len))**2
    rc_axis = sym.Matrix(tail_len, 1, list(range(tail_len)))
    rc_tail = rc_axis.applyfunc(rc_fun)
    rc_window = (
        rc_tail.col_join(sym.ones(idft_len+cp_len+cs_len-2*tail_len, 1))
        ).col_join(rc_tail.permute_rows(list(np.arange(tail_len-1, -1, -1))))
    
    return rc_window


def gen_add_redundancy_matrix(idft_len, cp_len, cs_len):
    """Method to generate matrix that adds redudancy.

    Parameters
    ----------
    idft_len : int
        Length of IDFT.
    cp_len : int
        Length of cyclic prefix (CP).
    cs_len : int
        Length of cyclic suffix (CS).

    Returns
    -------
    add_redundancy_mat : np.ndarray
        Matrix responsible for adding redundancy.
    """

    cp_row = np.hstack((np.zeros((cp_len, idft_len-cp_len)), np.eye(cp_len)))
    cs_row = np.hstack((np.eye(cs_len), np.zeros((cs_len, idft_len-cs_len))))
    add_redundancy_mat = np.vstack((cp_row, np.eye(idft_len), cs_row))

    return add_redundancy_mat


def gen_idft_matrix(idft_len):
    """Method to generate IDFT matrix.
    
    Parameters
    ----------
    idft_len : int
        Length of IDFT.
    
    Returns
    -------
    idft_mat : np.ndarray
        Matrix responsible for IDFT.
    """

    idft_mat = np.zeros((idft_len, idft_len), dtype=np.complex128)
    for row_idx in range(idft_len):
        for col_idx in range(idft_len):
            idft_mat[row_idx, col_idx] = \
                np.exp(-1j*2*np.pi*row_idx*col_idx/idft_len)
        
    return idft_mat/idft_len


def gen_rc_window_tx(idft_len, cp_len, cs_len, tail_len):
    """Method to generate the Raised Cosine window for transmitter.
    
    Parameters
    ----------
    idft_len : int
        Length of IDFT.
    cp_len : int
        Length of cyclic prefix.
    cs_len : int
        Length of cyclic suffix.
    tail_len : int
        Length of window tail.
    
    Returns
    -------
    rc_window : np.ndarray
        Raised cosine window for the transmitter.
    """

    tail_axis = np.arange(-(tail_len+1)/2 + 1, (tail_len+1)/2 -1, 1)
    rc_tail = np.sin(np.pi/2 * (.5 + tail_axis/tail_len))

    rc_window = np.hstack((rc_tail,
                           np.ones((idft_len+cp_len+cs_len-(2*tail_len),)),
                           rc_tail[::-1]))
    
    return rc_window


# EoF
