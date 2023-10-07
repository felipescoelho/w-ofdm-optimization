"""receiver.py

Script with functions for the receiver.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


import numpy as np
import sympy as sym


def gen_rm_redundancy_matrix(dft_len:int, tail_len:int, rm_len:int):
    """Method to generate the matrix that recovers the useful block.
    
    Parameters
    ----------
    idft_len : int
        Length od IDFT.
    tail_len : int
        Length of tail for receiver's window.
    rm_len : int
        Number of samples to skip before useful block.
    
    Returns
    -------
    rm_redundancy_mat : np.ndarray
    """

    rm_redundancy_mat = np.hstack((np.zeros((dft_len+tail_len, rm_len)),
                                   np.eye(dft_len+tail_len)))
    
    return rm_redundancy_mat


def gen_rc_window_rx_vector(dft_len:int, tail_len:int):
    """Method to generate the Raised Cosine window for receiver using sympy.
    
    Parameters
    ----------
    dft_len : int
        Length of DFT.
    tail_len : int
        Length of window tail.
    
    Returns
    -------
    rc_window : np.ndarray
        Raised cosine window for the receiver.
    """

    rc_axis = np.arange(0, 1, tail_len)
    rc_fun = lambda x: np.sin((sym.pi*x) / (2*tail_len))**2
    rc_tail = rc_axis.applyfunc(rc_fun)
    rc_window = (
        rc_tail.col_join(sym.ones(dft_len-tail_len, 1))
    ).col_join(rc_tail.permute_rows(list(np.arange(tail_len-1, -1, -1))))

    return rc_window


def gen_rm_redundancy_sym(dft_len:int, tail_len:int, rm_len:int):
    """Method to generate the matrix that recovers the useful block.
    
    Parameters
    ----------
    idft_len : int
        Length of DFT.
    tail_len : int
        Length of tail for receiver's window.
    rm_len : int
        Number of sampes to skip before useful block.
        
    Returns
    -------
    rm_redudnacy_mat : sym.Matrix
    """

    rm_redudancy_mat = sym.zeros(dft_len+tail_len, rm_len).row_join(
        sym.eye(dft_len+tail_len)
    )

    return rm_redudancy_mat


def gen_rc_window_rx_sym(dft_len:int, tail_len:int):
    """Method to generate the Raised Cosine window for receiver using sympy.
    
    Parameters
    ----------
    dft_len : int
        Length of DFT.
    tail_len : int
        Length of window tail.
    
    Returns
    -------
    rc_window : sym.Matrix
        Raised cosine window for the receiver.
    """

    rc_fun = lambda x: sym.sin((sym.pi*x) / (2*tail_len))**2
    rc_axis = sym.Matrix(tail_len, 1, list(range(tail_len)))
    rc_tail = rc_axis.applyfunc(rc_fun)
    rc_window = (
        rc_tail.col_join(sym.ones(dft_len-tail_len, 1))
    ).col_join(rc_tail.permute_rows(list(np.arange(tail_len-1, -1, -1))))

    return rc_window


def gen_overlap_and_add_sym(dft_len:int, tail_len:int):
    """Method to generate the overlap and add matrix using sympy.
    
    Parameters
    ----------
    dft_len : int
        Length of DFT.
    tail_len : int
        Length of widnow tail.

    Returns
    -------
    overlap_add_mat : sym.Matrix
    """

    overlap_add_mat = sym.Matrix((
        [sym.zeros(tail_len/2, tail_len/2), sym.eye(tail_len/2),
         sym.zeros(tail_len/2, dft_len-tail_len),
         sym.zeros(tail_len/2, tail_len/2), sym.eye(tail_len/2)],
        [sym.zeros(dft_len-tail_len, tail_len), sym.eye(dft_len-tail_len),
         sym.zeros(dft_len-tail_len, tail_len)],
        [sym.eye(tail_len/2), sym.zeros(tail_len/2, tail_len/2),
         sym.zeros(tail_len/2, dft_len-tail_len), sym.eye(tail_len/2),
         sym.zeros(tail_len/2, tail_len/2)]
    ))

    return overlap_add_mat


def gen_circ_shift_sym(dft_len:int, shift_len:int):
    """Method to generate circular shift matrix.
    
    Parameters
    ----------
    dft_len : int
        Length of DFT.
    shift_len : int
        Number of samples in circular shift.
        
    Returns
    -------
    circ_shift_mat : sym.Matrix
        Matrix responsible for circular shift.
    """

    circ_shift_mat = sym.zeros(dft_len-shift_len, shift_len).row_join(
        sym.eye(dft_len-shift_len)
    ).col_join(sym.eye(shift_len).row_join(
        sym.zeros(shift_len, dft_len-shift_len)
    ))

    return circ_shift_mat


def gen_dft_sym(dft_len:int):
    """Method to generate DFT matrix using sympy.
    
    Parameters
    ----------
    dft_len : int
        Length of DFT.

    Returns
    -------
    dft_mat : sym.Matrix
        Matrix responsible for DFT.
    """

    dft_mat = sym.Matrix(
        dft_len, dft_len, lambda i,j: sym.exp(-sym.I*2*sym.pi*i*j/dft_len)
    )

    return dft_mat


def gen_rm_redudancy_matrix(dft_len:int, tail_len:int, rm_len:int):
    """Method to generate the matrix that recovers the useful block.
    
    Parameters
    ----------
    idft_len : int
        Length of DFT.
    tail_len : int
        Length of tail for receiver's window.
    rm_len : int
        Number of sampes to skip before useful block.
        
    Returns
    -------
    rm_redudnacy_mat : sym.Matrix
    """

    rm_redundancy_mat = np.hstack((np.zeros((dft_len+tail_len, rm_len)),
                                   np.eye(rm_len)))
    
    return rm_redundancy_mat


# EoF
