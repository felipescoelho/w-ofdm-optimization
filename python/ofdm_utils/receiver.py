"""receiver.py

Script with functions for the receiver.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


import numpy as np


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


def gen_rc_window_rx(dft_len:int, tail_len:int):
    """Method to generate the Raised Cosine window for receiver.
    
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

    tail_axis = np.arange(-(tail_len+1)/2 + 1, (tail_len+1)/2, 1)
    rc_tail = np.sin(np.pi/2 * (.5 + tail_axis/tail_len))
    rc_window = np.hstack((rc_tail, np.ones((dft_len-tail_len,)), rc_tail[::-1]))

    return np.diag(rc_window)


def gen_overlap_and_add_matrix(dft_len:int, tail_len:int):
    """Method to generate the ovelap and add matrix.
    
    Parameters
    ----------
    dft_len : int
        Length of the DFT.
    tail_len : int
        Length of window tail.
    
    Returns
    -------
    overlap_and_add_matrix : np.ndarray
        Matrix responsible for the overlap-and-add process.
    """

    aux = int(tail_len/2)
    row0 = np.hstack((np.zeros((aux, aux)), np.eye(aux),
                      np.zeros((aux, dft_len-tail_len)), np.zeros((aux, aux)),
                      np.eye(aux)))
    row1 = np.hstack((np.zeros((dft_len-tail_len, tail_len)),
                      np.eye(dft_len-tail_len),
                      np.zeros((dft_len-tail_len, tail_len))))
    row2 = np.hstack((np.eye(aux), np.zeros((aux, aux)),
                      np.zeros((aux, dft_len-tail_len)), np.eye(aux),
                      np.zeros((aux, aux))))
    
    return np.vstack((row0, row1, row2))


def gen_circ_shift_matrix(dft_len:int, circ_shift_len:int):
    """Method to generate the circular shift matrix.
    
    Parameters
    ----------
    dft_len : int
        Length of the DFT.
    circ_shift_len : int
        Length of the circular shift.
    
    Returns
    -------
    circ_shift_matrix : np.ndarray
        Matrix responsible for the circular shift.
    """

    return np.vstack((
        np.hstack((np.zeros((dft_len-circ_shift_len, circ_shift_len)),
                   np.eye(dft_len-circ_shift_len))),
        np.hstack((np.eye(circ_shift_len),
                   np.zeros((circ_shift_len, dft_len-circ_shift_len))))
    ))


def gen_dft_matrix(dft_len:int):
    """Method to generate the DFT matrix.
    
    Parameters
    ----------
    dft_len : int
        Length of the DFT.
    
    Returns
    -------
    dft_mat : np.ndarray
        Matrix responsible for DFT.
    """

    dft_mat = np.zeros((dft_len, dft_len), dtype=np.complex128)
    for row_idx in range(dft_len):
        for col_idx in range(dft_len):
            dft_mat[row_idx, col_idx] = \
                np.exp(-1j*2*np.pi*row_idx*col_idx/dft_len)
            
    return dft_mat


if __name__ == '__main__':

    print(gen_overlap_and_add_matrix(15, 8))


# EoF
