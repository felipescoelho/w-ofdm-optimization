"""channel.py

Script with functions for the channel.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


import numpy as np
from numba import jit


@jit(nopython=True)
def gen_channel_tensor(channel_ir:np.ndarray, dft_len:int, tail_tx:int,
                       tail_rx:int, rm_len:int, cp_len:int, cs_len:int):
    """Method to generate channel tensor.
    
    Paramters
    ---------
    channel_ir: np.ndarray
        Array with channel's impulse response.
    dft_len : int
        Length of DFT.
    tail_tx : int
        Number of samples in transmitter's window tail.
    tail_rx : int
        Number of samples in receiver's window tail.
    rm_len : int
        Number of samples to skip before recovering useful block.
    cp_len : int
        Number of samples in cyclic prefix.
    cs_len : int
        Number of samples in cyclic suffix.

    Returns
    -------
    channel_tensor : np.ndarray
        Array with channel matrices for system.
    """

    depth = int(np.ceil((len(channel_ir)+tail_tx) / (dft_len+tail_rx+rm_len)))
    N_0 = dft_len+cp_len+cs_len-tail_tx
    channel_tensor = np.zeros((dft_len+tail_rx+rm_len, dft_len+cp_len+cs_len,
                               depth), dtype=np.complex128)
    for depth_idx in range(depth):
        for row_idx in range(dft_len+tail_rx+rm_len):
            for col_idx in range(dft_len+cp_len+cs_len):
                idx = depth_idx*N_0 + row_idx - col_idx
                if 0 <= idx <= len(channel_ir)-1:
                    channel_tensor[row_idx, col_idx, depth_idx] = \
                        channel_ir[depth_idx*N_0 + row_idx - col_idx]

    return channel_tensor


# EoF
