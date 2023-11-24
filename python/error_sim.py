import numpy as np
from optimization_tools import reduce_variable_rx, reduce_variable_tx
from numba import njit
from ofdm_utils import (
    gen_idft_matrix, gen_add_redundancy_matrix, gen_rm_redundancy_matrix,
    gen_overlap_and_add_matrix, gen_dft_matrix, gen_circ_shift_matrix,
    gen_channel_tensor, gen_rc_window_rx, gen_rc_window_tx
)

@njit(fastmath=False)
def my_fun(B, C):
    n_rows, n_cols = B.shape
    M_mat = np.zeros((n_cols, n_cols), dtype=np.complex128)
    for i in range(n_cols):
        for j in range(n_cols):
            for m in range(n_rows):
                for n in range(n_rows):
                    M_mat[i,j] += C[i, n]* B[m, i] * np.conj(C[j, n]) \
                        * np.conj(B[m, j])
                    
    return M_mat.real


if __name__ == '__main__':

    # Definitions:
    dft_len = 256
    cp_len = 16
    cs_len = 8
    tail_tx = 10
    tail_rx = 0
    rm_len = 0

    red_mat = reduce_variable_tx(dft_len, cp_len, cs_len, tail_tx)

    W = gen_dft_matrix(dft_len)
    K = gen_circ_shift_matrix(dft_len, 0)
    P = gen_overlap_and_add_matrix(dft_len, tail_rx)
    R = gen_rm_redundancy_matrix(dft_len, tail_rx, rm_len)
    Gamma = gen_add_redundancy_matrix(dft_len, cp_len, cs_len)
    W_inv = gen_idft_matrix(dft_len)

    random_matrix = (np.random.randn(dft_len+tail_rx, dft_len+cp_len+cs_len) \
        + 1j*np.random.randn(dft_len+tail_rx, dft_len+cp_len+cs_len))*1e-3
    
    B = W@K@P@R@random_matrix
    C = Gamma@W_inv
    # Experiment calculation:
    M_mat1 = my_fun(B, C)
    M_mat2 = (W@K@P@R@random_matrix@Gamma@W_inv) \
        @ np.conj(W@K@P@R@random_matrix@Gamma@W_inv).T
    M_mat3 = B@C @ np.conj(B@C).T

    Q_mat1 = red_mat.T @ M_mat1 @ red_mat
    # Q_mat2 = np.trace(M_mat)
    # Q_mat3 = red_mat.T @ np.diag(np.diag(M_mat3.real)) @ red_mat

    x = np.ones((6, 1), dtype=np.float64)

    print((x.T @ Q_mat1 @ x)[0])
    print(np.trace(M_mat2).real)
    print(np.trace(M_mat3).real)

    # First should be really different.


# EoF
