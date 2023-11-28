"""error_analysis.py

luizfelipe.coelho@smt.ufrj.br
Nov 25, 2023
"""


import numpy as np
from optimization_tools import reduce_variable_rx, reduce_variable_tx
from numba import njit
from ofdm_utils import (
    gen_idft_matrix, gen_add_redundancy_matrix, gen_rm_redundancy_matrix,
    gen_overlap_and_add_matrix, gen_dft_matrix, gen_circ_shift_matrix,
    gen_channel_tensor, gen_rc_window_rx, gen_rc_window_tx
)


def test_win_tx(dft_len:int, cp_len:int, cs_len:int, tail_tx:int):
    
    return np.eye(dft_len+cp_len+cs_len)


def test_win_rx(dft_len:int, tail_rx:int):

    return np.diag(np.hstack((
        np.zeros((int(tail_rx/2),)),
        np.ones((dft_len,)),
        np.zeros((int(tail_rx/2),))
    )))


@njit(fastmath=True)
def my_fun(B, C):
    n_rows, n_cols = B.shape
    M_mat = np.zeros((n_cols, n_cols), dtype=np.complex128)
    for i in range(n_cols):
        for j in range(n_cols):
            for m in range(n_rows):
                for n in range(n_rows):
                    if m != n:
                        M_mat[i,j] += C[i, n]* B[m, i] * np.conj(C[j, n]) \
                            * np.conj(B[m, j])
                    
    return M_mat.real


def my_fun2(B_mat:np.ndarray, C_mat:np.ndarray):
    """Method to generate the ISI part of the cost function.
    
    Parameters
    ----------
    B_mat : np.ndarray
    C_mat : np.ndarray
    
    Returns
    -------
    out_mat : np.ndarray
    """

    out_mat = (C_mat @ np.conj(C_mat.T)) * (np.conj(B_mat.T) @ B_mat)

    return out_mat.real


@njit(fastmath=True)
def my_fun3(B_mat:np.ndarray, C_mat:np.ndarray, D_mat:np.ndarray,
            red_tx:np.ndarray, red_rx:np.ndarray):
    """"""
    dft_len, _ = B_mat.shape
    len_rx, len_tx = C_mat.shape
    _, len_red_rx = red_rx.shape
    _, len_red_tx = red_tx.shape
    n_samples = int(len_red_tx*len_red_rx)
    mm = np.zeros((n_samples, int(dft_len*dft_len)), dtype=np.complex128)
    for l in range(dft_len):
        for k in range(dft_len):
            if l == k:
                continue
            M_mat = np.zeros((len_rx, len_tx), dtype=np.complex128)
            for m in range(len_rx):
                for n in range(len_tx):
                    M_mat[m, n] = B_mat[l, m]*C_mat[m, n]*D_mat[n, k]
            M_prim = red_rx.T @ M_mat @ red_tx
            m = M_prim.flatten()
            mm[:, int(dft_len*l + k)] = m

    Q1_mat = mm @ np.conj(mm.T)

    return Q1_mat.real


@njit(fastmath=True)
def my_fun4(B_mat:np.ndarray, C_mat:np.ndarray, D_mat:np.ndarray,
                red_tx:np.ndarray, red_rx:np.ndarray):
    """
    Method to generate Q1 using numba (which is way faster).
    
    Parameters
    ----------
    B_mat : np.ndarray
        Array with B matrix.
    C_mat : np.ndarray
        Array with C matrix.
    D_mat : np.ndarray
        Array with D matrix.
    red_tx : np.ndarray
        Array with reduction matrix for the transmitter.
    red_rx : np.ndarray
        Array with reduction matrix for the receiver.
    """

    dft_len, _ = B_mat.shape
    len_rx, len_tx = C_mat.shape
    _, len_red_rx = red_rx.shape
    _, len_red_tx = red_tx.shape
    n_samples = int(len_red_tx*len_red_rx)
    mm = np.zeros((n_samples, int(dft_len*dft_len)), dtype=np.complex128)
    for i in range(dft_len):
        for j in range(dft_len):
            M_mat = np.zeros((len_rx, len_tx), dtype=np.complex128)
            for m in range(len_rx):
                for n in range(len_tx):
                    M_mat[m, n] = B_mat[i, m] * D_mat[n, j]
            M_mat *= C_mat
            M_prim = red_rx.T @ M_mat @ red_tx
            m = M_prim.flatten()
            mm[:, int(dft_len*i + j)] = m
    Q2_mat = mm @ np.conj(mm.T)

    return Q2_mat.real


if __name__ == '__main__':

    # Definitions:
    dft_len = 256
    cp_len = 16
    cs_len = 8
    tail_tx = 10
    tail_rx = 0
    rm_len = 0

    red_mat_tx = reduce_variable_tx(dft_len, cp_len, cs_len, tail_tx)
    red_mat_rx = reduce_variable_rx(dft_len, tail_rx)
    W = gen_dft_matrix(dft_len)
    K = gen_circ_shift_matrix(dft_len, 0)
    P = gen_overlap_and_add_matrix(dft_len, tail_rx)
    R = gen_rm_redundancy_matrix(dft_len, tail_rx, rm_len)
    Gamma = gen_add_redundancy_matrix(dft_len, cp_len, cs_len)
    W_inv = gen_idft_matrix(dft_len)

    chann = np.load('channels/vehicularA.npy')
    chann_avg = np.mean(chann, 1)
    H_mat = gen_channel_tensor(chann_avg, dft_len, tail_tx, tail_rx, rm_len,
                               cp_len, cs_len)
    # Calculations:
    B = W@K@P@R@H_mat[:, :, 0]
    C = Gamma@W_inv
    D = W@K@P@R@np.sum(H_mat[:, :, 1:], 2)
    M_mat1 = my_fun(B, C)
    M_mat2 = my_fun2(D, C)
    A0 = B@C
    AICI1 = A0 - np.diag(np.diag(A0))
    PICI1 = np.trace(AICI1 @ np.conj(AICI1.T)).real
    Am = W@K@P@R@np.moveaxis(H_mat[:, :, 1:], [2], [0])@Gamma@W_inv
    PISI_aux = np.zeros((dft_len, dft_len), dtype=np.float64)
    for idx in range(Am.shape[0]):
        A_idx = Am[idx, :, :]
        PISI_aux += (A_idx @ np.conj(A_idx.T)).real
    PISI = np.trace(PISI_aux)

    # Both ends:
    B2 = W@K@P
    C2 = R@H_mat[:, :, 0]
    C2_2 = R@np.sum(H_mat[:, :, 1:], 2)
    D2 = Gamma@W_inv
    win_rx = test_win_rx(dft_len, tail_rx)
    win_tx = test_win_tx(dft_len, cp_len, cs_len, tail_tx)
    
    M_mat3 = my_fun3(B2, C2, D2, red_mat_tx.astype(np.complex128),
                     red_mat_rx.astype(np.complex128))
    M_mat4 = my_fun4(B2, C2_2, D2, red_mat_tx.astype(np.complex128),
                     red_mat_rx.astype(np.complex128))
    A0_2 = B2 @ win_rx @ C2 @ win_tx @ D2
    AICI1_2 = A0_2 - np.diag(np.diag(A0_2))
    PICI1_2 = np.trace(AICI1_2 @ np.conj(AICI1_2.T)).real
    Am_2 = W@K@P@win_rx@R@np.moveaxis(H_mat[:, :, 1:], [2], [0])@win_tx \
        @Gamma@W_inv
    PISI_2aux = np.zeros((dft_len, dft_len), dtype=np.float64)
    for idx in range(Am_2.shape[0]):
        A_2idx = Am_2[idx, :, :]
        PISI_2aux += (A_2idx @ np.conj(A_2idx.T)).real
    PISI_2 = np.trace(PISI_2aux)
    x = np.ones((dft_len+cp_len+cs_len, 1), dtype=np.float64)
    x_rx = np.ones((int(tail_rx/2)+1, 1), dtype=np.float64)
    x_tx = np.ones((tail_tx+1, 1), dtype=np.float64)
    K = int((tail_rx/2 + 1) * (tail_tx+1))
    xx = np.zeros((K, 1), dtype=np.float64)
    for k in range(K):
        xx[k] = x_rx[int(np.floor(k/(tail_tx + 1)))] \
            * x_tx[int(k - (tail_tx + 1)*np.floor(k/(tail_tx + 1)))]

    y1 = (x.T @ M_mat1 @ x)[0]
    y2 = PICI1.copy()

    y3 = (x.T @ M_mat2 @ x)[0]
    y4 = PISI.copy()
    
    y5 = (xx.T @ M_mat3 @ xx)[0]
    y6 = PICI1_2.copy()

    y7 = (xx.T @ M_mat4 @ xx)[0]
    y8 = PISI_2.copy()

    # Single window:
    e_abs = abs(y1-y2)
    e_rltv = e_abs/y2
    e_abs2 = abs(y3-y4)
    e_rltv2 = e_abs2/y4
    # Both windows:
    e_abs3 = abs(y5-y6)
    e_rltv3 = e_abs3/y6
    e_abs4 = abs(y7-y8)
    e_rltv4 = e_abs4/y8

    print('Functions of single window:')
    print('\n')
    print('PICI1 Error:')
    print(f'Absolute error: {e_abs[0]}.')
    print(f'Relative error: {100*e_rltv[0]:.2f}%.')
    print('\n')
    print('PICI2 + PISI Error:')
    print(f'Absolute error: {e_abs2[0]}.')
    print(f'Relative error: {100*e_rltv2[0]:.2f}%.')
    print('\n')
    print('---------------------------------')
    print('Functions of two windows:')
    print('\n')
    print('PICI1 Error:')
    print(f'Absolute error: {e_abs3[0]}.')
    print(f'Relative error: {100*e_rltv3[0]:.2f}%.')
    print('\n')
    print('PICI2 + PISI Error:')
    print(f'Absolute error: {e_abs4[0]}.')
    print(f'Relative error: {100*e_rltv4[0]:.2f}%.')



# EoF
