"""interf_calc.py

Script with method to calculate the interference based on the window
and the OFDM system.

luizfelipe.coelho@smt.ufrj.br
Nov 11, 2023
"""


import numpy as np
from .transmitter import (gen_add_redundancy_matrix, gen_idft_matrix,
                          gen_rc_window_tx)
from .receiver import (gen_circ_shift_matrix, gen_dft_matrix,
                       gen_overlap_and_add_matrix, gen_rc_window_rx,
                       gen_rm_redundancy_matrix)
from .channel import gen_channel_tensor


def interf_power(sys_design:str, window_data:list, dft_len:int, cp_len:int,
                 tail_tx:int, tail_rx:int):
    """Method to calculate the interference power for optimized windows.
    
    Parameters
    ----------
    sys_design : str
        Name of the system.
    window_data : list
        List with window data.
    
    Returns
    -------
    interf_power : float
        Estimated interference power.
    """

    if sys_design in ['wtx', 'CPwtx']:
        cs_len = tail_tx if sys_design == 'wtx' else 0
        circ_shift = tail_tx if sys_design == 'CPwtx' else 0
        rm_len = cp_len if sys_design == 'wtx' else cp_len - tail_tx
        Vtx_rc = gen_rc_window_tx(dft_len, cp_len, cs_len, tail_tx)
        Vrx_rc = np.eye(dft_len+tail_rx)
    elif sys_design in ['wrx', 'CPwrx']:
        cs_len = int(tail_rx/2) if sys_design == 'wrx' else 0
        circ_shift = int(tail_rx/2) if sys_design == 'CPwrx' else 0
        rm_len = cp_len - tail_rx if sys_design == 'CPwrx' else \
            int(cp_len - tail_rx/2)
        Vtx_rc = np.eye(dft_len+cp_len+cs_len)
        Vrx_rc = gen_rc_window_rx(dft_len, tail_rx)
    elif sys_design in ['CPW', 'WOLA']:
        cs_len = int(tail_tx + tail_rx/2) if sys_design == 'WOLA' else tail_tx
        circ_shift = int(tail_rx/2) if sys_design == 'CPW' else 0
        rm_len = cp_len - tail_rx if sys_design == 'CPW' else \
            int(cp_len - tail_rx/2)
        Vtx_rc = gen_rc_window_tx(dft_len, cp_len, cs_len, tail_tx)
        Vrx_rc = gen_rc_window_rx(dft_len, tail_rx)
    Vtx_opt, Vrx_opt = window_data
    # Transmitter:
    W_inv = gen_idft_matrix(dft_len)
    Gamma = gen_add_redundancy_matrix(dft_len, cp_len, cs_len)
    # Channel:
    chann = np.load('channels/vehicularA.npy')
    channel_ir = np.mean(chann, axis=1)
    H_m = gen_channel_tensor(channel_ir, dft_len, tail_tx, tail_rx, rm_len,
                             cp_len, cs_len)
    # Transmitter:
    R = gen_rm_redundancy_matrix(dft_len, tail_rx, rm_len)
    P = gen_overlap_and_add_matrix(dft_len, tail_rx)
    K = gen_circ_shift_matrix(dft_len, circ_shift)
    W = gen_dft_matrix(dft_len)
    # Calculation:
    # 1 - Optimized:
    A0 = W@K@P@Vrx_opt@R@H_m[:, :, 0]@Vtx_opt@Gamma@W_inv
    AICI1 = A0 - np.diagflat(np.diag(A0))
    Am = W@K@P@Vrx_opt@R@np.moveaxis(H_m[:, :, 1:], [2], [0])@Vtx_opt \
        @ Gamma@W_inv
    PISI = np.zeros((dft_len, dft_len))
    for idx in range(Am.shape[0]):
        PISI += (Am[idx, :, :] @ np.conjugate(np.transpose(Am[idx, :, :]))).real
    PICI1 = (AICI1 @ np.conjugate(np.transpose(AICI1))).real
    # import pdb; pdb.set_trace()
    P_tot = np.diag(PISI+PICI1)
    # import pdb;pdb.set_trace()
    # 2 - RC:
    A0_rc = W@K@P@Vrx_rc@R@H_m[:, :, 0]@Vtx_rc@Gamma@W_inv
    AICI1_rc = A0_rc - np.diagflat(np.diag(A0_rc))
    Am_rc = W@K@P@Vrx_rc@R@np.moveaxis(H_m[:, :, 1:], [2], [0])@Vtx_rc \
        @ Gamma@W_inv
    PISI_rc = np.zeros((dft_len, dft_len))
    for idx in range(Am_rc.shape[0]):
        PISI_rc += (Am_rc[idx, :, :] @ np.conjugate(np.transpose(Am_rc[idx, :, :]))).real
    PICI1_rc = (AICI1_rc @ np.conjugate(np.transpose(AICI1_rc))).real
    P_rc = np.diag(PISI_rc+PICI1_rc)

    return P_tot, P_rc


# EoF
