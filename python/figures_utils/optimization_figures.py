"""optimization_figures.py

Script with functions to create figures for analyzing the optimization
process.

luizfelipe.coelho@smt.ufrj.br
Oct 15, 2023
"""


import matplotlib.pyplot as plt
import numpy as np
import os
from optimization_tools import reduce_variable_tx, reduce_variable_rx
from ofdm_utils import interf_power


plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'helvet'
})


def gen_figures_opt(folder_path:str, figures_path:str, dft_len:int,
                    cp_list:list, sys_list:list):
    """Method to generate figures for optimization process.
    
    Parameters
    ----------
    folder_path : str
        Path to optimized windows.
    figures_path : str
        Path to save figures.
    dft_len : int
        Number of bins in the DFT.
    cp_list : list
        List with cyclic prefix lengths.
    sys_list : list
        List of system designs.
    """

    fig_path = os.path.join(figures_path, 'interference_power')
    fig_path_win = os.path.join(figures_path, 'window_tail')
    os.makedirs(fig_path_win, exist_ok=True)
    os.makedirs(fig_path, exist_ok=True)
    tail_tx_fun = lambda x,y: x if y in ['CPW', 'WOLA', 'CPwtx', 'wtx'] else 0
    tail_rx_fun = lambda x,y: x if y in ['CPW', 'WOLA', 'CPwrx', 'wrx'] else 0
    data = np.zeros((dft_len, len(cp_list), len(sys_list)), dtype=np.float64)
    data_rc = np.zeros((dft_len, len(cp_list), len(sys_list)), dtype=np.float64)
    data_cp = np.zeros((dft_len, len(cp_list)), dtype=np.float64)
    for i, cp_len in enumerate(cp_list):
        for j, sys in enumerate(sys_list):
            win_tx, win_rx = read_window(folder_path, cp_len, sys, dft_len,
                                         tail_tx_fun(8, sys),
                                         tail_rx_fun(10, sys))
            plot_window([win_tx, win_rx], sys, cp_len, fig_path_win)
            P_opt, P_rc = interf_power(sys, [win_tx, win_rx], dft_len, cp_len,
                                       tail_tx_fun(8, sys), tail_rx_fun(10, sys))    
            data[:, i, j] = P_opt
            data_rc[:, i, j] = P_rc
        data_cp[:, i] = interf_power('CP', None, dft_len, cp_len, None, None)
    plot_interference(data, data_rc, data_cp, cp_list, sys_list, fig_path)


def read_window(folder_path:str, cp_len:int, sys_design:str, dft_len:int,
                tail_tx:int, tail_rx:int):
    """Method to recover/read windows for the different systems."""
 
    file_path = os.path.join(folder_path, f'{sys_design}_{cp_len}.npy')
    if sys_design in ['wtx', 'CPwtx']:
        win_data_tx = np.load(file_path)
        win_data_rx = np.array([1], ndmin=2, dtype=np.float64)
        cs_len = tail_tx if sys_design == 'wtx' else 0
    elif sys_design in ['wrx', 'CPwrx']:
        win_data_tx = np.array([1], ndmin=1)
        win_data_rx = np.load(file_path)
        cs_len = int(tail_rx/2) if sys_design == 'wrx' else 0
    elif sys_design in ['WOLA', 'CPW']:
        win_data = np.load(file_path)
        win_data_tx = win_data[:tail_tx+1].reshape(tail_tx+1, 1)
        win_data_rx = win_data[tail_tx+1:].reshape(int(tail_rx/2+1), 1)
        cs_len = int(tail_tx + tail_rx/2) if sys_design == 'CPW' else tail_tx
    win_tx = np.diagflat(reduce_variable_tx(dft_len, cp_len, cs_len, tail_tx)
                         @ win_data_tx)
    win_rx = np.diagflat(reduce_variable_rx(dft_len, tail_rx)@win_data_rx)

    return win_tx, win_rx


def plot_window(win_data, sys, cp_len, fig_path):
    """I assume beta = 8, delta = 10."""

    width = 5.93 / 3.1
    height = width * ((1 + 5**.5)/2)
    font_size = 12

    fig0 = plt.figure(figsize=(width, height))
    ax0 = fig0.add_subplot(1, 1, 1)
    if sys in ['wtx', 'CPwtx']:
        win = win_data[0]
        tail_tx = np.flipud(np.diag(win)[:8])
        ax0.plot(tail_tx, c='tab:blue')
    elif sys in ['wrx', 'CPwrx']:
        win = win_data[1]
        tail_rx = np.flipud(np.diag(win)[:10])
        ax0.plot(tail_rx, c='tab:orange')
    else:
        win_tx, win_rx = win_data
        tail_tx = np.flipud(np.diag(win_tx)[:8])
        tail_rx = np.flipud(np.diag(win_rx)[:10])
        ax0.plot(tail_tx, label='Tx', c='tab:blue')
        ax0.plot(tail_rx, label='Rx', c='tab:orange')
        ax0.legend(fontsize=font_size)

    ax0.set_ylabel('Value', fontsize=font_size)
    ax0.set_xlabel('Sample', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid()
    fig0.tight_layout()
    path = os.path.join(fig_path, f'{sys}_{cp_len}.eps')
    fig0.savefig(path, bbox_inches='tight')

    plt.close()


def plot_interference(data:np.ndarray, data_rc:np.ndarray, data_cp:np.ndarray,
                      cp_list:list, sys_list:list, fig_path:str):
    """Method to plot interference power.
    
    Parameters
    ----------
    data : np.ndarray
        Interference power data.
    cp_list : list
        List with CP length.
    sys_list : list
        List of each system.
    """
    
    data_sum = np.sum(data, axis=0)
    data_rc_sum = np.sum(data_rc, axis=0)
    data_cp_sum = np.sum(data_cp, axis=0)
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                  'tab:purple', 'tab:brown']
    markers = ['P', '*', 's', 'd', 'o', 'X']
    
    width = 5.93
    height = width/((1 + 5**.5)/2)
    font_size = 12

    fig0 = plt.figure(figsize=(width, height))
    ax0 = fig0.add_subplot(1, 1, 1)
    ax0.plot(cp_list, data_cp_sum, marker='h', mfc='none', ms=10, label='CP',
             c='k')
    for idx, sys in enumerate(sys_list):
        ax0.plot(cp_list, data_sum[:, idx], marker=markers[idx], mfc='none',
                 ms=10, label=sys, c=color_list[idx])
        ax0.plot(cp_list, data_rc_sum[:, idx], '--', marker=markers[idx],
                 mfc='none', ms=10, c=color_list[idx])
    ax0.legend(fontsize=font_size)
    ax0.set_yscale('log')
    ax0.set_ylim([1e-12, 1])
    ax0.set_xlim([10, 32])
    ax0.set_xticks([10, 14, 18, 22, 26, 30])
    ax0.set_xlabel('CP Length, $\mu$', fontsize=font_size)
    ax0.set_ylabel('Total Interf. Power', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid()
    fig0.tight_layout()

    fig1 = plt.figure(figsize=(width, height))
    ax0 = fig1.add_subplot(1, 1, 1)
    tail_tx = 8
    tail_rx = 10
    ax0.plot(cp_list, data_cp_sum, marker='h', mfc='none', ms=10, label='CP',
             c='k')
    for idx, sys in enumerate(sys_list):
        if sys == 'CPwrx':
            cp_list_corrected = cp_list
        elif sys == 'wrx':
            cp_list_corrected = [cp + int(tail_rx/2) for cp in cp_list]
        else:
            if sys in ['wtx', 'WOLA']:
                cs = tail_tx
            elif sys == 'CPwtx':
                cs = 0
            else:
                cs = tail_tx + int(tail_rx/2)
            cp_list_corrected = [cp+cs-tail_tx for cp in cp_list]
        ax0.plot(cp_list_corrected, data_sum[:, idx], marker=markers[idx],
                 mfc='none',
                 ms=10, label=sys, c=color_list[idx])
        ax0.plot(cp_list_corrected, data_rc_sum[:, idx], '--', marker=markers[idx],
                 mfc='none', ms=10, c=color_list[idx])
    ax0.legend(fontsize=font_size)
    ax0.set_yscale('log')
    ax0.set_ylim([1e-12, 1])
    ax0.set_xlim([2, 32+int(tail_rx/2)])
    ax0.set_xlabel('Total Overhead, $\mu + \\rho - \\beta$', fontsize=font_size)
    ax0.set_ylabel('Total Interf. Power', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid()
    fig1.tight_layout()

    path0 = os.path.join(fig_path, 'by_cp_len.eps')
    path1 = os.path.join(fig_path, 'by_overhead.eps')
    fig0.savefig(path0, bbox_inches='tight')
    fig1.savefig(path1, bbox_inches='tight')

    plt.close()


# EoF
