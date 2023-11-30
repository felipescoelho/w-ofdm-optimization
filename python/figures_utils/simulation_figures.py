"""simulation_figures.py

Script with functions to create figures for analyzing the simulation
results. Mainly SER, Symbol Error Rate.

luizfelipe.coelho@smt.ufrj.br
Nov 27, 2023
"""


import matplotlib.pyplot as plt
import numpy as np
import os


def gen_figures_sim(folder_path:str, cp_list:list, sys_list:list,
                    snr_arr:np.ndarray):
    """Method to generate figures for the simulation.
    
    Parameters
    ----------
    folder_path : str
        Path to simulation results.
    dft_len : int
        Number of bins in the DFT.
    cp_list : list
        List with cyclic prefix lengths.
    sys_list : list
        List of system designs.
    snr_arr : np.ndarray
        Array with values for the SNR.
    """

    data_opt = np.zeros((len(snr_arr), len(sys_list)),
                    dtype=np.float64)
    data_opt_corrected = np.zeros((len(snr_arr), len(sys_list)),
                                  dtype=np.float64)
    data_rc = np.zeros((len(snr_arr), len(sys_list)),
                       dtype=np.float64)
    data_rc_corrected = np.zeros((len(snr_arr), len(sys_list)),
                                dtype=np.float64)
    tail_tx = 8
    tail_rx = 10
    for cp_len in cp_list:
        cp_ser = read_results_ser(folder_path, cp_len, 'CP')
        for j, sys in enumerate(sys_list):
            if sys == 'CPwrx':
                cp_corrected = cp_len
            elif sys == 'wrx':
                cp_corrected = cp_len - int(tail_rx/2)
            else:
                if sys in ['wtx', 'CPW']:
                    cs_len = tail_tx
                elif sys == 'CPwtx':
                    cs_len = 0
                else:
                    cs_len = tail_tx + int(tail_rx/2)
                cp_corrected = cp_len-cs_len+tail_tx
            opt_ser, rc_ser = read_results_ser(folder_path, cp_len, sys)
            opt_ser_corrected, rc_ser_corrected = read_results_ser(
                folder_path, cp_corrected, sys
            )
            data_opt[:, j] = opt_ser
            data_opt_corrected[:, j] = opt_ser_corrected
            data_rc[:, j] = rc_ser
            data_rc_corrected[:, j] = rc_ser_corrected
        plot_ser_single_cp(data_opt, data_opt_corrected, data_rc,
                           data_rc_corrected, cp_ser, snr_arr, sys_list)
        

def read_results_ser(folder_path:str, cp_len:int, sys:str):
    """Reads the saved results."""

    if sys == 'CP':
        file_path = os.path.join(folder_path, 'ser', f'{sys}_{cp_len}.npy')

        return np.load(file_path)

    file_path_opt = os.path.join(folder_path, 'ser', f'opt_{sys}_{cp_len}.npy')
    file_path_rc = os.path.join(folder_path, 'ser', f'rc_{sys}_{cp_len}.npy')

    return np.load(file_path_opt), np.load(file_path_rc)

def plot_ser_single_cp(data_opt:np.ndarray, data_opt_corrected:np.ndarray,
                       data_rc:np.ndarray, data_rc_corrected:np.ndarray,
                       data_cp:np.ndarray, snr_arr:np.ndarray, sys_list:list):
    
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                  'tab:purple', 'tab:brown']
    markers = ['P', '*', 's', 'd', 'o', 'X']

    fig0 = plt.figure()
    ax0 = fig0.add_subplot(1, 1, 1)
    ax0.plot(snr_arr, data_cp, marker='h', mfc='none', ms=10, label='CP',
             c='k')
    for idx, sys in enumerate(sys_list):
        ax0.plot(snr_arr, data_opt[:, idx], marker=markers[idx], mfc='none',
                 ms=10, label=sys, c=color_list[idx])
        ax0.plot(snr_arr, data_rc[:, idx], '--', marker=markers[idx],
                 mfc='none', ms=10, c=color_list[idx])
    ax0.legend()
    ax0.set_yscale('log')
    ax0.set_xlabel('SNR, dB')
    ax0.set_ylabel('SER')
    ax0.grid()

    fig1 = plt.figure()
    ax0 = fig1.add_subplot(1, 1, 1)
    ax0.plot(snr_arr, data_cp, marker='h', mfc='none', ms=10, label='CP',
             c='k')
    for idx, sys in enumerate(sys_list):
        ax0.plot(snr_arr, data_opt_corrected[:, idx], marker=markers[idx], mfc='none',
                 ms=10, label=sys, c=color_list[idx])
        ax0.plot(snr_arr, data_rc_corrected[:, idx], '--', marker=markers[idx],
                 mfc='none', ms=10, c=color_list[idx])
    ax0.legend()
    ax0.set_yscale('log')
    ax0.set_xlabel('SNR, dB')
    ax0.set_ylabel('SER')
    ax0.grid()

    plt.show()

    

# EoF
