"""simulation_figures.py

Script with functions to create figures for analyzing the simulation
results. Mainly SER, Symbol Error Rate.

luizfelipe.coelho@smt.ufrj.br
Nov 27, 2023
"""


import matplotlib.pyplot as plt
import numpy as np
import os


plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'helvet'
})


def gen_figures_sim(folder_path:str, figures_path:str, cp_list:list,
                    sys_list:list, snr_arr:np.ndarray):
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

    fig_path = os.path.join(figures_path, 'ser')
    os.makedirs(fig_path, exist_ok=True)
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
                           data_rc_corrected, cp_ser, snr_arr, sys_list,
                           fig_path, cp_len)
        

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
                       data_cp:np.ndarray, snr_arr:np.ndarray, sys_list:list,
                       fig_path:str, cp_len:int):
    
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                  'tab:purple', 'tab:brown']
    markers = ['P', '*', 's', 'd', 'o', 'X']
    width = 5.93
    height = width/((1 + 5**.5)/2)
    font_size = 12
    fig0 = plt.figure(figsize=(width, height))
    ax0 = fig0.add_subplot(1, 1, 1)
    ax0.plot(snr_arr, data_cp, marker='h', mfc='none', ms=10, label='CP',
             c='k')
    for idx, sys in enumerate(sys_list):
        ax0.plot(snr_arr, data_opt[:, idx], marker=markers[idx], mfc='none',
                 ms=10, label=sys, c=color_list[idx])
        ax0.plot(snr_arr, data_rc[:, idx], '--', marker=markers[idx],
                 mfc='none', ms=10, c=color_list[idx])
    ax0.legend(fontsize=font_size)
    ax0.set_yscale('log')
    ax0.set_xlabel('SNR, dB', fontsize=font_size)
    ax0.set_ylabel('SER', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid(axis='both', which='both')
    fig0.tight_layout()

    fig1 = plt.figure(figsize=(width, height))
    ax0 = fig1.add_subplot(1, 1, 1)
    ax0.plot(snr_arr, data_cp, marker='h', mfc='none', ms=10, label='CP',
             c='k')
    for idx, sys in enumerate(sys_list):
        ax0.plot(snr_arr, data_opt_corrected[:, idx], marker=markers[idx],
                 mfc='none', ms=10, label=sys, c=color_list[idx])
        ax0.plot(snr_arr, data_rc_corrected[:, idx], '--', marker=markers[idx],
                 mfc='none', ms=10, c=color_list[idx])
    ax0.legend(fontsize=font_size)
    ax0.set_yscale('log')
    ax0.set_xlabel('SNR, dB', fontsize=font_size)
    ax0.set_ylabel('SER', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid(axis='both', which='both')
    fig1.tight_layout()

    if cp_len == 22:
        fig2 = plt.figure(figsize=(width, height))
        ax0 = fig2.add_subplot(1, 1, 1)
        x0, x1, y0, y1 = 40, 50, 5*1e-4, 2*1e-3
        axins = ax0.inset_axes(
            [.1, .1, .5, .5], xlim=(x0, x1), ylim=(y0, y1),
            yticklabels=[]
        )
        ax0.plot(snr_arr, data_cp, marker='h', mfc='none', ms=10, label='CP',
                c='k')
        axins.plot(snr_arr, data_cp, marker='h', mfc='none', ms=10, c='k')
        for idx, sys in enumerate(sys_list):
            ax0.plot(snr_arr, data_opt_corrected[:, idx], marker=markers[idx],
                     mfc='none', ms=10, label=sys, c=color_list[idx])
            ax0.plot(snr_arr, data_rc_corrected[:, idx], '--',
                     marker=markers[idx], mfc='none', ms=10, c=color_list[idx])
            axins.plot(snr_arr, data_opt_corrected[:, idx], marker=markers[idx],
                       mfc='none', ms=10, c=color_list[idx])
            axins.plot(snr_arr, data_rc_corrected[:, idx], '--',
                       marker=markers[idx], mfc='none', ms=10,
                       c=color_list[idx])
        ax0.indicate_inset_zoom(axins, edgecolor='black')
        ax0.legend(fontsize=10)
        ax0.set_yscale('log')
        ax0.set_xlabel('SNR, dB', fontsize=font_size)
        ax0.set_ylabel('SER', fontsize=font_size)
        ax0.tick_params(axis='both', which='major', labelsize=font_size)
        ax0.grid(axis='both', which='both')
        axins.set_yscale('log')
        axins.grid(axis='both', which='both')
        axins.set_yticklabels([], minor=True)
        plt.setp(axins.get_xticklabels(), backgroundcolor='white')
        plt.setp(axins.get_yticklabels(), backgroundcolor='white')
        fig2.tight_layout()
        path2 = os.path.join(fig_path, f'{cp_len}_corrected_zoom.eps')
        fig2.savefig(path2, bbox_inches='tight')

    path0 = os.path.join(fig_path, f'{cp_len}.eps')
    path1 = os.path.join(fig_path, f'{cp_len}_corrected.eps')

    fig0.savefig(path0, bbox_inches='tight')
    fig1.savefig(path1, bbox_inches='tight')

    plt.close()

    

# EoF
