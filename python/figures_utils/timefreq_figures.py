"""timefreq_figures.py

Script with functions to create figures for analyzing the out of band
radiation figures.

luizfelipe.coelho@smt.ufrj.br
Dec 3, 2023
"""


import matplotlib.pyplot as plt
import numpy as np
import os

# os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin/'
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'helvet'
})

def gen_figures_timefreq(folder_path:str, figures_path:str, cp_list:list,
                         sys_list:list):
    """Method to generate figures for our simulations.
    
    Parameters
    ----------
    folder_path : str
        Path to resutls.
    cp_list : list
        CP values to plot.
    sys_list : list
        List of systems to plot.
    """

    tail_tx = 8
    tail_rx = 10
    fig_path = os.path.join(figures_path, 'timefreq')
    os.makedirs(fig_path, exist_ok=True)
    for cp_len in cp_list:
        cp_data = read_results(folder_path, cp_len, 'CP')
        cp_tf = cp_data['X_est_cp']
        f_axis = cp_data['f_axis']
        for sys in sys_list:
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
            opt_timefreq, rc_timefreq = read_results(folder_path, cp_len, sys)
            opt_corrected, rc_corrected = read_results(folder_path,
                                                       cp_corrected, sys)
            # plot_timefreq_single_cp(
            #     cp_tf, opt_timefreq['X_est_opt'], opt_corrected['X_est_opt'],
            #     rc_timefreq['X_est_rc'], rc_corrected['X_est_rc'], f_axis, sys,
            #     cp_len, fig_path
            # )
            estimate_papr(
                cp_data['mf_band_cp'], opt_timefreq['mf_band_opt'],
                opt_corrected['mf_band_opt'], rc_timefreq['mf_band_rc'],
                rc_corrected['mf_band_rc'], sys, cp_len
            )


def read_results(folder_path:str, cp_len:int, sys:str):
    """"""

    if sys == 'CP':
        file_path = os.path.join(
            folder_path, 'timefreq', f'{sys}_{cp_len}.npz'
        )

        return np.load(file_path)
    
    filepath_opt = os.path.join(
        folder_path, 'timefreq', f'opt_{sys}_{cp_len}.npz'
    )
    filepath_rc = os.path.join(
        folder_path, 'timefreq', f'rc_{sys}_{cp_len}.npz'
    )
    
    return np.load(filepath_opt), np.load(filepath_rc)


def plot_timefreq_single_cp(cp_tf:np.ndarray, opt_tf:np.ndarray,
                            opt_corr:np.ndarray, rc_tf:np.ndarray,
                            rc_corr:np.ndarray, f_axis:np.ndarray, sys:str,
                            cp_len:int, fig_path:str):
    """"""

    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                  'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
                  'tab:olive', 'tab:cyan']
    width = 5.93
    height = width/((1 + 5**.5)/2)
    font_size = 12
    fig0 = plt.figure(figsize=(width, height))
    ax0 = fig0.add_subplot(1, 1, 1)
    ax0.plot(f_axis, 10*np.log10(cp_tf), label='CP', c=color_list[-3])
    ax0.plot(f_axis, 10*np.log10(rc_tf), label='RC', c=color_list[3])
    ax0.plot(f_axis, 10*np.log10(opt_tf), label='Opt.', c=color_list[-1])
    ax0.legend(fontsize=font_size)
    ax0.set_ylabel('PSD, dB', fontsize=font_size)
    ax0.set_xlabel('Frequency, Hz', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid()
    fig0.tight_layout()

    fig1 = plt.figure(figsize=(width, height))
    ax0 = fig1.add_subplot(1, 1, 1)
    ax0.plot(f_axis, 10*np.log10(cp_tf), label='CP', c=color_list[-3])
    ax0.plot(f_axis, 10*np.log10(rc_corr), label='RC', c=color_list[3])
    ax0.plot(f_axis, 10*np.log10(opt_corr), label='Opt.', c=color_list[-1])
    ax0.legend(fontsize=font_size)
    ax0.set_ylabel('PSD, dB', fontsize=font_size)
    ax0.set_xlabel('Frequency, Hz', fontsize=font_size)
    ax0.tick_params(axis='both', which='major', labelsize=font_size)
    ax0.grid()
    fig1.tight_layout()

    path0 = os.path.join(fig_path, f'{sys}_{cp_len}.eps')
    path1 = os.path.join(fig_path, f'{sys}_{cp_len}_corrected.eps')
    fig0.savefig(path0, bbox_inches='tight')
    fig1.savefig(path1, bbox_inches='tight')


def estimate_papr(cp_tf, opt_tf, opt_corr, rc_tf, rc_corr, sys, cp_len):
    """"""

    cp_papr = max(cp_tf) / np.mean(cp_tf)
    opt_papr = max(opt_tf) / np.mean(opt_tf)
    opt_corr_papr = max(opt_corr) / np.mean(opt_corr)
    rc_papr = max(rc_tf) / np.mean(rc_tf)
    rc_corr_papr = max(rc_corr) / np.mean(rc_corr)

    print(f'CP-OFDM - CP {cp_len}: {cp_papr}.')
    print(f'{sys} - CP {cp_len} - OPT: {opt_papr}.')
    print(f'{sys} - CP {cp_len} - RC: {rc_papr}.')
    print(f'{sys} - CP {cp_len} - OPT Corrected: {opt_corr_papr}.')
    print(f'{sys} - CP {cp_len} - RC Corrected: {rc_corr_papr}.')
    print('\n')


# EoF
