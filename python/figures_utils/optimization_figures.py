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


def gen_figures_opt(folder_path:str, dft_len:int, cp_list:list, sys_list:list):
    """Method to generate figures for optimization process.
    
    Parameters
    ----------
    folder_path : str
        Path to optimized windows.
    dft_len : int
        Number of bins in the DFT.
    cp_list : list
        List with cyclic prefix lengths.
    sys_list : list
        List of system designs.
    """
    tail_tx_fun = lambda x,y: x if y in ['CPW', 'WOLA', 'CPwtx', 'wtx'] else 0
    tail_rx_fun = lambda x,y: x if y in ['CPW', 'WOLA', 'CPwrx', 'wrx'] else 0
    data = np.zeros((len(cp_list), int(2*len(sys_list))), dtype=np.float64)
    for i, cp_len in enumerate(cp_list):
        for j, sys in enumerate(sys_list):
            win_tx, win_rx = read_window(folder_path, cp_len, sys, dft_len,
                                         tail_tx_fun(8, sys),
                                         tail_rx_fun(10, sys))
            print(win_tx)
            print(win_rx)


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
        cs_len = int(tail_tx + tail_rx/2) if sys_design == 'WOLA' else tail_tx
    win_tx = np.diagflat(reduce_variable_tx(dft_len, cp_len, cs_len, tail_tx)
                         @ win_data_tx)
    win_rx = np.diagflat(reduce_variable_rx(dft_len, tail_rx)@win_data_rx)

    return win_tx, win_rx


def read_cond_number(folder_path:str):
    """Method to read condition numbers."""

    folder_name = os.path.join(folder_path, 'condition_number')
    file_list = [f.path for f in os.scandir(folder_name) if
                 f.name.endswith('.npy')]
    cp_list = list(set([i.split('/').split('_')[1] for i in file_list]))
    sys_list = list(set([i.split('/').split('_')[0] for i in file_list]))
    reg_val = list(set([i.split()]))


def plot_interference():
    """Method to plot interference power."""


def plot_condition_number():
    """Method to plot condition number."""


# EoF
