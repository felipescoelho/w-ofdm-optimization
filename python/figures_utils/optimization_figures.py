"""optimization_figures.py

Script with functions to create figures for analyzing the optimization
process.

luizfelipe.coelho@smt.ufrj.br
Oct 15, 2023
"""


import matplotlib.pyplot as plt
import numpy as np
import os


def gen_figures_opt(folder_path:str):
    """Method to generate figures for optimization process.
    
    Parameters
    ----------
    folder_path : str
        Path to optimized windows.
    """

    data, cp_list, sys_list = read_windows(folder_path)

def read_windows(folder_path:str):
    """Method to read windows from files."""

    file_list = [file.path for file in os.scandir(folder_path) if
                 file.name.endswith('.npy')]
    cp_list = list(set([i.split('/')[-1].split('_')[-1].split('.')[0] for i in
                        file_list]))
    sys_list = list(set([i.split('/')[-1].split('_')[0] for i in file_list]))
    data = np.zeros((len(cp_list), len(sys_list)))

    return data, cp_list, sys_list

def read_cond_number(folder_path:str):
    """Method to read condition numbers."""

def plot_interference():
    """Method to plot interference power."""

def plot_condition_number():
    """Method to plot condition number."""


# EoF
