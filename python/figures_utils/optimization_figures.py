"""optimization_figures.py

Script with functions to create figures for analyzing the optimization
process.

luizfelipe.coelho@smt.ufrj.br
Oct 15, 2023
"""


import matplotlib.pyplot as plt
import os


def gen_figures():
    """"""


def read_windows(folder_path:str):
    """Method to read windows from files."""

    file_list = [file.name for file in os.scandir(folder_path) if
                 file.name.endswith('.npy')]


def plot_interference():
    """"""

def plot_condition_number():
    """"""


# EoF
