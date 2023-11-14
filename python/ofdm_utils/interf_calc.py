"""interf_calc.py

Script with method to calculate the interference based on the window
and the OFDM system.

luizfelipe.coelho@smt.ufrj.br
Nov 11, 2023
"""


import numpy as np
from optimization_tools import reduce_variable_tx, reduce_variable_rx


def interf_power_opt(sys_design:str, window_data:list):
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
        opt_window = reduce_variable_tx()
    elif sys_design in ['wrx', 'CPwrx']:
        pass
    elif sys_design in ['CPW', 'WOLA']:
        pass


def interf_power(sys_design:str):
    """Method to calculate the interference power for RC windows.
    
    Parameters
    ----------
    sys_design : str
        Name of the system.
        
    Returns
    -------
    interf_power : float
        Estimated interference power.
    """

    if sys_design in ['wtx', 'CPwtx']:
        pass
    elif sys_design in ['wrx', 'CPwrx']:
        pass
    elif sys_design in ['CPW', 'WOLA']:
        pass


# EoF
