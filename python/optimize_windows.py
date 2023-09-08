"""optimize_windows.py

Script for window optmization according to specifications.

Jul 22, 2023
"""

import numpy as np
import matplotlib.pyplot as plt
from optimization_tools import OptimizerTx
from utils import gen_channel_tensor_sym


if __name__ == '__main__':
    dft_len = 92
    cp_len = 10
    cs_len = 8
    tail_len = 8

    kwargs = {'system_design': 'wtx-OFDM', 'cp_len': cp_len, 'cs_len': cs_len,
              'dft_len': dft_len, 'tail_len': tail_len, 'rm_len': 20,
              'shift_len': 0}
    
    wtx = OptimizerTx(**kwargs)
    chan = np.random.randn(21) + 1j*np.random.randn(21)
    test = wtx.calculate_chann_matrices(chan)
    wtx.optimize_window(test)