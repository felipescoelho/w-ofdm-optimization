"""wofdm_simulation.py

Script with w-OFDM simulation methods.

luizfelipe.coelho@smt.ufrj.br
Oct 10, 2023
"""


import os
import numpy as np
from numba import jit
from .transmitter import (gen_add_redundancy_matrix, gen_idft_matrix,
                          gen_rc_window_tx)
from .receiver import (gen_circ_shift_matrix, gen_dft_matrix,
                       gen_overlap_and_add_matrix, gen_rc_window_rx,
                       gen_rm_redundancy_matrix)
from optimization_tools import reduce_variable_tx, reduce_variable_rx


def simulation_fun(data:tuple):
    """Method to perform the simulation process.
    
    Parameters
    ----------
    data : tuple
        Tuple containing (system_design, dft_len, cp_len, tail_tx,
        tail_rx), where
            system_design : str
                The w-OFDM system design.
            dft_len : int
                Number of bins in the DFT.
            cp_len : int
                Length of the cyclic prefix.
            tail_tx : int
                Length of the tail in the transmitter window.
            tail_rx : int
                Length of the tail in the receiver window.
            window_path : str
                Path to load windows.
    """

    system_design, dft_len, cp_len, tail_tx, tail_rx = data[0:5]
    channel_path, window_path, ensemble, snr_arr, no_symbols = data[5:]

    window_file_path = os.path.join(
        window_path, f'{system_design}_{cp_len}.npy'
    )
    channel_models = np.load(channel_path)
    if system_design in ['WOLA', 'CPW']:
        win_tail_tx, win_tail_rx = np.load(window_file_path)
        cs_len = int(tail_tx + tail_rx/2) if system_design == 'WOLA' else tail_tx
    elif system_design in ['wtx', 'CPwtx']:
        win_tail_tx = np.load(window_file_path)
        win_tail_rx = np.array([1], ndmin=2)
        cs_len = tail_tx if system_design == 'wtx' else 0
    elif system_design in ['wrx', 'CPwrx']:
        win_tail_rx = np.load(window_file_path)
        win_tail_tx = np.array([1], ndmin=1)
        cs_len = int(tail_rx/2) if system_design == 'wrx' else 0
    win_tx = np.diagflat(reduce_variable_tx(dft_len, cp_len, cs_len,
                                            tail_tx)@win_tail_tx)
    win_rx = np.diagflat(reduce_variable_rx(dft_len, tail_rx)@win_tail_rx)
    sim_model = wOFDMSystem(system_design, dft_len, cp_len, tail_tx, tail_rx)
    sim_model.run_simulation(channel_models, win_tx, win_rx, ensemble, snr_arr,
                             no_symbols)


class wOFDMSystem:
    """
    A class to simulate each system according to certain specifications.
    
    Attributes
    ----------

    """

    @staticmethod
    @jit(nopython=True)
    def __run_sim_mc(tx_mat:np.ndarray, rx_mat:np.ndarray, no_symbols:int,
                     channel_models:np.ndarray, ensemble:int, snr_arr:float,
                     tail_tx_len:int, ser_flag=True):
        """Method to run simulation using numba.

        In this simulation, we transmit a certain number of w-OFDM
        symbols, `no_symbols`, where the first symbol is our pilot.
        Furthermore, we estimate the symbol error rate (SER) for each
        signal to noise ratio (SNR) in `snr_arr` using a Monte Carlo
        process for each channel model in `channel_models`, and
        averaging its results.
        
        Parameters
        ----------
        tx_mat : np.ndarray
            Array representing transmitter.
        rx_mat : np.ndarray
            Array representing receiver.
        no_symbols : int
            Number of symbols in a single frame.
        channel_models : np.ndarray
            Array with channel models.
        ensemble : int
            Number of repetitions in Monte Carlo process.
        snr_arr : np.ndarray
            Signal to noise ratio in dB.
        ser_flag : bool (default=True)
            Symbol error rate flag, if true we evaluate the SER and not
            the BER (bit error rate).
        """

        def awgn(x, snr) -> np.ndarray:
            """Method to add white Gaussian noise.
            
            Parameters
            ----------
            x : np.ndarray
                Signal to be poluted with noise.
            snr : float
                Signal to noise ratio in dB.
            
            Returns
            -------
            y : np.ndarray
                Noisy signal.
            """

            Px = np.vdot(x, x)
            n = np.random.randn(len(x),) + 1j*np.random.randn(len(x),)
            Pn = np.vdot(n, n)
            n_adjusted = np.sqrt(Px*10**(-.1*snr)/Pn)*n

            return x + n_adjusted
        
        def decision(x:np.ndarray, symbols:np.ndarray):
            """Method to decide which symbol is received.
            
            Criterion: Shortest distance.
            
            Parameters
            ----------
            x : np.ndarray
                Received signal.
            symbols : np.ndarray
                Possible values for reception.
            
            Returns
            -------
            y : np.ndarray

            """
            n_rows, n_cols = x.shape
            y = np.zeros((n_rows, n_cols), dtype=np.complex128)
            for i in range(n_rows):
                for j in range(n_cols):
                    min_idx = np.argmin(np.abs(symbols - x[i, j]))
                    y[i, j] = symbols[min_idx]

            return y
        
        n_elements = len(snr_arr)
        ser = np.zeros((n_elements,), dtype=np.float64)
        for idx, snr in enumerate(snr_arr):
            result = 0
            for ch_idx in range(channel_models.shape[1]):
                result_mc = 0
                for _ in range(ensemble):
                    if ser_flag:
                        symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j,
                                            -1-1j, -1+1j, -1+3j, 1-3j, 1-1j,
                                            1+1j, 1+3j, 3-3j, 3-1j, 3+1j,
                                            3+3j))  # 16-QAM
                        signal_digmod = np.random.choice(
                            symbols, size=(tx_mat.shape[1], no_symbols),
                            replace=True
                        )
                        frame_tx = (tx_mat @ signal_digmod).T
                        # overlap-and-add -> serialize
                        signal_ov = frame_tx[:, tail_tx_len:]
                        signal_ov[:-1, -tail_tx_len:] = \
                            frame_tx[1:, :tail_tx_len] \
                            + frame_tx[:-1, -tail_tx_len:]
                        signal_tx = np.hstack((frame_tx[0, :tail_tx_len],
                                               signal_ov.flatten()))
                        # Convolution 'valid' + noise
                        channel_model = channel_models[:, ch_idx]
                        signal_conv = np.convolve(channel_model, signal_tx)
                        signal_rx = awgn(
                            signal_conv[:-(len(channel_model)+tail_tx_len-1)], snr
                        )
                        # parallelize -> preprocess -> recover
                        frame_rx = signal_rx.reshape((no_symbols,
                                                      rx_mat.shape[1]))
                        signal_preproc = rx_mat @ frame_rx.T
                        channel_est = signal_preproc[:, 0]/signal_digmod[:, 0]
                        channel_est_mat = np.repeat(
                            channel_est, no_symbols-1
                        ).reshape((tx_mat.shape[1], no_symbols-1))
                        recovered_signal = signal_preproc[:, 1:]/channel_est_mat
                        symbols_est = decision(recovered_signal, symbols)
                        result_mc += np.mean((symbols_est != signal_digmod[:, 1:]))
                result += result_mc/ensemble
            ser[idx] = result/channel_models.shape[1]

        return ser
        
    def __init__(self, system_design:str, dft_len:int, cp_len:int, tail_tx:int,
                 tail_rx:int):
        """
        Parameters
        ----------
        system_design : str
            w-OFDM system design.
        dft_len : int
            Length of DFT.
        cp_len : int
            Length of cyclic prefix
        tail_tx : int
            Number of samples in transmitter window tail.
        tail_rx : int
            Number of samples in receiver window tail.
        """

        self.name = system_design
        self.dft_len = dft_len
        self.cp_len = cp_len
        self.tail_tx = tail_tx
        self.tail_rx = tail_rx
        if self.name == 'wtx':
            self.cs_len = self.tail_tx
            self.rm_len = self.cp_len
            self.shift_len = 0
        elif self.name == 'CPwtx':
            self.cs_len = 0
            self.rm_len = self.cp_len - self.tail_tx
            self.shift_len = self.tail_tx
        elif self.name == 'wrx':
            self.cs_len = int(self.tail_rx/2)
            self.rm_len = int(self.cp_len - self.tail_rx/2)
            self.shift_len = 0
        elif self.name == 'CPwrx':
            self.cs_len = 0
            self.rm_len = self.cp_len - self.tail_rx
            self.shift_len = int(self.tail_rx/2)
        elif self.name == 'CPW':
            self.cs_len = self.tail_tx
            self.rm_len = self.cp_len - self.tail_rx
            self.shift_len = int(self.tail_rx/2)
        elif self.name == 'WOLA':
            self.cs_len = int(self.tail_tx + self.tail_rx/2)
            self.rm_len = int(self.cp_len - self.tail_rx/2)
            self.shift_len = 0
        # Set matrices:
        # Transmitter:
        self.idft_mat = gen_idft_matrix(self.dft_len)
        self.add_red_mat = gen_add_redundancy_matrix(self.dft_len, self.cp_len,
                                                     self.cs_len)
        # Receiver:
        self.rm_red_mat = gen_rm_redundancy_matrix(self.dft_len, self.tail_rx,
                                                   self.rm_len)
        self.overlap_add_mat = gen_overlap_and_add_matrix(self.dft_len,
                                                          self.tail_rx)
        self.circ_shift_mat = gen_circ_shift_matrix(self.dft_len, self.shift_len)
        self.dft_mat = gen_dft_matrix(self.dft_len)

    def run_simulation(self, channel_models:np.ndarray, window_tx:np.ndarray,
                       window_rx:np.ndarray, ensemble:int, snr_arr:np.ndarray,
                       no_symbols:int):
        """
        Method to run simulation of `channel_models` over ensemble in
        Monte Carlo.

        Parameters
        ----------
        channel_models : np.ndarray
            Array containing TDL channel models.
        window_tx : np.ndarray
            Array with transmitter window.
        window_rx : np.ndarray
            Array with receiver window. 
        ensemble : int
            Number of repetitions in Monte Carlo process.
        snr_arr : np.ndarray
            Signal to noise ratio in dB.
        no_symbols : int
            Number of symbols per transmitted frame.
        """

        tx_mat = window_tx @ self.add_red_mat @ self.idft_mat
        # tx_mat = self.add_red_mat@self.idft_mat
        rx_mat = self.dft_mat @ self.circ_shift_mat @ self.overlap_add_mat \
            @ window_rx @ self.rm_red_mat
        ser = self.__run_sim_mc(tx_mat, rx_mat, no_symbols, channel_models,
                                ensemble, snr_arr, self.tail_tx, ser_flag=True)
        print(ser)

# EoF
