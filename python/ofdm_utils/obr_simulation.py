"""obr_simulation.py

Script with simulation methods to estimate the out of band radiation
(OBR).

luizfelipe.coelho@smt.ufrj.br
Nov 28, 2023
"""


import os
import numpy as np
from .transmitter import (gen_idft_matrix, gen_add_redundancy_matrix,
                          gen_rc_window_tx)


def obr_fun(data:tuple):
    """Method to perform the OBR estimation.
    
    Parameters
    ----------
    data : tuple
        Tuple containing (system_design, dft_len, cp_len, tail_tx,
        tail_rx, window_path, folder_path), where
            system_design : str
                Name of the w-OFDM system.
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
            folder_path : str
                Path where we save the results.
    """
    from optimization_tools import reduce_variable_tx

    system_design, dft_len, cp_len, tail_tx, tail_rx = data[0:5]
    window_path, folder_path = data[5:]

    window_file_path = os.path.join(
        window_path, f'{system_design}_{cp_len}.npy'
    )
    if system_design in ['WOLA', 'CPW']:
        win_tail = np.load(window_file_path)
        win_tail_tx = win_tail[:tail_tx+1].reshape(tail_tx+1, 1)
        cs_len = int(tail_tx + tail_rx/2) if system_design == 'WOLA' else tail_tx
    elif system_design in ['wtx', 'CPwtx']:
        win_tail_tx = np.load(window_file_path)
        cs_len = tail_tx if system_design == 'wtx' else 0
    elif system_design in ['wrx', 'CPwrx']:
        win_tail_tx = np.array([1], ndmin=1)
        cs_len = int(tail_rx/2) if system_design == 'wrx' else 0
    elif system_design == 'CP':
        win_tail_tx = np.array([[1.]], ndmin=2, dtype=np.float64)
        cs_len = 0
    win_tx = np.diagflat(reduce_variable_tx(dft_len, cp_len, cs_len,
                                            tail_tx)@win_tail_tx)
    obr_model = wOFDMSystem(system_design, dft_len, cp_len, tail_tx, tail_rx,
                            folder_path)
    obr_model.estimate_obr(win_tx, 200*1e-9)


class wOFDMSystem:
    """"""

    @staticmethod
    def overlap_and_add(x:np.ndarray, beta:int):
        """"""
        x_aux = x[:, beta:]
        x_aux[:-1, -beta:] += x[1:, :beta]
        y = np.hstack((x[0, :beta], x_aux.flatten()))

        return y
    
    @staticmethod
    def __psd_estimate(x:np.ndarray, fft_len:int):
        """Method to estimate the PSD.
        
        Parameters
        ----------
        x : np.ndarray
            Input signal.
        fft_len : int
            Length of the FFT for the analysis.
        """

        x_len = len(x)
        no_slices = int(np.floor(x_len/fft_len))
        X = np.zeros((fft_len,), dtype=np.float64)
        for no_slice in range(no_slices):
            x_sliced = x[int(no_slice*fft_len):int((no_slice+1)*fft_len)]
            X += np.abs(np.fft.fftshift(np.fft.fft(x_sliced, fft_len)))**2

        return X/no_slices

    def __init__(self, system_design:str, dft_len:int, cp_len:int, tail_tx:int,
                 tail_rx:int, folder_path:str):
        """"""

        self.name = system_design
        self.dft_len = dft_len
        self.cp_len = cp_len
        self.tail_tx = tail_tx
        self.tail_rx = tail_rx
        self.folder_path = folder_path
        if self.name == 'wtx':
            self.cs_len = self.tail_tx
        elif self.name == 'CPwtx':
            self.cs_len = 0
        elif self.name == 'wrx':
            self.cs_len = int(self.tail_rx/2)
        elif self.name == 'CPwrx':
            self.cs_len = 0
        elif self.name == 'CPW':
            self.cs_len = self.tail_tx
        elif self.name == 'WOLA':
            self.cs_len = int(self.tail_tx + self.tail_rx/2)
        elif self.name == 'CP':
            self.cs_len = 0
        self.idft_mat = gen_idft_matrix(self.dft_len)
        self.add_red_mat = gen_add_redundancy_matrix(self.dft_len, self.cp_len,
                                                     self.cs_len)

    def analytical_psd(self, sampling_period, win_tx_mat, guard_band):
        """"""

        from scipy.signal import firwin
        # Definitions:
        f_axis = np.linspace(-self.dft_len/2, self.dft_len/2-1, 2*self.dft_len) \
            / (sampling_period*self.dft_len)
        delta_f = 1/(self.dft_len*sampling_period)
        win_tx = np.diag(win_tx_mat)
        sigma_x_sqrd = ((self.dft_len-guard_band) / self.dft_len)**2
        # Interpolation Filter:
        g_i = firwin(21, f_axis[-guard_band], window='boxcar', fs=1/sampling_period)
        G_i = np.abs(np.fft.fftshift(np.fft.fft(g_i, 2*self.dft_len)))**2
        print(G_i)
        # Pre-calculation:
        scale_coeff = G_i*(self.dft_len*sigma_x_sqrd
                       / ((self.dft_len+self.cp_len+self.cs_len)*sampling_period))
        cp_corr_win = np.sum(
            win_tx[:self.cp_len]*win_tx[self.dft_len:self.dft_len+self.cp_len]
        )
        cs_corr_win = np.sum(
            win_tx[self.dft_len+self.cp_len:self.dft_len+self.cp_len+self.cs_len]
            * win_tx[self.cp_len:self.cp_len+self.cs_len]
        )
        corr_win = np.sum(win_tx**2)
        # Calculation
        S =  scale_coeff*(corr_win + 2*(cp_corr_win + cs_corr_win)
                          * np.cos(f_axis/delta_f))

        return S

    def estimate_obr(self, win_tx_mat, samp_period):
        """"""

        no_symbols = 128
        guard_band = 48
        fft_len = 2*self.dft_len

        subcar_alloc_mat = np.vstack((
            np.zeros((1, self.dft_len-int(2*(guard_band)+1))),
            np.hstack((np.eye(int(self.dft_len/2)-guard_band),
                       np.zeros((int(self.dft_len/2)-guard_band, int(self.dft_len/2)-(guard_band+1))))),
            np.zeros((int(2*guard_band), self.dft_len-int(2*guard_band + 1))),
            np.hstack((np.zeros((int(self.dft_len/2)-(guard_band+1), int(self.dft_len/2)-guard_band)),
                       np.eye(int(self.dft_len/2)-(guard_band+1))))
        ))
        symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j, -1+1j,
                            -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j, 3-1j, 3+1j,
                            3+3j))  # 16-QAM
        X = np.random.choice(symbols, size=(self.dft_len-int(2*(guard_band)+1),
                                            no_symbols), replace=True)
        x_opt = win_tx_mat@self.add_red_mat@self.idft_mat@subcar_alloc_mat@X
        x_rc = gen_rc_window_tx(
            self.dft_len, self.cp_len, self.cs_len, self.tail_tx
        )@self.add_red_mat@self.idft_mat@subcar_alloc_mat@X
        if self.name in ['wtx', 'CPwtx', 'WOLA', 'CPW']:
            x_opt = self.overlap_and_add(x_opt.T, self.tail_tx)
            x_rc = self.overlap_and_add(x_rc.T, self.tail_tx)
        else:
            x_opt = x_opt.T.flatten()
            x_rc = x_rc.T.flatten()
        x = (self.idft_mat@subcar_alloc_mat@X).T.flatten()
        x_cp = (self.add_red_mat@self.idft_mat@subcar_alloc_mat@X).T.flatten()

        X_est_opt = self.__psd_estimate(x_opt, fft_len)
        X_est_rc = self.__psd_estimate(x_rc, fft_len)
        X_est_cp = self.__psd_estimate(x_cp, fft_len)
        X_est = self.__psd_estimate(x, fft_len)

        interp_factor = fft_len/self.dft_len
        gb_corrected = int(interp_factor*guard_band)

        # OBR Calculation:
        obr_opt = np.mean(np.hstack((
            X_est_opt[:int(gb_corrected/2)],
            X_est_opt[-int(gb_corrected/2):]
        )))
        obr_rc = np.mean(np.hstack((
            X_est_rc[:int(gb_corrected/2)],
            X_est_rc[-int(gb_corrected/2):]
        )))
        obr_cp = np.mean(np.hstack((
            X_est_cp[:int(gb_corrected/2)],
            X_est_cp[-int(gb_corrected/2):]
        )))
        
        # Main freq. band Calculation:
        # mf_band_opt = np.hstack((
        #     X_est_opt_avg[
        #         int(guardband_wofdm/2):int((self.dft_len/2)*interp_factor_wofdm)
        #     ],
        #     X_est_opt_avg[:-int(guardband_wofdm/2)-int(interp_factor_wofdm)]
        # ))
        
        print(obr_opt)
        print(obr_cp)
        print(obr_rc)
        
        S_f = self.analytical_psd(samp_period, win_tx_mat, guard_band)
        # return X_est_avg, X_est_rc_avg, X_est_cp_avg

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(10*np.log10(S_f))
        # ax.plot(10*np.log10(X_est_cp))
        # ax.plot(10*np.log10(X_est_opt))
        # ax.plot(10*np.log10(X_est_rc))
        # ax.plot(10*np.log10(X_est), ':', c='k')

        plt.show()


# EoF
