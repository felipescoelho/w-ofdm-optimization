"""timefreq_analysis.py

Script with method to perfm the time-frequency analysis of the w-OFDM
systems.

luizfelipe.coelho@smt.ufrj.br
Nov 28, 2023
"""


import os
import numpy as np
from .transmitter import (gen_idft_matrix, gen_add_redundancy_matrix,
                          gen_rc_window_tx)


def timefreq_fun(data:tuple):
    """Method to perform the time-frequency analysis.
    
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
        cs_len = int(tail_tx + tail_rx/2) if system_design == 'CPW' else tail_tx
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
    opt_dict, rc_dict, cp_dict = obr_model.estimate_obr(win_tx, 200*1e-9)
    
    path_to_save = os.path.join(folder_path, 'timefreq')
    os.makedirs(path_to_save, exist_ok=True)
    path_opt = os.path.join(path_to_save, f'opt_{system_design}_{cp_len}.npz')
    path_rc = os.path.join(path_to_save, f'rc_{system_design}_{cp_len}.npz')
    path_cp = os.path.join(path_to_save, f'CP_{cp_len}.npz')
    
    np.savez(path_opt, **opt_dict)
    np.savez(path_rc, **rc_dict)
    np.savez(path_cp, **cp_dict)


class wOFDMSystem:
    """"""

    @staticmethod
    def overlap_and_add(x:np.ndarray, beta:int):
        """Method to perform the overlap and add operation.
        
        Parameters
        ----------
        x : np.ndarray
            OFDM frame (simbols x samples)
        beta : int
            Samples to overlap.
        
        Returns
        -------
        y : np.ndarray
            Vectorized OFDM signal.
        """
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
        x_sliced = x[int((no_slice+1)*fft_len):]
        X += np.abs(np.fft.fftshift(np.fft.fft(x_sliced, fft_len)))**2

        return X/(no_slices+1)

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
        elif self.name == 'WOLA':
            self.cs_len = self.tail_tx
        elif self.name == 'CPW':
            self.cs_len = int(self.tail_tx + self.tail_rx/2)
        elif self.name == 'CP':
            self.cs_len = 0
        self.idft_mat = gen_idft_matrix(self.dft_len)
        self.add_red_mat = gen_add_redundancy_matrix(self.dft_len, self.cp_len,
                                                     self.cs_len)

    def analytical_psd(self, sampling_period:float, win_tx_mat:np.ndarray,
                       guard_band:int, fft_len:int):
        """"""

        def calculate_correlations(win:np.ndarray):
            cp_corr = np.sum(
                win[:self.cp_len]*win[self.dft_len:self.dft_len+self.cp_len]
            )
            cs_corr = np.sum(
                win[self.dft_len+self.cp_len:self.dft_len+self.cp_len+self.cs_len]
                * win[self.cp_len:self.cp_len+self.cs_len]
            )
            whole_corr = np.sum(win**2)

            return whole_corr, cp_corr, cs_corr

        from scipy.signal import firwin
        # Definitions:
        upsampling_factor = fft_len/self.dft_len
        f_axis = np.linspace(-.5, .5-(1/fft_len), fft_len)/sampling_period
        delta_f = 1/(self.dft_len*sampling_period)
        win_tx = np.diag(win_tx_mat)
        sigma_x_sqrd = (self.dft_len/(self.dft_len-guard_band))**2
        # Interpolation Filter:
        filt_len = int(self.dft_len+self.cp_len+self.cs_len)
        g_i = firwin(filt_len, [f_axis[int(fft_len/2 + 1*upsampling_factor)],
                                f_axis[-int(guard_band*upsampling_factor)]],
                     window='boxcar', fs=1/sampling_period, pass_zero=False)
        # Pre-calculation:
        # Optimal:
        g_i_opt = g_i*win_tx
        G_i_opt = np.abs(np.fft.fftshift(np.fft.fft(g_i_opt, fft_len)))**2
        scale_coeff_opt = G_i_opt*(self.dft_len*sigma_x_sqrd
                       / ((self.dft_len+self.cp_len+self.cs_len)))/upsampling_factor
        corr_opt, cp_corr_opt, cs_corr_opt = calculate_correlations(win_tx)
        # Raised Cosine:
        win_rc = np.diag(gen_rc_window_tx(self.dft_len, self.cp_len,
                                             self.cs_len, self.tail_tx))
        g_i_rc = g_i*win_rc
        G_i_rc = np.abs(np.fft.fftshift(np.fft.fft(g_i_rc, fft_len)))**2
        scale_coeff_rc = G_i_rc*(self.dft_len*sigma_x_sqrd
                       / ((self.dft_len+self.cp_len+self.cs_len)))/upsampling_factor
        corr_rc, cp_corr_rc, cs_corr_rc = calculate_correlations(win_rc)
        # CP:
        G_i_cp = np.abs(np.fft.fftshift(np.fft.fft(g_i, fft_len)))**2
        scale_coeff = G_i_cp*(self.dft_len*sigma_x_sqrd
                       / ((self.dft_len+self.cp_len)))/upsampling_factor
        cp_corr_cp = self.cp_len
        corr_cp = self.dft_len + self.cp_len
        # Calculation
        S_opt = scale_coeff_opt*(
            corr_opt + 2*(cp_corr_opt+cs_corr_opt) * np.cos(f_axis/delta_f)
        )
        S_rc = scale_coeff_rc*(
            corr_rc + 2*(cp_corr_rc+cs_corr_rc) * np.cos(f_axis/delta_f)
        )
        S_cp = scale_coeff*(
            corr_cp + 2*(cp_corr_cp)*np.cos(f_axis/delta_f)
        )

        return S_opt, S_rc, S_cp

    def estimate_obr(self, win_tx_mat:np.ndarray, samp_period:float):
        """"""

        no_symbols = 256
        guard_band = 48
        fft_len = 8*self.dft_len

        subcar_alloc_mat = np.vstack((
            np.zeros((1, self.dft_len-int(2*guard_band))),
            np.hstack((np.eye(int(self.dft_len/2)-guard_band),
                       np.zeros((int(self.dft_len/2)-guard_band,
                                 int(self.dft_len/2)-guard_band)))),
            np.zeros((int(2*guard_band - 1), self.dft_len-int(2*guard_band))),
            np.hstack((np.zeros((int(self.dft_len/2)-guard_band,
                                 int(self.dft_len/2)-guard_band)),
                       np.eye(int(self.dft_len/2)-guard_band)))
        ))
        symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j, -1+1j,
                            -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j, 3-1j, 3+1j,
                            3+3j))  # 16-QAM
        # symbols /= np.sqrt(np.mean(np.abs(symbols)**2))
        X = np.random.choice(symbols, size=(self.dft_len-int(2*guard_band),
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
        # x = (self.idft_mat@subcar_alloc_mat@X).T.flatten()
        x_cp = (self.add_red_mat@self.idft_mat@subcar_alloc_mat@X).T.flatten()
        X_est_opt = self.__psd_estimate(x_opt, fft_len)
        X_est_rc = self.__psd_estimate(x_rc, fft_len)
        X_est_cp = self.__psd_estimate(x_cp, fft_len)
        # X_est = self.__psd_estimate(x, fft_len)

        f_axis = np.linspace(-.5, .5-(1/fft_len), fft_len) / (200*1e-9)

        interp_factor = fft_len/self.dft_len
        gb_corrected = int(interp_factor*guard_band)

        # OBR Calculation:
        obr_opt = np.mean(np.hstack((
            X_est_opt[:gb_corrected], X_est_opt[-gb_corrected:]
        )))
        obr_rc = np.mean(np.hstack((
            X_est_rc[:gb_corrected], X_est_rc[-gb_corrected:]
        )))
        obr_cp = np.mean(np.hstack((
            X_est_cp[:gb_corrected], X_est_cp[-gb_corrected:]
        )))
        
        # Main freq. band Calculation:
        mf_band_opt = np.hstack((
            X_est_opt[gb_corrected:int(fft_len/2)],
            X_est_opt[-int(fft_len/2 - interp_factor):-gb_corrected]
        ))
        mf_band_rc = np.hstack((
            X_est_rc[gb_corrected:int(fft_len/2)],
            X_est_rc[-int(fft_len/2 - interp_factor):-gb_corrected]
        ))
        mf_band_cp = np.hstack((
            X_est_cp[gb_corrected:int(fft_len/2)],
            X_est_cp[-int(fft_len/2 - interp_factor):-gb_corrected]
        ))

        S_opt, S_rc, S_cp = self.analytical_psd(samp_period, win_tx_mat,
                                                guard_band, fft_len)
        
        opt_out = {'X_est_opt': X_est_opt, 'S_opt': S_opt, 'f_axis': f_axis,
                   'obr_opt': obr_opt, 'mf_band_opt': mf_band_opt}
        rc_out = {'X_est_rc': X_est_rc, 'S_rc': S_rc, 'f_axis': f_axis,
                  'obr_rc': obr_rc, 'mf_band_rc': mf_band_rc}
        cp_out = {'X_est_cp': X_est_cp, 'S_cp': S_cp, 'f_axis': f_axis,
                  'obr_cp': obr_cp, 'mf_band_cp': mf_band_cp}

        return opt_out, rc_out, cp_out


# EoF
