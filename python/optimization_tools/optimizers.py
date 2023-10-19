"""optimizers.py

Script with class for the window optimization process.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


import os
import numpy as np
from numba import jit
from scipy.linalg import issymmetric
from ofdm_utils import (
    gen_idft_matrix, gen_add_redundancy_matrix, gen_rm_redundancy_matrix,
    gen_overlap_and_add_matrix, gen_dft_matrix, gen_circ_shift_matrix,
    gen_channel_tensor
)
from .utils import (reduce_variable_tx, reduce_variable_rx, gen_constraints_tx,
                    gen_constraints_rx, gen_costriants_q)
from .quadratic_programming import quadratic_solver


def optimization_fun(data:tuple):
    """Method to perform the optimization process.
    
    Parameters
    ----------
    data : tuple
        Tuple containing (system_design, dft_len, cp_len), where
            system_desing : str
                The w-OFDM system design.
            dft_len : int
                Number of bins in the DFT.
            cp_len : int
                Length of the cyclic prefix.
    """
    
    system_design, dft_len, cp_len, channel_path, window_path = data
    channel = np.load(channel_path)
    channel_avg = np.mean(channel, 1)
    if system_design in ['wtx', 'CPwtx']:
        tail_tx = 8
        opt_model = OptimizerTx(system_design, dft_len, cp_len, tail_tx)
        var_len = tail_tx+1
    elif system_design in ['wrx', 'CPwrx']:
        tail_rx = 10
        opt_model = OptimizerRx(system_design, dft_len, cp_len, tail_rx)
        var_len = int(tail_rx/2 + 1)
    elif system_design in ['CPW', 'WOLA']:
        tail_tx = 8
        tail_rx = 10
        opt_model = OptimizerTxRx(system_design, dft_len, cp_len, tail_tx, tail_rx)
        var_len = int((1+tail_rx/2)*(tail_tx+1))
    
    chann_ten = opt_model.calculate_chann_matrices(channel_avg)
    reg = 0
    H_mat = opt_model.gen_hessian(chann_ten) + reg*np.eye(var_len)
    if not issymmetric(H_mat):
        print('Not symmetric! Applying correction H = (H + H.T)/2.')
        H_mat = .5*(H_mat + H_mat.T)
    import pdb; pdb.set_trace()
    x, _ = opt_model.optimize(H_mat)
    
    os.makedirs(os.path.join(window_path, 'condition_number'), exist_ok=True)
    window_file_name = os.path.join(window_path, f'{system_design}_{cp_len}.npy')
    condition_number_file_name = os.path.join(
        window_path, 'condition_number', f'{system_design}_{cp_len}_{reg}.npy'
    )
    np.save(window_file_name, x)
    np.save(condition_number_file_name, np.linalg.cond(H_mat))

    return 'ok'


class OptimizerTx:
    """
    A class for the optmization process of wtx- and CPwtx-OFDM systems.

    Attributes
    ----------
    name : str
        The name of the system in question 'wtx' or 'CPwtx'.
    cp_len : int
        The number of samples in the cyclic prefix.
    cs_len : int
        The number of samples in the cyclic suffix.
    tail_len : int
        The number of samples in the window's tails.
    rm_len : int
        The number of samples removed at the receiver.
    shift_len : int
        The number of samples in the circular shift at the receiver.
    idft_mat : np.ndarray
        The matrix resposible for the IDFT.
    add_red_mat : np.ndarray
        The matrix responsible for adding redundancy to the system.
    rm_red_mat : np.ndarray
        The matrix responsible for removing redundancy from the system.
    circ_shift_mat : np.ndarray
        The matrix responsible for the circular shift.
    dft_mat : np.ndarray
        The matrix responsible for the 
    """

    @staticmethod
    @jit(nopython=True)
    def __gen_Q1(B_mat, C_mat):
        """Method to generate Q1 using numba (which is way faster)."""

        dft_len, n_samples = B_mat.shape
        Q1_mat = np.zeros((n_samples, n_samples), dtype=np.complex128)
        for i in range(n_samples):
            for j in range(i+1):
                for m in range(dft_len):
                    logic_aux = np.ones((dft_len,), dtype=np.bool_)
                    logic_aux[m] = 0
                    Q1_mat[i, j] += \
                        np.transpose(np.conj(B_mat[logic_aux, i])) \
                        * np.conj(C_mat[i, m]) @ B_mat[logic_aux, j] \
                        * C_mat[j, m]
                    
        return np.real(Q1_mat + np.transpose(Q1_mat) - np.diag(np.diag(Q1_mat)))

    def __init__(self, system_design, dft_len, cp_len, tail_len):
        """Constructor method.
        
        Parameters
        ----------
        system_desing : str
            w-OFDM system design.
        dft_len : int
            Length of DFT.
        cp_len : int
            Number of samples in cyclic prefix (CP).
        tail_len : int
            Number of samples in window tail.
        """

        self.name = system_design
        self.dft_len = dft_len
        self.cp_len = cp_len
        self.tail_len = tail_len
        if self.name == 'wtx':
            self.cs_len = self.tail_len
            self.rm_len = self.cp_len
            self.shift_len = 0
        elif self.name == 'CPwtx':
            self.cs_len = 0
            self.rm_len = self.cp_len - self.tail_len
            self.shift_len = self.tail_len
        # Set matrices:
        # Transmitter:
        self.idft_mat = gen_idft_matrix(self.dft_len)
        self.add_red_mat = gen_add_redundancy_matrix(self.dft_len, self.cp_len,
                                                     self.cs_len)
        # Receiver:
        self.rm_red_mat = gen_rm_redundancy_matrix(self.dft_len, 0, self.rm_len)
        self.overlap_add_mat = gen_overlap_and_add_matrix(self.dft_len, 0)
        self.circ_shift_mat = gen_circ_shift_matrix(self.dft_len, self.shift_len)
        self.dft_mat = gen_dft_matrix(self.dft_len)
        # Variable manipulation:
        self.reduce_var_mat = reduce_variable_tx(self.dft_len, self.cp_len,
                                                 self.cs_len, self.tail_len)

    def calculate_chann_matrices(self, channel_ir:np.ndarray):
        """Method to calculate channel matrices for the system.
        
        Parameters
        ----------
        channel_ir : np.ndarray
            Channel impulse response.
        
        Returns
        -------
        channel_tensor : np.ndarray
            Tensor with channel matrix, where the depth is the 3rd dim.
        """
        channel_tensor = gen_channel_tensor(channel_ir, self.dft_len,
                                            self.tail_len, 0, self.rm_len,
                                            self.cp_len, self.cs_len)
        
        return channel_tensor

    def gen_hessian(self, channel_tensor:np.ndarray):
        """Method to calculate the matrices necessary for the Hessian.
        
        Parameters
        ----------
        channel_tensor : np.ndarray
        
        Returns
        -------
        H_mat : np.ndarray
            Hessian matrix for the optimization.
        """

        B_mat = self.dft_mat@self.circ_shift_mat@self.overlap_add_mat \
            @ self.rm_red_mat@channel_tensor[:, :, 0]
        C_mat = self.add_red_mat@self.idft_mat
        D_mat = self.add_red_mat@self.circ_shift_mat@self.overlap_add_mat \
            @ self.rm_red_mat@np.sum(channel_tensor[:, :, 1:], 2)
        Q1_mat = self.__gen_Q1(B_mat, C_mat)
        Q2_mat = np.real(np.diag(np.diag(D_mat @ np.transpose(D_mat))))

        return np.transpose(self.reduce_var_mat) @ (2*(Q1_mat + Q2_mat)) \
            @ self.reduce_var_mat

    def optimize(self, H_mat):
        """
        Method to run optimization using our optimization algorithms.
        """

        A, b, C, d = gen_constraints_tx(self.tail_len)

        p = np.zeros((1+self.tail_len, 1), dtype=np.float64)
        x, n_iter = quadratic_solver(H_mat, p, A, b, C, d, epsilon=1e-9)

        return x, n_iter


class OptimizerRx:
    """
    A class for the optimization of the wrx- and CPwrx-OFDM systems.
    
    """

    @staticmethod
    @jit(nopython=True)
    def __gen_Q1(B_mat, C_mat):
        """Method to generate Q1 using numba (which is way faster)."""

        dft_len, n_samples = B_mat.shape
        Q1_mat = np.zeros((n_samples, n_samples), dtype=np.complex128)
        for i in range(n_samples):
            for j in range(i+1):
                for m in range(dft_len):
                    logic_aux = np.ones((dft_len,), dtype=np.bool_)
                    logic_aux[m] = 0
                    Q1_mat[i, j] += \
                        np.transpose(np.conj(B_mat[logic_aux, i])) \
                        * np.conj(C_mat[i, m]) @ B_mat[logic_aux, j] \
                        * C_mat[j, m]
                    
        return np.real(Q1_mat + np.transpose(Q1_mat) - np.diag(np.diag(Q1_mat)))

    def __init__(self, system_design:str, dft_len:int, cp_len:int,
                 tail_len:int):
        """
        Parameters
        ----------
        system_design : str
            w-OFDM system design.
        dft_len : int
            Length of DFT.
        cp_len : int
            Number of samples in cyclic prefix (CP).
        tail_len : int
            Number of samples in window tail.
        """

        self.name = system_design
        self.dft_len = dft_len
        self.cp_len = cp_len
        self.tail_len = tail_len
        if self.name == 'wrx':
            self.cs_len = int(self.tail_len/2)
            self.rm_len = int(self.cp_len - self.tail_len/2)
            self.shift_len = 0
        elif self.name == 'CPwrx':
            self.cs_len = 0
            self.rm_len = self.cp_len - self.tail_len
            self.shift_len = int(self.tail_len/2)
        # Set matrices:
        # Transmitter:
        self.idft_mat = gen_idft_matrix(self.dft_len)
        self.add_red_mat = gen_add_redundancy_matrix(self.dft_len, self.cp_len,
                                                     self.cs_len)
        # Receiver:
        self.rm_red_mat = gen_rm_redundancy_matrix(self.dft_len, self.tail_len,
                                                   self.rm_len)
        self.overlap_add_mat = gen_overlap_and_add_matrix(self.dft_len,
                                                          self.tail_len)
        self.circ_shift_mat = gen_circ_shift_matrix(self.dft_len, self.shift_len)
        self.dft_mat = gen_dft_matrix(self.dft_len)
        # Variable manipulation:
        self.reduce_var_mat = reduce_variable_rx(self.dft_len, self.tail_len)
        
    def calculate_chann_matrices(self, channel_ir:np.ndarray):
        """Method to calculate channel matrices for the system.
        
        Parameters
        ----------
        channel_ir : np.ndarray
            Channel impulse response.
        
        Returns
        -------
        channel_tensor : np.ndarray
            Tensor with channel matrix, where the depth is the 3rd dim.
        """
        channel_tensor = gen_channel_tensor(channel_ir, self.dft_len, 0,
                                            self.tail_len, self.rm_len,
                                            self.cp_len, self.cs_len)
        
        return channel_tensor
        
    def gen_hessian(self, channel_tensor:np.ndarray):
        """Method to calculate the matrices necessary for the Hessian.
        
        Parameters
        ----------
        channel_tensor : np.ndarray
        
        Returns
        -------
        H_mat : np.ndarray
            Hessian matrix for the optimization.
        """

        B_mat = self.dft_mat@self.circ_shift_mat@self.overlap_add_mat
        C_mat = self.rm_red_mat@channel_tensor[:, :, 0]@self.add_red_mat \
            @ self.idft_mat
        D_mat = self.rm_red_mat@np.sum(channel_tensor[:, :, 1:], 2) \
            @ self.add_red_mat@self.circ_shift_mat@self.overlap_add_mat
        Q1_mat = self.__gen_Q1(B_mat, C_mat)
        Q2_mat = np.real(np.diag(np.diag(D_mat @ np.transpose(D_mat))))

        return np.transpose(self.reduce_var_mat) @ (2*(Q1_mat + Q2_mat)) \
            @ self.reduce_var_mat
        
    def optimize(self, H_mat):
        """
        Method to run optimization using our algorithms.
        """

        A, b, C, d = gen_constraints_rx(self.tail_len)
        p = np.zeros((int(1 + self.tail_len/2), 1), dtype=np.float64)
        x, n_iter = quadratic_solver(H_mat, p, A, b, C, d, epsilon=1e-9)

        return x, n_iter


class OptimizerTxRx:
    """"""

    @staticmethod
    @jit(nopython=True)
    def __gen_Q1(B_mat:np.ndarray, C_mat:np.ndarray, D_mat:np.ndarray,
                 red_tx:np.ndarray, red_rx:np.ndarray):
        """
        Method to generate Q1 using numba (which is way faster).
        
        Parameters
        ----------
        B_mat : np.ndarray
            Array with B matrix.
        C_mat : np.ndarray
            Array with C matrix.
        D_mat : np.ndarray
            Array with D matrix.
        red_tx : np.ndarray
            Array with reduction matrix for the transmitter.
        red_rx : np.ndarray
            Array with reduction matrix for the receiver.
        """

        dft_len, _ = B_mat.shape
        len_rx, len_tx = C_mat.shape
        _, len_red_rx = red_rx.shape
        _, len_red_tx = red_tx.shape
        n_samples = int(len_red_tx*len_red_rx)
        mm = np.zeros((n_samples, int(dft_len*dft_len)), dtype=np.complex128)
        for i in range(dft_len):
            for j in range(dft_len):
                if i == j:
                    continue
                M_mat = np.zeros((len_rx, len_tx), dtype=np.complex128)
                for m in range(len_rx):
                    for n in range(len_tx):
                        M_mat[m, n] = B_mat[i, m] * D_mat[n, j]
                M_mat *= C_mat
                M_prim = red_rx.T @ M_mat @ red_tx
                m = M_prim.flatten()
                mm[:, int(dft_len*i + j)] = m
        Q1_mat = mm @ np.conj(mm.T)

        return Q1_mat

    @staticmethod
    @jit(nopython=True)
    def __gen_Q2(B_mat:np.ndarray, C_mat:np.ndarray, D_mat:np.ndarray,
                 red_tx:np.ndarray, red_rx:np.ndarray):
        """
        Method to generate Q1 using numba (which is way faster).
        
        Parameters
        ----------
        B_mat : np.ndarray
            Array with B matrix.
        C_mat : np.ndarray
            Array with C matrix.
        D_mat : np.ndarray
            Array with D matrix.
        red_tx : np.ndarray
            Array with reduction matrix for the transmitter.
        red_rx : np.ndarray
            Array with reduction matrix for the receiver.
        """

        dft_len, _ = B_mat.shape
        len_rx, len_tx = C_mat.shape
        _, len_red_rx = red_rx.shape
        _, len_red_tx = red_tx.shape
        n_samples = int(len_red_tx*len_red_rx)
        mm = np.zeros((n_samples, int(dft_len*dft_len)), dtype=np.complex128)
        for i in range(dft_len):
            for j in range(dft_len):
                M_mat = np.zeros((len_rx, len_tx), dtype=np.complex128)
                for m in range(len_rx):
                    for n in range(len_tx):
                        M_mat[m, n] = B_mat[i, m] * D_mat[n, j]
                M_mat *= C_mat
                M_prim = red_rx.T @ M_mat @ red_tx
                m = M_prim.flatten()
                mm[:, int(dft_len*i + j)] = m
        Q2_mat = mm @ np.conj(mm.T)

        return Q2_mat

    def __init__(self, system_design:str, dft_len:int, cp_len:int,
                 tail_tx_len:int, tail_rx_len:int):
        """
        Parameters
        ----------
        system_design : str
            w-OFDM system design.
        dft_len : int
            Length of DFT.
        cp_len : int
            Number of samples in cyclic prefix (CP).
        tail_tx_len : int
            Number of samples in transmitter window tail.
        tail_rx_len : int
            Number of samples in receiver window tail.
        """

        self.name = system_design
        self.dft_len = dft_len
        self.cp_len = cp_len
        self.tail_tx_len = tail_tx_len
        self.tail_rx_len = tail_rx_len
        if self.name == 'WOLA':
            self.cs_len = int(self.tail_tx_len + self.tail_rx_len/2)
            self.rm_len = int(self.cp_len - self.tail_rx_len/2)
            self.shift_len = 0
        elif self.name == 'CPW':
            self.cs_len = self.tail_tx_len
            self.rm_len = self.cp_len - self.tail_rx_len
            self.shift_len = int(self.tail_rx_len/2)
        # Set matrices:
        # Transmitter:
        self.idft_mat = gen_idft_matrix(self.dft_len)
        self.add_red_mat = gen_add_redundancy_matrix(self.dft_len, self.cp_len,
                                                     self.cs_len)
        # Receiver:
        self.rm_red_mat = gen_rm_redundancy_matrix(
            self.dft_len, self.tail_rx_len, self.rm_len
        )
        self.overlap_add_mat = gen_overlap_and_add_matrix(self.dft_len,
                                                          self.tail_rx_len)
        self.circ_shift_mat = gen_circ_shift_matrix(self.dft_len,
                                                    self.shift_len)
        self.dft_mat = gen_dft_matrix(self.dft_len)
        # Sparsity exploration:
        self.reduce_var_tx_mat = reduce_variable_tx(
            self.dft_len, self.cp_len, self.cs_len, self.tail_tx_len
        )
        self.reduce_var_rx_mat = reduce_variable_rx(self.dft_len,
                                                    self.tail_rx_len)

    def calculate_chann_matrices(self, channel_ir:np.ndarray):
        """Method to calculate channel matrices for the system.
        
        Parameters
        ----------
        channel_ir : np.ndarray
            Channel impulse response.
        
        Returns
        -------
        channel_tensor : np.ndarray
            Tensor with channel matrix, where the depth is the 3rd dim.
        """

        channel_tensor = gen_channel_tensor(
            channel_ir, self.dft_len, self.tail_tx_len, self.tail_rx_len,
            self.rm_len, self.cp_len, self.cs_len
        )

        return channel_tensor
    
    def gen_hessian(self, channel_tensor:np.ndarray):
        """Method to calculate the matrices necessary for the Hessian.
        
        Parameters
        ----------
        channel_tensor : np.ndarray
        
        Returns
        -------
        H_mat : np.ndarray
            Hessian matrix for the optimization.
        """

        B_mat = self.dft_mat@self.circ_shift_mat@self.overlap_add_mat
        C_mat = self.rm_red_mat@channel_tensor[:, :, 0]
        C2_mat = self.rm_red_mat@np.sum(channel_tensor[:, :, 1:], 2)
        D_mat = self.add_red_mat@self.idft_mat
        Q1_mat = self.__gen_Q1(B_mat, C_mat, D_mat,
                               self.reduce_var_tx_mat.astype(np.complex128),
                               self.reduce_var_rx_mat.astype(np.complex128))

        Q2_mat = self.__gen_Q2(B_mat, C2_mat, D_mat,
                               self.reduce_var_tx_mat.astype(np.complex128),
                               self.reduce_var_rx_mat.astype(np.complex128))

        return 2*(Q1_mat + Q2_mat)


# EoF
