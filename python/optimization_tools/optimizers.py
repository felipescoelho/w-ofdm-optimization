"""optimizers.py

Script with class for the window optimization process.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


import numpy as np
from numba import jit
from scipy.linalg import issymmetric
from ofdm_utils import (
    gen_idft_matrix, gen_add_redundancy_matrix, gen_rm_redundancy_matrix,
    gen_overlap_and_add_matrix, gen_dft_matrix, gen_circ_shift_matrix,
    gen_channel_tensor
)
from .utils import reduce_variable_tx, gen_constraints_tx
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
    
    system_design, dft_len, cp_len, channel_path = data
    channel = np.load(channel_path)
    channel_avg = np.mean(channel, 1)
    if system_design in ['wtx', 'CPwtx']:
        tail_tx = 8
        opt_model = OptimizerTx(system_design, dft_len, cp_len, tail_tx)
        chann_ten = opt_model.calculate_chann_matrices(channel_avg)
        reg = 1e-4
        # H_mat = (1-reg)*opt_model.gen_hessian(chann_ten) + reg*np.eye(9)
        H_mat = opt_model.gen_hessian(chann_ten) + reg*np.diagflat(np.linspace(
            1, 2, 9
        ))
        # print(np.linalg.eigvals(H_mat))
        # print(np.linalg.cond(H_mat))
        if not issymmetric(H_mat):
            print('Not symmetric!')
            H_mat = .5 * (H_mat + H_mat.T)
            # print(np.linalg.eigvals(H_mat))
            # print(np.linalg.cond(H_mat))
    elif system_design in ['wrx', 'CPwrx']:
        tail_len = 10
        opt_model = OptimizerRx(system_design, dft_len, cp_len, tail_len)
    elif system_design in ['CPW', 'WOLA']:
        opt_model = OptimizerTxRx()

    x, _ = opt_model.optimize(H_mat)
    print(x)
    np.save(f'./optimized_windows/{system_design}_{cp_len}.npy', x)

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
    """"""


class OptimizerTxRx:
    """"""

# EoF

