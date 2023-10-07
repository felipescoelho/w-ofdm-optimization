"""optimizers.py

Script with class for the window optimization process.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


import sympy as sym
from ofdm_utils import (
    gen_idft_matrix, gen_rc_window_tx, gen_add_redundancy_matrix,
    gen_rm_redundancy_matrix, gen_overlap_and_add, gen_dft_matrix,
    gen_circ_shift_matrix, gen_channel_tensor
)


class OptimizerTx:
    """
    A class for the optmization process
    
    """

    def __init__(self, **kwargs):
        """Constructor method.
        
        Parameters
        ----------
        Keyword args:
            sysmte_desing : str
                w-OFDM system design.
            dft_len : int
                Length of DFT.
            cp_len : int
                Number of samples in cyclic prefix (CP).
            cs_len : int
                Number of samples in cyclic suffix (CS).
            tail_len : int
                Number of samples in window tail.
            rm_len : int
                Number of samples to skip before recovering useful block.
            shift_len : int
                Number of samples for the circular shift.
        """

        self.name = kwargs['system_design']
        self.dft_len = kwargs['dft_len']
        self.cp_len = kwargs['cp_len']
        self.cs_len = kwargs['cs_len']
        self.tail_len = kwargs['tail_len']
        if self.name == 'wtx-OFDM':
            self.rm_len = self.cp_len
            self.shift_len = 0
        elif self.name == 'CPwtx-OFDM':
            self.rm_len = self.cp_len - self.tail_len
            self.shift_len = self.tail_len
        # Set matrices:
        # Transmitter:
        self.idft_mat = gen_idft(self.dft_len)
        self.add_red_mat = gen_add_redundancy(self.dft_len, self.cp_len,
                                                  self.cs_len)
        # Receiver:
        self.rm_red_mat = gen_rm_redundancy(self.dft_len, 0, self.rm_len)
        self.circ_shift_mat = gen_circ_shift(self.dft_len, self.shift_len)
        self.dft_mat = gen_dft(self.dft_len)

    def calculate_chann_matrices(self, channel_ir):
        """Method to calculate channel matrices for the system.
        
        Parameters
        ----------
        channel_ir : np.ndarray
            Channel's impulse response.
        
        Returns
        -------
        channel_list : list
            List with channel matrices.
        """
        channel_list = gen_channel_tensor(channel_ir, self.dft_len,
                                              self.tail_len, 0, self.rm_len,
                                              self.cp_len, self.cs_len)
        
        return channel_list

    def optimize_window(self, channel_list:list):
        """Method to calculate the window optmization."""
        mat_B = self.dft_mat*self.circ_shift_mat*self.rm_red_mat \
            * channel_list[0]
        mat_C = self.add_red_mat*self.idft_mat
        mat_D = self.add_red_mat*self.circ_shift_mat*self.rm_red_mat \
            * sum(channel_list[1:])
        mat_Q1 = sym.Matrix(
            self.dft_len+self.cp_len+self.cs_len,
            self.dft_len+self.cp_len+self.cs_len,
            lambda i, j: (mat_B[:, i].T * mat_B[:, j].C) \
                * (mat_C[i, :] * mat_C[j, :].H)
        )
        aux_Q2 = mat_D * mat_D.H
        mat_Q2 = sym.Matrix(
            self.dft_len+self.cp_len+self.cs_len,
            self.dft_len+self.cp_len+self.cs_len,
            lambda i, j: aux_Q2[i, j] if i == j else 0
        )

        mat_H = 2*(mat_Q1+mat_Q2)

        print(mat_H)


# EoF

