"""wofdm_simulation.py

Script with w-OFDM simulation methods.

luizfelipe.coelho@smt.ufrj.br
Oct 10, 2023
"""


from .transmitter import (gen_add_redundancy_matrix, gen_idft_matrix,
                          gen_rc_window_tx)
from .receiver import (gen_circ_shift_matrix, gen_dft_matrix,
                       gen_overlap_and_add_matrix, gen_rc_window_rx,
                       gen_rm_redundancy_matrix)


class wOFDMSystem:
    """
    A class to simulate each system according to certain specifications.
    
    Attributes
    ----------

    """

    def __init__(self, system_design, dft_len, cp_len, tail_tx, tail_rx):
        """
        Parameters
        ----------
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
            self._shift_len = self.tail_tx
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
        self.rm_red_mat = gen_rm_redundancy_matrix(self.dft_len, 0, self.rm_len)
        self.overlap_add_mat = gen_overlap_and_add_matrix(self.dft_len, 0)
        self.circ_shift_mat = gen_circ_shift_matrix(self.dft_len, self.shift_len)
        self.dft_mat = gen_dft_matrix(self.dft_len)


