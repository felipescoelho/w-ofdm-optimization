"""Module with optmization tools for optimizal window design in w-OFDM.

luizfelipe.coelho@smt.ufrj.br
Jul 23, 2023
"""


__all__ = ['optimization_fun', 'reduce_variable_tx', 'reduce_variable_rx']


from .optimizers import optimization_fun
from .utils import reduce_variable_tx, reduce_variable_rx


# EoF
