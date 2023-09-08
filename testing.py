"""Script for general code tests.

luizfeliep.coelho@smt.ufrj.br
Aug 2, 2023
"""


import numpy as np
import time


if __name__ == '__main__':
    H = np.array(((1, 0, -1, 0), (0, 1, 0, -1), (-1, 0, 1, 0), (0, -1, 0, 1)),
                 dtype=np.float64, ndmin=2)
    x = np.array((.400002, .799999, 1.000001, 2.000003), dtype=np.float64,
                 ndmin=2).T
    print(np.sqrt(x.T@H@x))