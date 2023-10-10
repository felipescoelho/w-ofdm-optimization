"""wofdm_optimization.py

Main script for w-ofdm optimization.

Oct 7, 2023
"""


import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from multiprocessing import cpu_count, Pool
from scipy.linalg import issymmetric
from optimization_tools import OptimizerTx
from channel_model import gen_chan


def parallel_optimization(data:tuple):
    """Method to parallellize the optimization process.
    
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
    
    system_design, dft_len, cp_len = data
    channel = np.load('./channels/vehicularA.npy')
    # channel_avg = np.mean(channel, 1)
    channel_avg = channel[:, 0]
    if system_design in ['wtx', 'CPwtx']:
        # print(channel_avg)
        tail_len = 8
        opt_model = OptimizerTx(system_design, dft_len, cp_len, tail_len)
        # print(opt_model.overlap_add_mat)
        chann_ten = opt_model.calculate_chann_matrices(channel_avg)
        reg = 1e-1
        # H_mat = (1-reg)*opt_model.gen_hessian(chann_ten) + reg*np.eye(9)
        H_mat = opt_model.gen_hessian(chann_ten) + reg*np.diagflat(np.linspace(
            1, 2, 9
        ))
        print(np.linalg.eigvals(H_mat))
        print(np.linalg.cond(H_mat))
        if not issymmetric(H_mat):
            print('Not symmetric!')
            H_mat = .5 * (H_mat + H_mat.T)
            print(np.linalg.eigvals(H_mat))
            print(np.linalg.cond(H_mat))

        opt_model.optimize(H_mat)

    return 'ok'


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', type=str,
                        help='Operation mode. Can be gen_chan, run_opt, or run_sim.')
    parser.add_argument('--no_channels', type=int, default=250,
                        help='Number of independent channels generated.')
    parser.add_argument('-cp', '--cp_length', type=str,
                        default='10,12,14,16,18,20,22,24,26,28,30,32',
                        help='Lengths of cyclic-prefix to be tested.')
    parser.add_argument('-s', '--systems', type=str,
                        default='wtx,CPwtx,wrx,CPwrx,CPW,WOLA',
                        help='w-OFDM systems that will be tested.')
    parser.add_argument('--parallel', action=argparse.BooleanOptionalAction,
                        help='Add the argument parallel to run processes in parallel.')

    return parser


if __name__ == '__main__':

    parser = arg_parser()
    args = parser.parse_args()

    channel_data_folder = './channels/'
    channel_standard = 'vehicularA'
    channel_path = os.path.join(channel_data_folder, channel_standard+'.npy')

    if args.mode == 'gen_chan':
        from scipy.constants import speed_of_light
        # Simulation settings:
        carrier_frequency = 2*1e9
        sample_period = 200*1e-9
        velocity = 100/3.6
        no_samples = 21  # Number of samples in tapped-delay line
        no_symbols = 16  # Number of symbols per frame
        dft_length = 256
        # Calculations:
        symbol_duration = dft_length*sample_period
        frame_duration = no_symbols*symbol_duration
        doppler_freq = (velocity/speed_of_light)*carrier_frequency
        # Generate channels:
        channel_tdl = np.zeros((no_samples, args.no_channels),
                               dtype=np.complex128)
        for idx in range(args.no_channels):
            channel_tdl[:, idx] = gen_chan(
                channel_standard, no_samples, doppler_freq, 1/sample_period,
                frame_duration, 1
            )[:, 0]
        # Save channels:
        os.makedirs(channel_data_folder, exist_ok=True)
        np.save(channel_path, channel_tdl, fix_imports=False)

    elif args.mode == 'run_opt':
        # Simulation settings:
        dft_len = 256
        cp_list = [int(cp_len) for cp_len in args.cp_length.split(',')]
        sys_list = args.systems.split(',')
        data_list = [(sys, dft_len, cp) for sys in sys_list for cp in cp_list]
        if args.parallel:
            with Pool(cpu_count()) as pool:
                a = pool.map(parallel_optimization, data_list)
        else:
            a = [parallel_optimization(data) for data in data_list]
            print([b for b in a])
    elif args.mode == 'run_sim':
        pass
    else:
        print('Mode is not defined, be sure to use gen_chan, run_opt, or run_sim.')

    # dft_len = 92
    # cp_len = 10
    # cs_len = 8
    # tail_len = 8

    # kwargs = {'system_design': 'wtx-OFDM', 'cp_len': cp_len, 'cs_len': cs_len,
    #           'dft_len': dft_len, 'tail_len': tail_len, 'rm_len': 20,
    #           'shift_len': 0}
    
    # wtx = OptimizerTx(**kwargs)
    # chan = np.random.randn(21) + 1j*np.random.randn(21)
    # test = wtx.calculate_chann_matrices(chan)
    # wtx.optimize_window(test)