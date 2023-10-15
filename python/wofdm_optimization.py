"""wofdm_optimization.py

Main script for w-ofdm optimization.

Oct 7, 2023
"""


import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from multiprocessing import cpu_count, Pool
from channel_model import gen_chan
from ofdm_utils import simulation_fun


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', type=str,
                        help='Operation mode. Can be gen_chan, run_opt or ' \
                            + 'run_sim.')
    parser.add_argument('--no_channels', type=int, default=250,
                        help='Number of independent channels generated.')
    parser.add_argument('-cp', '--cp_length', type=str,
                        default='10,12,14,16,18,20,22,24,26,28,30,32',
                        help='Lengths of cyclic-prefix to be tested.')
    parser.add_argument('-s', '--systems', type=str,
                        default='wtx,CPwtx,wrx,CPwrx,CPW,WOLA',
                        help='w-OFDM systems that will be tested.')
    parser.add_argument('--parallel', action=argparse.BooleanOptionalAction,
                        help='Add the argument parallel to run processes in ' \
                            + 'parallel.')
    parser.add_argument('--channel_path', type=str, default='channels',
                        help='Path to folder with channel. (save or load)')
    parser.add_argument('--window_path', type=str, default='optimized_windows',
                        help='Path to folder with optimized windows.')
    parser.add_argument('--channel_standard', type=str, default='vehicularA',
                        help='Channel model standard.')
    parser.add_argument('-mc', '--monte_carlo', type=int, default=1,
                        help='Number of repetitions in Monte Carlo process.')
    parser.add_argument('--snr', type=str, default='-21,51,3',
                        help='SNR for the simulation ([start, ]stop, [step).')
    parser.add_argument('--no_symbols', type=int, default=16,
                        help='Number of symbols in a frame.')

    return parser


if __name__ == '__main__':

    parser = arg_parser()
    args = parser.parse_args()
    channel_data_folder = args.channel_path
    channel_standard = args.channel_standard
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
        from optimization_tools import optimization_fun
        os.makedirs(args.window_path, exist_ok=True)
        # Simulation settings:
        dft_len = 256
        cp_list = [int(cp_len) for cp_len in args.cp_length.split(',')]
        sys_list = args.systems.split(',')
        data_list = [(sys, dft_len, cp, channel_path) for sys in sys_list for
                     cp in cp_list]
        if args.parallel:
            with Pool(cpu_count()) as pool:
                a = pool.map(optimization_fun, data_list)
        else:
            a = [optimization_fun(data) for data in data_list]
    elif args.mode == 'run_sim':
        # Simulation settings:
        dft_len = 256
        cp_list = [int(cp_len) for cp_len in args.cp_length.split(',')]
        sys_list = args.systems.split(',')
        snr_arr = np.arange(*[int(val) for val in args.snr.split(',')])
        tail_tx_fun = lambda x,y: x if y in ['CPW', 'WOLA', 'CPwtx', 'wtx'] \
            else 0
        tail_rx_fun = lambda x,y: x if y in ['CPW', 'WOLA', 'CPwrx', 'wrx'] \
            else 0
        data_list = [(sys, dft_len, cp, tail_tx_fun(8, sys),
                      tail_rx_fun(10, sys), channel_path, args.window_path,
                      args.monte_carlo, snr_arr, args.no_symbols)
                     for sys in sys_list for cp in cp_list]
        if args.parallel:
            with Pool(cpu_count()) as pool:
                pool.map(simulation_fun, data_list)
        else:
            [simulation_fun(data) for data in data_list]
    elif args.mode == 'gen_figs.':
        pass
    else:
        print('Mode is not defined, be sure to use gen_chan, run_opt, or ' \
              + 'run_sim.')


# EoF
