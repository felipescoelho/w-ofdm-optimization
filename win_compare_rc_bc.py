"""Script for general code tests.

luizfeliep.coelho@smt.ufrj.br
Aug 2, 2023
"""


import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'helvet'
})


def rc_fun(n, fs):
    return np.sin(n*np.pi*fs)**2


if __name__ == '__main__':
    # Definitions
    fs = 1000
    n = np.linspace(0, 1, fs)
    width = 5.93
    height = width / ((1 + 5**.5)/2)
    font_size = 12
    fig_path = 'python/figures/window_compare_example_rc_bc.eps'

    # Calculate:
    y_rc = rc_fun(n, fs)
    y_bc = np.ones((fs,))
    Y_rc = np.fft.fftshift(np.fft.fft(y_rc, 2**11))
    Y_bc = np.fft.fftshift(np.fft.fft(y_bc, 2**11))
    f = np.fft.fftshift(np.fft.fftfreq(2**11, 1/fs))

    # Plot:
    # Time axis
    fig = plt.figure(figsize=(width, height))
    ax1 = fig.add_subplot(121)
    ax1.plot(n, y_rc, label='RC')
    ax1.plot(n, y_bc, label='boxcar')
    ax1.set_ylabel('Amplitude', fontsize=font_size)
    ax1.set_xlabel('$n$', fontsize=font_size)
    ax1.legend(fontsize=font_size)
    ax1.grid()
    # Freq axis
    ax2 = fig.add_subplot(122)
    ax2.plot(f, 10*np.log10(np.abs(Y_rc)), label='RC')
    ax2.plot(f, 10*np.log10(np.abs(Y_bc)), label='boxcar')
    ax2.set_ylim(-80, ax2.get_ylim()[1])
    ax2.set_ylabel('Magnitude, dB', fontsize=font_size)
    ax2.set_xlabel('Frequency, Hz', fontsize=font_size)
    ticklabels = ['', '$-f_s/2$', '$-f_s/4$', '0', '$f_s/4$', '$f_s/2$']
    ax2.set_xticklabels(ticklabels)
    ax2.grid()
    
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches='tight')

    plt.close()