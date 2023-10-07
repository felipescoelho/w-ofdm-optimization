"""Script to generate Rayleigh fading distribution.

Author: Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
Jul 3, 2023
"""


import numpy as np


def rayleigh_fading_jakes(no_oscillators:int, doppler_freq:float,
                          sampling_freq : float):
    """This function generates a Rayleigh fading waveform. We use
    Jakes' model for fading channel.
    
    1. Each oscilator corresponds has a frequency at the first quadrant
    of the unitary circle. The frequencies are obtained by the cosinus
    of evenly spaced angles in [0, pi/2).

    Parameters
    ----------
    no_oscillators : int
        Number of oscillators (>= 15 is good practice).
    doppler_freq : float
        Maximum Doppler frequency.
    sampling_freq : float
        Sampling rate.
    """

    time_index = np.arange(0, 1, 1/sampling_freq)
    real_part = np.sqrt(2)*np.cos(2*np.pi*doppler_freq*time_index)
    imag_part = 0*time_index
    for oscillator_idx in range(no_oscillators):
        # Generate oscillators:
        oscillator_freq = doppler_freq*np.cos(
            2*np.pi*(oscillator_idx+1)/(4*no_oscillators+2)
        )
        random_phase = np.pi*np.random.randn(1)
        oscillator = np.cos(2*np.pi*oscillator_freq*time_index + random_phase)
        # Gain:
        gain_angle = (oscillator_idx+1)*np.pi/(no_oscillators+1)
        # Calculate and sum oscillators:
        real_part += 2*np.cos(gain_angle)*oscillator
        imag_part += 2*np.sin(gain_angle)*oscillator
    # Scale:
    real_part /= np.sqrt(2*no_oscillators)
    imag_part /= np.sqrt(2*(no_oscillators+1))

    return real_part + 1j*imag_part


def rayleigh_fading_gmeds_1(no_oscillators:int, doppler_freq:float,
                            no_channels:int, sampling_freq:float,
                            no_waveforms:int):
    """Method to generate Rayleigh fading waveform using the GMEDS_1
    algorithm described in:
    
    - M. Patzold, C. -X. Wang and B. O. Hogstad, "Two new
    sum-of-sinusoids-based methods for the efficient generation of
    multiple uncorrelated rayleigh fading waveforms," in IEEE
    Transactions on Wireless Communications, vol. 8, no. 6,
    pp. 3122-3131, June 2009, doi: 10.1109/TWC.2009.080769.
    
    Parameters
    ----------
    no_oscillators : int
        Number of oscillators (= 20 is good enough).
    doppler_freq : float
        Maximum Doppler frequency.
    sampling_freq : float
        Sampling rate.
    no_channels : int
        Number of channels in simulation.
    no_waveforms : int
        Total number of uncorrelated waveforms.
    
    Returns
    -------
    rayleigh_fading_waveforms : np.ndarray
        Rayleigh fading waveform.
    """

    time_axis = np.arange(0, no_channels, 1)/sampling_freq
    rayleigh_fading_waveforms = np.zeros((len(time_axis), no_waveforms),
                                         dtype=np.complex128)
    for wave_idx in range(no_waveforms):
        real_wave = np.zeros((len(time_axis),))
        imag_wave = np.zeros((len(time_axis),))
        for oscil_idx in range(no_oscillators):
            angle_rotation = (np.pi/(4*no_oscillators)) \
                * (wave_idx/(no_waveforms+2))  # Must be negative for imaginary
            angle_arrival = (np.pi/(2*no_oscillators))*(oscil_idx+.5)
            oscil_freq_real = doppler_freq*np.cos(angle_arrival+angle_rotation)
            oscil_freq_imag = doppler_freq*np.cos(angle_arrival-angle_rotation)
            real_wave += np.cos(2*np.pi*oscil_freq_real*time_axis
                                + np.pi*np.random.randn(1))
            imag_wave += np.cos(2*np.pi*oscil_freq_imag*time_axis
                                + np.pi*np.random.randn(1))
        rayleigh_wave = np.sqrt(2/no_oscillators)*(real_wave + 1j*imag_wave)
        rayleigh_fading_waveforms[:, wave_idx] = rayleigh_wave

    return rayleigh_fading_waveforms


if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    from scipy.constants import speed_of_light

    kmh2ms = lambda x: (1000*x) / (60*60)
    v = 100  # Speed in km/h
    Fs = 1/(200*1e-9)  # Sampling rate
    Fc = 2*1e9  # Carrier frequency
    doppler_freq = kmh2ms(v) * Fc / speed_of_light
    m = rayleigh_fading_gmeds_1(20, doppler_freq, Fs, 4)

    mm = m[0:2000, 0]
    N = 2**12
    MM = np.fft.fftshift(np.fft.fft(mm, N))
    ff = np.linspace(-Fs/2, Fs/2, N)
    print(doppler_freq)

    print(m[:, 0])
    print(m.shape)

    plt.figure()
    plt.plot(ff, 20*np.log10(abs(MM)))
    plt.show()

# EoF
