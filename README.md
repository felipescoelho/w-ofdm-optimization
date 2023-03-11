# w-ofdm-optimization

## Content:
- `channels/` foder contains the considered channels, previously generated
in MATLAB.
- `main_window_optimization.m` runs the window optimization process and
saves the results in the `optimized_windows/` folder. This will also generate
the `settingsData.mat` file, with general settings for the simulation and
specific settings for each system.
- `main_BER_calculation.m` calculates the BER for each system.
- `main_interference_calculation.m` calculates the interference power for
each system.
- 

## OBS.: Code was written using MATLAB R2021a to run on MATLAB R2017a.
