% window_optimization_new_new.m
%   This script performs the window optmization of all w-OFDM systems.
%
%   luizfelipe.coelho@smt.ufrj.br
%   Set 7, 2023
%


clc
clear
close all


addpath('functions')  % Best feature in MATLAB

global savePath channelPath
savePath = 'optimized_windows_new_new';
channelPath = './channels/vehA200channel2.mat';
if ~isdir(savePath)  %#ok
    mkdir(savePath)
end
assert(isfile(channelPath), 'Missing channel folder.\n')
settingsFileName = 'settingsData.mat';
fprintf('Generating settingsData structure array.\n')
settingsData = struct();
settingsData.generalSettings = struct('numberSubcarriers', 256, ...
   'bitsPerSubcarrier', 4, 'cyclicPrefix', 10:2:32, 'symbolsPerTx', ...
   16, 'ensemble', 50, 'snrValues', linspace(-20, 50, 30));
settingsData.wtx = struct('tailTx', 8, 'tailRx', 0);
settingsData.wrx = struct('tailTx', 0, 'tailRx', 10);
settingsData.WOLA = struct('tailTx', 8, 'tailRx', 10);
settingsData.CPW = struct('tailTx', 8, 'tailRx', 10);
settingsData.CPwtx = struct('tailTx', 8, 'tailRx', 0);
settingsData.CPwrx = struct('tailTx', 0, 'tailRx', 10);
fprintf('Saving settingsData structure array in %s\n', settingsFileName)
save(settingsFileName, 'settingsData')
fprintf('Successfully saved settingsData.\n')
fprintf('------------------------------------\n')
fprintf('Starting run ...\n\n')
fields = fieldnames(settingsData);
clear settingsData

for fieldIndex = 1:length(fields)
    fieldName = fields{fieldIndex};
    dataLoader = load(settingsFileName);
    dftLength = dataLoader.settingsData.generalSettings.numberSubcarriers;
    cpLengthVect = dataLoader.settingsData.generalSettings.cyclicPrefix;
    fprintf('Working on %s-OFDM ...\n', fieldName)
    for cpIndex = 1:length(cpLengthVect)
        cpLength = cpLengthVect(cpIndex);
        switch fieldName
            case 'wtx'
                tailTx = dataLoader.settingsData.(fieldName).tailTx;
                tailRx = dataLoader.settingsData.(fieldName).tailRx;
                csLength = tailTx;
                prefixRemovalLength = cpLength;
                circularShiftLength = 0;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', dftLength)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    dftLength, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength)
            case 'CPwtx'
                tailTx = dataLoader.settingsData.(fieldName).tailTx;
                tailRx = dataLoader.settingsData.(fieldName).tailRx;
                csLength = 0;
                prefixRemovalLength = cpLength - tailTx;
                circularShift = tailTx;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', dftLength)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    dftLength, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength)
            case 'wrx'
                tailTx = dataLoader.settingsData.(fieldName).tailTx;
                tailRx = dataLoader.settingsData.(fieldName).tailRx;
                csLength = tailRx/2;
                prefixRemovalLength = cpLength - tailRx/2;
                circularShiftLength = 0;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', dftLength)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    dftLength, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength)
            case 'CPwrx'
                tailTx = dataLoader.settingsData.(fieldName).tailTx;
                tailRx = dataLoader.settingsData.(fieldName).tailRx;
                csLength = 0;
                prefixRemovalLength = cpLength - tailRx;
                circularShiftLength = tailRx/2;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', dftLength)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    dftLength, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength)
            case 'CPW'
                tailTx = dataLoader.settingsData.(fieldName).tailTx;
                tailRx = dataLoader.settingsData.(fieldName).tailRx;
                csLength = tailTx + tailRx/2;
                prefixRemovalLength = cpLength - tailRx/2;
                circularShiftLength = 0;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', dftLength)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    dftLength, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength)
            case 'WOLA'
                tailTx = dataLoader.settingsData.(fieldName).tailTx;
                tailRx = dataLoader.settingsData.(fieldName).tailRx;
                csLength = tailTx;
                prefixRemovalLength = cpLength - tailRx;
                circularShiftLength = tailRx/2;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', dftLength)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    dftLength, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength)
        end
    end
end


% EoF
