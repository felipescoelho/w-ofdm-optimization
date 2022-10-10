% main_BER_calculation.m
%
% Luiz Felipe Coelho -- luizfelipe.coelho@smt.ufrj.br
% Out 7, 2022
%


clc
clear
close all


fprintf('Staritng main_BER_calculation.m ...\n\n')
fprintf('This script runs a BER calculation for the proposed systems,\n')
fprintf('comparing the results for systems using the RC window with\n')
fpritnf('those using the optimized window.\n\n')

% Definitions
optimizedWindowsFolder = 'optimized_windows';
windowFiles = dir(fullfile(optimizedWindowsFolder));
snrValues = linspace(-20, 50, 30);

for fileIndex = 1:length(windowFiles)
    if windowFiles(fileIndex).isdir
        continue
    end
    fprintf('Working on file %s.\n', windowFiles(fileIndex).name)
    windowInfo = split(windowFiles(fileIndex).name, '_');
    typeOFDM = windowInfo{3};
    cyclicPrefix = str2double(windowInfo{end}(1:end-6));
    fprintf('This file contains the optmized window for %s-OFDM.\n', ...
        typeOFDM)
    fprintf('With %u samples in cyclicprefix.\n', cyclicPrefix)
    load([optimizedWindowsFolder '/' windowFiles(fileIndex).name])
    switch typeOFDM
        case {'wtx', 'CPwtx'}
        case {'wrx', 'CPwrx'}
        case {'CPW', 'WOLA'}
    end
end


function run_simulation(ensemble, symbolsPerTx, bitsPerSubcarrier, ...
    numberSubcarriers)
% Function that performs the transmission and reception of the OFDM
% symbols.
%
% - Input:
%   . ensemble: Number of repetition for Monte Carlo
%   . symbolsPerTx: Number of symbols per transmission
%   . bitsPerSubcarrier: Number of bits in a subcarrier
%   . numberSubcarriers: Number of subcarriers in OFDM symbol
%   . tailTx: Number of samples in rise and fall tail for Tx
%

for runID = 1:ensemble
    transmittedSymbols = qammod(randi([0 1], numberSubcarriers,...
        symbolsPerTx), 2^bitsPerSubcarrier, 'InputType', 'bit', ...
        'UnitAveragePower', true);
    pilotSymbols = transmittedSymbols(1, :);
    windowedOFDM = wofdm_tx(transmittedSymbols, numberSubcarriers, ...
        cyclicPrefix, cyclicSuffix, windowTx);
    matrixOverlap = windowedOFDM(:, tailTx+1:end);
    matrixOverlap(1:symbolsPerTx-1, end-tailTx+1:end) ...
        = windowedOFDM(2:symbolsPerTx, 1:tailTx) ...
        + windowedOFDM(1:symbolsPerTx-1, end-tailTx+1:end);
    
end
end


function windowedOFDM = wofdm_tx(transmittedSymbols, numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix, windowTx)
% Function that prepares the symbols from the digital modulation to be
% transmitted as OFDM symbols.
%
% - Input:
%   . transmittedSymbols: Symbols from digital modulation
%   . numberSubcarriers: Number of subcarriers in OFDM symbol
%   . cyclicPrefix: Number of elements in cyclic prefix
%   . cyclicSuffix: Number of elements in cyclic suffix
%   . windowTx: Transmitter window
%

invertTransformMatrix = dftmtx(numberSubcarriers)'/numberSubcarriers;
ofdmSymbols = invertTransformMatrix*transmittedSymbols.';
redundancyMatrix = add_redundancy_matrix(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix);
ofdmSymbolsWithRedundancy = redundancyMatrix*ofdmSymbols;
windowedOFDM = (windowTx*ofdmSymbolsWithRedundancy).';
end


function transmitterWindow = transmitter_rc_window(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix, tailTx)
% Function to generate the Raised Cosine transmitter window.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cyclicPrefix: Number of samples in cyclic prefix.
%   . cyclicSuffix: Number of samples in cyclic suffix.
%   . tailTx: Number of samples in rise and fall tail for the transmitter
%   window.
%
% - Output:
%   . transmitterWindow: Matrix with diagonal equivalent to the transmitted
%   window.

raisedCosineAxis = (-(tailTx+1)/2 + 1):1:((tailTx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailTx)).^2;
onesTransmitted = ones(1, numberSubcarriers+cyclicPrefix+cyclicSuffix ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
transmitterWindow = diag(windowVector);
end


function receiverWindow = receiver_rc_window(numberSubcarriers, tailRx)
% Function to generate the Raised Cosine receiver window.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail for the receiver
%   window.

raisedCosineAxis = (-(tailRx+1)/2+1):1:((tailRx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailRx)).^2;
onesReceived = ones(1, numberSubcarriers-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
receiverWindow = diag(windowVector);
end


function redundancyMatrix = add_redundancy_matrix(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix)
% Function to generate the matrix that adds redundancy.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cyclicPrefix: Number of samples in cyclic prefix.
%   . cyclicSuffix: Number of samples in cyclic suffix.
%
% - Output:
%   . redundancyMatrix: Matrix that adds the cyclic redundancy to the
%   system.

identityPrefix = eye(cyclicPrefix);
zerosPrefix = zeros(cyclicPrefix, (numberSubcarriers-cyclicPrefix));
identitySuffix = eye(cyclicSuffix);
zerosSuffix = zeros(cyclicSuffix, (numberSubcarriers-cyclicSuffix));
identitySubcarriers = eye(numberSubcarriers);
redundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


% EoF

