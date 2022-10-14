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
    numberSubcarriers, cyclicPrefix, cyclicSuffix, windowTx, channel, ...
    snr, tailRx, windowRx)
% Function that performs the transmission and reception of the OFDM
% symbols.
%
% - Input:
%   . ensemble: Number of repetition for Monte Carlo
%   . symbolsPerTx: Number of symbols per transmission
%   . bitsPerSubcarrier: Number of bits in a subcarrier
%   . numberSubcarriers: Number of subcarriers in OFDM symbol
%   . tailTx: Number of samples in rise and fall tail for Tx
%   . channel: Channel for transmitted signal
%   . snr: Power for additive white Gaussian noise (AWGN), in dB

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
    serialized = reshape(matrixOverlap.', 1, symbolsPerTx ...
        * (numberSubcarriers+cyclicPrefix+cyclicSuffix-tailTx));
    transmittedSignal = [windowedOFDM(1, 1:tailTx) serialized];
    receivedSignal = add_wgn(conv(channel, transmittedSignal), snr);
    validReceivedSignal = receivedSignal(1:end-length(channel)-tailTx+1);
    parallelized = reshape(validReceivedSignal.', ...
        (numberSubcarriers+tailRx+prefixRemoval), symbolsPerTx);
    receivedSymbols = wofdm_rx(parallelized, numberSubcarriers, tailRx, ...
        prefixRemoval, circularShift, windowRx);
end
end


function noisySignal = add_wgn(inputSignal, snr)
% Function to add white Gaussian noise (WGN) to a signal, according to a SNR.
%
% - Input:
%   . inputSignal: Signal to add WGN
%   . snr: Signal to noise ration to be adjusted, in dB.
%

signalPower = norm(inputSignal, 2)/length(inputSignal);
unadjustedNoise = randn(1, length(signalPower)) ...
    + 1j*randn(1, length(signalPower));
adjustedNoise = sqrt((signalPower * 10^(-.1*snr) ...
    / (norm(unadjustedNoise, 2)/length(unadjustedNoise)))) ...
    * unadjustedNoise;
noisySignal = inputSignal + adjustedNoise;
end


function receivedSymbols = wofdm_rx(receivedSignal, numberSubcarriers, ...
    tailRx, prefixRemoval, circularShift, windowRx)
% Function that performs the reception of OFDM symbols
%

redundancyRemovalMatrix = remove_redundancy_matrix(numberSubcarriers, ...
    tailRx, prefixRemoval);
overlapAdd = overlap_and_add_matrix(numberSubcarriers, tailRx);
circularShiftMatrix = circular_shift_matrix(numberSubcarriers, ...
    circularShift);
transformMatrix = dftmtx(numberSubcarriers);
receivedSymbols = transformMatrix*circularShiftMatrix ...
    * overlapAdd*windowRx*redundancyRemovalMatrix*(receivedSignal.');

end


function circularShiftMatrix = circular_shift_matrix(numberSubcarriers, ...
    circularShift)
% Function to generate matrix that operates the circular shift.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cicularShift: Samples in the circular shift.
%
% - Output:
%   . circularShiftMatrix: A matrix capable of performing the circular
%   shift operation.

identityCircularShift = eye(circularShift);
identitySubcarriersMinusCircularShift = eye(numberSubcarriers- ...
    circularShift);
zerosSubcarriersMinusCircularShift = zeros(circularShift, ...
    numberSubcarriers-circularShift);
circularShiftMatrix = [zerosSubcarriersMinusCircularShift' ...
    identitySubcarriersMinusCircularShift; identityCircularShift ...
    zerosSubcarriersMinusCircularShift];
end


function overlapAdd = overlap_and_add_matrix(numberSubcarriers, tailRx)
% Function to generate matrix that operates the overlap and add in the
% samples.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tails for the receiver
%   window.

identityHalfTail = eye(tailRx/2);
identitySubcarriers = eye(numberSubcarriers-tailRx);
zerosHalfTail = zeros(tailRx/2);
zerosHalfTailSubcarriers = zeros(tailRx/2, numberSubcarriers-tailRx);
zerosTailSubcarriers = zeros(numberSubcarriers-tailRx, tailRx);
overlapAdd = [zerosHalfTail identityHalfTail zerosHalfTailSubcarriers ...
    zerosHalfTail identityHalfTail; zerosTailSubcarriers ...
    identitySubcarriers zerosTailSubcarriers; ...
    identityHalfTail zerosHalfTail zerosHalfTailSubcarriers ...
    identityHalfTail zerosHalfTail];
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


function redundancyRemove = remove_redundancy_matrix(numberSubcarriers, ...
    tailRx, prefixRemoval)
% Function to generate matrix to remove redundancy.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail at the receiver
%   window.
%   . prefixRemoval: Number of samples to be removed at the reception.

redundancyRemove = [zeros(numberSubcarriers+tailRx, prefixRemoval) ...
    eye(numberSubcarriers+tailRx)];
end


% EoF

