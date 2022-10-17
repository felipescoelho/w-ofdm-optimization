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
fprintf('those using the optimized window.\n\n')

% Definitions
settingsFileName = 'settingsData.mat';
optimizedWindowsFolder = 'optimized_windows';
channelsFilePath = './channels/vehA200channel2.mat';
windowFiles = dir(fullfile(optimizedWindowsFolder));
berResultsFolder = 'ber_results';
dataLoader = load(settingsFileName);
snrValues = dataLoader.settingsData.generalSettings.snrValues;
if ~isfolder(berResultsFolder)
    mkdir(berResultsFolder)
end

for fileIndex = 1:length(windowFiles)
    if windowFiles(fileIndex).isdir
        continue
    end
    fprintf('Working on file %s.\n', windowFiles(fileIndex).name)
    windowInfo = split(windowFiles(fileIndex).name, '_');
    typeOFDM = windowInfo{3};
    cpLength = str2double(windowInfo{end}(1:end-6));
    fprintf('This file contains the optmized window for %s-OFDM.\n', ...
        typeOFDM)
    load([optimizedWindowsFolder '/' windowFiles(fileIndex).name])
    load(channelsFilePath, 'vehA200channel2')
    dataLoader = load(settingsFileName);
    numSubcar = dataLoader.settingsData.generalSettings.numberSubcarriers;
    bitsPerSubcar = dataLoader.settingsData.generalSettings.bitsPerSubcarrier;
    symbolsPerTx = dataLoader.settingsData.generalSettings.symbolsPerTx;
    ensemble = dataLoader.settingsData.generalSettings.ensemble;
    tailTx = dataLoader.settingsData.(typeOFDM).tailTx;
    tailRx = dataLoader.settingsData.(typeOFDM).tailRx;
    [csLength, prefixRemoval, circularShift] = calculate_parameters( ...
        typeOFDM, cpLength, tailTx, tailRx);
    receiverRCWindow = receiver_rc_window(numSubcar, tailRx);
    transmitterRCWindow = transmitter_rc_window(numSubcar, cpLength, ...
        csLength, tailTx);
    fileName = strcat('ber_', typeOFDM, '_', num2str(cpLength), 'CP');
    switch typeOFDM
        case {'wtx', 'CPwtx'}
            berSNR = zeros(length(snrValues), 1);
            berRCSNR = zeros(length(snrValues), 1);
            for snrIndex = 1:length(snrValues)
                ber = 0;
                berRC = 0;
                for channelIndex = 1:length(vehA200channel2)
                    ber = ber + run_simulation(ensemble, symbolsPerTx, ...
                        bitsPerSubcar, numSubcar, cpLength, ...
                        csLength, optimizedTransmitterWindow, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        receiverRCWindow, prefixRemoval, circularShift);
                    berRC = berRC + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, transmitterRCWindow, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        receiverRCWindow, prefixRemoval, circularShift);
                end
                berSNR(snrIndex) = ber/length(vehA200channel2);
                berRCSNR(snrIndex) = berRC/length(vehA200channel2);
            end
            save([berResultsFolder '/optimized_' fileName], 'berSNR')
            save([berResultsFolder '/rc_' fileName], 'berRCSNR')
        case {'wrx', 'CPwrx'}
            berSNR = zeros(length(snrValues), 1);
            berRCSNR = zeros(length(snrValues), 1);
            for snrIndex = 1:length(snrValues)
                ber = 0;
                berRC = 0;
                for channelIndex = 1:length(vehA200channel2)
                    ber = ber + run_simulation(ensemble, symbolsPerTx, ...
                        bitsPerSubcar, numSubcar, cpLength, ...
                        csLength, transmitterRCWindow, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        optimizedReceiverWindow, prefixRemoval, ...
                        circularShift);
                    berRC = berRC + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, transmitterRCWindow, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        receiverRCWindow, prefixRemoval, circularShift);
                end
                berSNR(snrIndex) = ber/length(vehA200channel2);
                berRCSNR(snrIndex) = berRC/length(vehA200channel2);
            end
            save([berResultsFolder '/optimized_' fileName], 'berSNR')
            save([berResultsFolder '/rc_' fileName], 'berRCSNR')
        case {'CPW', 'WOLA'}
            berSNRStep1A = zeros(length(snrValues), 1);
            berSNRStep1B = zeros(length(snrValues), 1);
            berSNRStep2A = zeros(length(snrValues), 1);
            berSNRStep2B = zeros(length(snrValues), 1);
            berSNRStep3A = zeros(length(snrValues), 1);
            berSNRStep3B = zeros(length(snrValues), 1);
            berRCSNR = zeros(length(snrValues), 1);
            for snrIndex = 1:length(snrValues)
                ber1A = 0;
                ber1B = 0;
                ber2A = 0;
                ber2B = 0;
                ber3A = 0;
                ber3B = 0;
                berRC = 0;
                for channelIndex = 1:length(vehA200channel2)
                    berRC = berRC + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, transmitterRCWindow, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        receiverRCWindow, prefixRemoval, circularShift);
                    % Case A:
                    ber1A = ber1A + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, ...
                        optimizedTransmitterStep1, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        receiverRCWindow, prefixRemoval, circularShift);
                    ber2A = ber2A + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, optimizedTransmitterStep1, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        optimizedReceiverStep2, prefixRemoval, ...
                        circularShift);
                    ber3A = ber3A + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, optimizedTransmitterStep3, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        optimizedReceiverStep2, prefixRemoval, ...
                        circularShift);
                    % Case B
                    ber1B = ber1B + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, transmitterRCWindow, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        optimizedReceiverStep1, prefixRemoval, ...
                        circularShift);
                    ber2B = ber2B + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, optimizedTransmitterStep2, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        optimizedReceiverStep1, prefixRemoval, ...
                        circularShift);
                    ber3B = ber3B + run_simulation(ensemble, ...
                        symbolsPerTx, bitsPerSubcar, numSubcar, ...
                        cpLength, csLength, optimizedTransmitterStep2, ...
                        vehA200channel2(channelIndex, :), ...
                        snrValues(snrIndex), tailTx, tailRx, ...
                        optimizedReceiverStep3, prefixRemoval, ...
                        circularShift);
                end
                berSNRStep1A(snrIndex) = ber1A/length(vehA200channel2);
                berSNRStep1B(snrIndex) = ber1B/length(vehA200channel2);
                berSNRStep2A(snrIndex) = ber2A/length(vehA200channel2);
                berSNRStep2B(snrIndex) = ber2B/length(vehA200channel2);
                berSNRStep3A(snrIndex) = ber3A/length(vehA200channel2);
                berSNRStep3B(snrIndex) = ber3B/length(vehA200channel2);
                berRCSNR(snrIndex) = berRC/length(vehA200channel2);
            end
            save([berResultsFolder '/optimized_' fileName], ...
                'berSNRStep1A', 'berSNRStep1B', 'berSNRStep2A', ...
                'berSNRStep2B', 'berSNRStep3A', 'berSNRStep3B')
            save([berResultsFolder '/rc_' fileName], 'berRCSNR')
        otherwise
            continue
    end
end


function ber = run_simulation(ensemble, symbolsPerTx, bitsPerSubcarrier, ...
    numberSubcarriers, cyclicPrefix, cyclicSuffix, windowTx, channel, ...
    snr, tailTx, tailRx, windowRx, prefixRemoval, circularShift)
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
    transmittedBits = randi([0 1], numberSubcarriers*bitsPerSubcarrier, ...
        symbolsPerTx);
    transmittedSymbols = qammod(transmittedBits, 2^bitsPerSubcarrier, ...
        'InputType', 'bit', 'UnitAveragePower', true).';
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
        prefixRemoval, circularShift, windowRx).';
    estimatedChannel = receivedSymbols(1, :)./pilotSymbols;
    estimatedChannelMatrix = repmat(estimatedChannel, symbolsPerTx-1, 1);
    estimatedSymbols = receivedSymbols(2:end, :)./estimatedChannelMatrix;
    recoveredBits = qamdemod(estimatedSymbols.', 2^bitsPerSubcarrier, ...
        'OutputType', 'bit', 'UnitAveragePower', true);
    
    [~, ber] = biterr(transmittedBits(:, 2:end), recoveredBits);
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
    * overlapAdd*windowRx*redundancyRemovalMatrix*receivedSignal;

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


function [cyclicSuffix, prefixRemoval, circularShift] = ...
    calculate_parameters(typeOFDM, cyclicPrefix, tailTx, tailRx)
% Function to calculate parameters for each w-OFDM system
%
% - Input:
%   . typeOFDM: w-OFDM system
%   . cyclicPrefix: Number of samples in the cyclic prefix
%   . tailTx: Number of samples in rise and fall tail for Tx
%   . tailRX: Number of samples in rise and fall tail for Rx

switch typeOFDM
    case 'wtx'
        cyclicSuffix = tailTx;
        prefixRemoval = cyclicPrefix;
        circularShift = 0;
    case 'wrx'
        cyclicSuffix = tailRx/2;
        prefixRemoval = cyclicPrefix - tailRx/2;
        circularShift = 0;
    case 'WOLA'
        cyclicSuffix = tailTx;
        prefixRemoval = cyclicPrefix - tailRx;
        circularShift = tailRx/2;
    case 'CPW'
        cyclicSuffix = tailTx + tailRx/2;
        prefixRemoval = cyclicPrefix - tailRx/2;
        circularShift = 0;
    case 'CPwtx'
        cyclicSuffix = 0;
        prefixRemoval = cyclicPrefix - tailTx;
        circularShift = tailTx;
    case 'CPwrx'
        cyclicSuffix = 0;
        prefixRemoval = cyclicPrefix - tailRx;
        circularShift = tailRx/2;
end
end


% EoF

