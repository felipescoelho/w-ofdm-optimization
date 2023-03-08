% main_channel_mask.m
%   This script will simulate a communication system using the channel mask
%   for out-of-band randiation (OOB) reduction, and will compare results
%   from systems using the optimized windows and the RC windows.
%
% Author: Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Mar 2, 2023
%

clc
clear
close all

fprintf('Starting main_channel_mask.m ... \n')

global numSubcar bitsPerSubcar symbolsPerTx offset

folderName = 'optimized_windows';
assert(isdir(folderName), 'Missing optimized_windows folder:', ...
    'Run main_window_optimization.m.\n')  %#ok
settingsFileName = 'settingsData.mat';
channelsFile = 'channels/vehA200channel2.mat';
channelLoader = load(channelsFile);
channels = channelLoader.vehA200channel2;
settingsLoader = load(settingsFileName);
numSubcar = settingsLoader.settingsData.generalSettings.numberSubcarriers;
bitsPerSubcar = settingsLoader.settingsData.generalSettings.bitsPerSubcarrier;
symbolsPerTx = settingsLoader.settingsData.generalSettings.symbolsPerTx;
ensemble = settingsLoader.settingsData.generalSettings.ensemble;
typeOFDM = 'wtx';
cpLength = 22;
snr = 20;
tailTx = settingsLoader.settingsData.(typeOFDM).tailTx;
tailRx = settingsLoader.settingsData.(typeOFDM).tailRx;
switch typeOFDM
    case 'wtx'
        csLength = tailTx;
        prefixRemovalLength = cpLength;
        circularShiftLength = 0;
end
rollOff = 10;
windowRCTx = transmitter_rc_window(numSubcar, cpLength, csLength, tailTx);
windowRCRx = receiver_rc_window(numSubcar, tailRx);
windowPath = strcat(folderName, '/', 'optimal_win_', typeOFDM, ...
    '_VehA200_', num2str(cpLength), 'CP.mat');
windowOptTx = load(windowPath).optimizedWindow;
offset = numSubcar/4;

for idx = 1:ensemble
    transmittedBits = randi([0 1], numSubcar*bitsPerSubcar/2, symbolsPerTx);
    transmittedSymbols = qammod(transmittedBits, 2^bitsPerSubcar, ...
        'InputType', 'bit', 'UnitAveragePower', true).';
    pilot = transmittedSymbols(1, :);
    [maskedRCwOFDM, RCwOFDM] = gen_tx_ofdm(transmittedSymbols, ...
        cpLength, csLength, windowRCTx, rollOff);
    [maskedOptwOFDM, optwOFDM] = gen_tx_ofdm(transmittedSymbols, ...
        cpLength, csLength, windowOptTx, rollOff);
    [numChannels, channelLength] = size(channels);
    for channel = 1:numChannels
        % RC window
        BERRC = run_sim(RCwOFDM, tailTx, cpLength, csLength, ...
            channels(channel, :), snr, tailRx, prefixRemovalLength, ...
            circularShiftLength, windowRCRx, pilot, transmittedBits);
        % Masked RC window
        BERmaskedRC = run_sim(maskedRCwOFDM, tailTx, cpLength, csLength, ...
            channels(channel, :), snr, tailRx, prefixRemovalLength, ...
            circularShiftLength, windowRCRx, pilot, transmittedBits);
        % Optimal window
        BEROpt = run_sim(optwOFDM, tailTx, cpLength, csLength, ...
            channels(channel, :), snr, tailRx, prefixRemovalLength, ...
            circularShiftLength, windowRCRx, pilot, transmittedBits);
        % Masked optimal window
        BERmaskedOpt = run_sim(maskedOptwOFDM, tailTx, cpLength, csLength, ...
            channels(channel, :), snr, tailRx, prefixRemovalLength, ...
            circularShiftLength, windowRCRx, pilot, transmittedBits);
    end
end


function [ber] = run_sim(signalTx, tailTx, cpLength, csLength, channel, ...
    snr, tailRx, prefixRemovalLength, circularShiftLength, windowRx, ...
    pilot, transmittedBits)
% Function to run simulation and estimate BER.
%

global symbolsPerTx bitsPerSubcar offset

signalRx = tx2rx(signalTx, tailTx, cpLength, csLength, channel, snr, ...
    tailRx, prefixRemovalLength);
symbolsRx = wofdm_rx(signalRx, tailRx, prefixRemovalLength, ...
    circularShiftLength, windowRx).';
channelEstimate = symbolsRx(1, offset+1:end-offset)./pilot;
estimatedChannelMatrix = repmat(channelEstimate, symbolsPerTx-1, 1);
estimatedSymbols = symbolsRx(2:end, offset+1:end-offset)./estimatedChannelMatrix;
recoveredBits = qamdemod(estimatedSymbols.', 2^bitsPerSubcar, ...
    'OutputType', 'bit', 'UnitAveragePower', true);
[~, ber] = biterr(transmittedBits(:, 2:end), recoveredBits);
end


function [maskedwOFDMTx, wOFDMTx] = gen_tx_ofdm(transmittedSymbols, ...
    cpLength, csLength, windowTx, rollOff)
% Function to generate transmitted OFDM symbols
%

global offset numSubcar symbolsPerTx

symbolsInOFDM = [zeros(symbolsPerTx, offset) transmittedSymbols ...
    zeros(symbolsPerTx, offset)];
invertTransformMatrix = dftmtx(numSubcar)'/numSubcar;
OFDMSymbols = invertTransformMatrix*symbolsInOFDM.';
redundancyMatrix = add_redundancy_matrix(cpLength, csLength);
OFDMSymbolsWithRedundancy = redundancyMatrix*ifftshift(OFDMSymbols, 1);
wOFDMTx = (windowTx*OFDMSymbolsWithRedundancy).';
transformMatrix = dftmtx(numSubcar+cpLength+csLength);
wOFDMSpectrum = transformMatrix*(wOFDMTx.');
windowRC = gen_raised_cosine(length(transmittedSymbols), rollOff, ...
    length(wOFDMSpectrum));
maskedSpectrum = diag(windowRC)*wOFDMSpectrum;
invertTransformMatrix = dftmtx(numSubcar+cpLength+csLength)' ...
    /(numSubcar+cpLength+csLength); 
maskedwOFDMTx = (invertTransformMatrix*maskedSpectrum).';
end


function [signalRx] = tx2rx(signalTx, tailTx, cpLength, csLength, ...
    channel, snr, tailRx, prefixRemovalLength)
% Function to perfrom the signal transmission and reception.
%

global symbolsPerTx numSubcar

matrixOverlap = signalTx(:, tailTx+1:end);
matrixOverlap(1:symbolsPerTx-1, end-tailTx+1:end) ...
    = signalTx(2:symbolsPerTx, 1:tailTx) ...
    + signalTx(1:symbolsPerTx-1, end-tailTx+1:end);
serialized = reshape(matrixOverlap.', 1, symbolsPerTx ...
    * (numSubcar+cpLength+csLength-tailTx));
transmittedSignal = [signalTx(1, 1:tailTx) serialized];
receivedSignal = add_wgn(conv(channel, transmittedSignal), snr);
validReceivedSignal = receivedSignal(1:end-length(channel)-tailTx+1);
signalRx = reshape(validReceivedSignal.', ...
    (numSubcar+tailRx+prefixRemovalLength), symbolsPerTx);
end


function [addRedundancyMatrix] = add_redundancy_matrix(cpLength, csLength)
% Function to generate the matrix that adds redundancy.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cpLength: Number of samples in cyclic prefix.
%   . csLength: Number of samples in cyclic suffix.
%
% - Output:
%   . redundancyMatrix: Matrix that adds the cyclic redundancy to the
%   system.
%

global numSubcar

identityPrefix = eye(cpLength);
zerosPrefix = zeros(cpLength, (numSubcar-cpLength));
identitySuffix = eye(csLength);
zerosSuffix = zeros(csLength, (numSubcar-csLength));
identitySubcarriers = eye(numSubcar);
addRedundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


function [windowRC] = gen_raised_cosine(windowLength, tailLength, ...
    totalLength)
% Function to generate the raised cosine window for the channel mask.
% This considers a centralized spectrum for the OFDM.
%

raisedCosineAxis = (-(tailLength+1)/2+1):1:((tailLength+1)/2-1);
raisedCosine = sin(pi/2*(.5+raisedCosineAxis/tailLength)).^2;
zeroPadd = zeros(1, (totalLength-windowLength-2*tailLength)/2);
windowRC = [zeroPadd raisedCosine ones(1, windowLength) ...
    fliplr(raisedCosine) zeroPadd];
end


function receivedSymbols = wofdm_rx(rxSignal, tailRx, ...
    prefixRemovalLength, circularShiftLength, windowRx)
% Function that performs the reception of w-OFDM symbols
%

global numSubcar

removeRedundancyMatrix = remove_redundancy_matrix(numSubcar, tailRx, ...
    prefixRemovalLength);
overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx);
circularShiftMatrix = circular_shift_matrix(numSubcar, circularShiftLength);
transformMatrix = dftmtx(numSubcar);
receivedSymbols = transformMatrix*circularShiftMatrix ...
    * overlapAddMatrix*windowRx*removeRedundancyMatrix*rxSignal;
end


function [overlapAddMatrix] = overlap_and_add_matrix(numSubcar, tailRx)
% Function to generate matrix that operates the overlap and add in the
% samples.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tails for the receiver
%   window.

identityHalfTail = eye(tailRx/2);
identitySubcarriers = eye(numSubcar-tailRx);
zerosHalfTail = zeros(tailRx/2);
zerosHalfTailSubcarriers = zeros(tailRx/2, numSubcar-tailRx);
zerosTailSubcarriers = zeros(numSubcar-tailRx, tailRx);
overlapAddMatrix = [zerosHalfTail identityHalfTail zerosHalfTailSubcarriers ...
    zerosHalfTail identityHalfTail; zerosTailSubcarriers ...
    identitySubcarriers zerosTailSubcarriers; ...
    identityHalfTail zerosHalfTail zerosHalfTailSubcarriers ...
    identityHalfTail zerosHalfTail];
end


function [circularShiftMatrix] = circular_shift_matrix(numSubcar, ...
    circularShiftLength)
% Function to generate matrix that operates the circular shift.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cicularShiftLength: Samples in the circular shift.
%
% - Output:
%   . circularShiftMatrix: A matrix capable of performing the circular
%   shift operation.

identityCircularShift = eye(circularShiftLength);
identitySubcarriersMinusCircularShift = eye(numSubcar- ...
    circularShiftLength);
zerosSubcarriersMinusCircularShift = zeros(circularShiftLength, ...
    numSubcar-circularShiftLength);
circularShiftMatrix = [zerosSubcarriersMinusCircularShift' ...
    identitySubcarriersMinusCircularShift; identityCircularShift ...
    zerosSubcarriersMinusCircularShift];
end


function [removeRedundancyMatrix] = remove_redundancy_matrix(numSubcar, ...
    tailRx, prefixRemovalLength)
% Function to generate matrix to remove redundancy.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail at the receiver
%   window.
%   . prefixRemoval: Number of samples to be removed at the reception.
%

removeRedundancyMatrix = [zeros(numSubcar+tailRx, prefixRemovalLength) ...
    eye(numSubcar+tailRx)];
end


function [windowTxRC] = transmitter_rc_window(numSubcar, cpLength, ...
    csLength, tailTx)
% Function to generate the Raised Cosine transmitter window.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cpLength: Number of samples in cyclic prefix.
%   . csLength: Number of samples in cyclic suffix.
%   . tailTx: Number of samples in rise and fall tail for the transmitter
%   window.
%
% - Output:
%   . windowTxRC: Matrix with diagonal equivalent to the transmitted raised
%   cosine window.
%

raisedCosineAxis = (-(tailTx+1)/2 + 1):1:((tailTx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailTx)).^2;
onesTransmitted = ones(1, numSubcar+cpLength+csLength ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
windowTxRC = diag(windowVector);
end


function [noisySignal] = add_wgn(inputSignal, snr)
% Function to add white Gaussian noise (WGN) to a signal, according to a SNR.
%
% - Input:
%   . inputSignal: Signal to add WGN
%   . snr: Signal to noise ration to be adjusted, in dB.
%

if iscolumn(inputSignal)
    inputSignal = inputSingal.';
end

signalPower = inputSignal*inputSignal' / length(inputSignal);
noise = randn(1, length(inputSignal)) + 1j*randn(1, length(inputSignal));
noisePower = noise*noise' / length(noise);
adjustedNoise = sqrt((signalPower * 10^(-.1*snr) / noisePower))*noise;
noisySignal = inputSignal + adjustedNoise;
end


function [windowRx] = receiver_rc_window(numSubcar, tailRx)
% Function to generate the Raised Cosine receiver window.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail for the receiver
%   window.
%

raisedCosineAxis = (-(tailRx+1)/2+1):1:((tailRx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailRx)).^2;
onesReceived = ones(1, numSubcar-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
windowRx = diag(windowVector);
end


% EoF
