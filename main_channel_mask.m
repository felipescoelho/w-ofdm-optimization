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

folderName = 'optimized_windows';
assert(isdir(folderName), 'Missing optimized_windows folder:', ...
    'Run main_window_optimization.m.\n')  %#ok
windowFiles = dir(fullfile(folderName));
settingsFileName = 'settingsData.mat';
channelsFile = 'channels/vehA200channel2.mat';
berResultsFolder = 'ber_results';
if ~isdir(berResultsFolder)  %#ok
    mkdir(berResultsFolder)
end
channelMaskFolder = 'simulation_with_channel_mask';
resultsPath = [berResultsFolder '/' channelMaskFolder];
if ~isdir(resultsPath)  %#ok
    mkdir(resultsPath)
end

for fileIndex = 1:length(windowFiles)
    if windowFiles(fileIndex).isdir
        continue
    else
        dummy = split(windowFiles(fileIndex).name, '.');
        if strcmp(dummy{end}, 'log')
            continue
        end
    end
    fprintf('Working on file %s.\n', windowFiles(fileIndex).name)
    windowInfo = split(windowFiles(fileIndex).name, '_');
    typeOFDM = windowInfo{3};
    cpLength = str2double(windowInfo{end}(1:end-6));
    fprintf('%s-OFDM with %d samples in CP.\n', typeOFDM, cpLength)
    windowLoader = load([folderName '/' windowFiles(fileIndex).name]);
    channelLoader = load(channelsFile);
    channels = channelLoader.vehA200channel2;
    [numChannels, ~] = size(channels);
    dataLoader = load(settingsFileName);
    snrValues = dataLoader.settingsData.generalSettings.snrValues;
    numSubcar = dataLoader.settingsData.generalSettings.numberSubcarriers;
    offset = numSubcar/4;  % We add the double of this number as zeros in the DFT.
    rollOff = 10;
    bitsPerSubcar = dataLoader.settingsData.generalSettings.bitsPerSubcarrier;
    symbolsPerTx = dataLoader.settingsData.generalSettings.symbolsPerTx;
    ensemble = dataLoader.settingsData.generalSettings.ensemble;
    tailTx = dataLoader.settingsData.(typeOFDM).tailTx;
    tailRx = dataLoader.settingsData.(typeOFDM).tailRx;
    [csLength, prefixRemovalLength, circularShiftLength] ...
        = calculate_parameters(typeOFDM, cpLength, tailTx, tailRx);
    receiverRCWindow = rx_rc_window(numSubcar, tailRx);
    transmitterRCWindow = tx_rc_window(numSubcar, cpLength, csLength, ...
        tailTx);
    fileName = strcat('ber_', typeOFDM, '_', num2str(cpLength), 'CP');
    switch typeOFDM
        case {'wtx', 'CPwtx'}
            berOptSNR = zeros(length(snrValues), 1);
            berMaskedOptSNR = zeros(length(snrValues), 1);
            berRCSNR = zeros(length(snrValues), 1);
            berMaskedRCSNR = zeros(length(snrValues), 1);
            for snrIndex = 1:length(snrValues)
                berOpt = 0;
                berMaskedOpt = 0;
                berRC = 0;
                berMaskedRC = 0;
                for channelIndex = 1:numChannels
                    % Optimal
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, ...
                        windowLoader.optimizedWindow, receiverRCWindow, ...
                        channels(channelIndex, :), snrValues(snrIndex), ...
                        offset, prefixRemovalLength, ...
                        circularShiftLength, numSubcar, bitsPerSubcar, ...
                        symbolsPerTx, rollOff);
                    berOpt = berOpt + ber;
                    berMaskedOpt = berMaskedOpt + berMasked;
                    % RC
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, transmitterRCWindow, ...
                        receiverRCWindow, channels(channelIndex, :), ...
                        snrValues(snrIndex), offset, ...
                        prefixRemovalLength, circularShiftLength, ...
                        numSubcar, bitsPerSubcar, symbolsPerTx, rollOff);
                    berRC = berRC + ber;
                    berMaskedRC = berMaskedRC + berMasked;
                end
                berOptSNR(snrIndex) = berOpt/numChannels;
                berMaskedOptSNR(snrIndex) = berMaskedOpt/numChannels;
                berRCSNR(snrIndex) = berRC/numChannels;
                berMaskedRCSNR(snrIndex) = berMaskedRC/numChannels;
            end
            save_opt([resultsPath '/optimized_' fileName], berOptSNR)
            save_masked_opt([resultsPath '/masked_optimized_' fileName], ...
                berMaskedOptSNR)
            save_rc([resultsPath '/rc_' fileName], berRCSNR)
            save_masked_rc([resultsPath '/masked_rc_' fileName], ...
                berMaskedRCSNR)
        case {'wrx', 'CPwrx'}
            berOptSNR = zeros(length(snrValues), 1);
            berMaskedOptSNR = zeros(length(snrValues), 1);
            berRCSNR = zeros(length(snrValues), 1);
            berMaskedRCSNR = zeros(length(snrValues), 1);
            for snrIndex = 1:length(snrValues)
                berOpt = 0;
                berMaskedOpt = 0;
                berRC = 0;
                berMaskedRC = 0;
                for channelIndex = 1:numChannels
                    % Optimal
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, transmitterRCWindow, ...
                        windowLoader.optimizedWindow, ...
                        channels(channelIndex, :), snrValues(snrIndex), ...
                        offset, prefixRemovalLength, ...
                        circularShiftLength, numSubcar, bitsPerSubcar, ...
                        symbolsPerTx, rollOff);
                    berOpt = berOpt + ber;
                    berMaskedOpt = berMaskedOpt + berMasked;
                    % RC
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, transmitterRCWindow, ...
                        receiverRCWindow, channels(channelIndex, :), ...
                        snrValues(snrIndex), offset, ...
                        prefixRemovalLength, circularShiftLength, ...
                        numSubcar, bitsPerSubcar, symbolsPerTx, rollOff);
                    berRC = berRC + ber;
                    berMaskedRC = berMaskedRC + berMasked;
                end
                berOptSNR(snrIndex) = berOpt/numChannels;
                berMaskedOptSNR(snrIndex) = berMaskedOpt/numChannels;
                berRCSNR(snrIndex) = berRC/numChannels;
                berMaskedRCSNR(snrIndex) = berMaskedRC/numChannels;
            end
            save_opt([resultsPath '/optimized_' fileName], berOptSNR)
            save_masked_opt([resultsPath '/masked_optimized_' fileName], ...
                berMaskedOptSNR)
            save_rc([resultsPath '/rc_' fileName], berRCSNR)
            save_masked_rc([resultsPath '/masked_rc_' fileName], ...
                berMaskedRCSNR)
        case {'CPW', 'WOLA'}
            berA1SNR = zeros(length(snrValues), 1);
            berMaskedA1SNR = zeros(length(snrValues), 1);
            berA2SNR = zeros(length(snrValues), 1);
            berMaskedA2SNR = zeros(length(snrValues), 1);
            berA3SNR = zeros(length(snrValues), 1);
            berMaskedA3SNR = zeros(length(snrValues), 1);
            berB1SNR = zeros(length(snrValues), 1);
            berMaskedB1SNR = zeros(length(snrValues), 1);
            berB2SNR = zeros(length(snrValues), 1);
            berMaskedB2SNR = zeros(length(snrValues), 1);
            berB3SNR = zeros(length(snrValues), 1);
            berMaskedB3SNR = zeros(length(snrValues), 1);
            berRCSNR = zeros(length(snrValues), 1);
            berMaskedRCSNR = zeros(length(snrValues), 1);
            for snrIndex = 1:length(snrValues)
                berA1 = 0;
                berMaskedA1 = 0;
                berA2 = 0;
                berMaskedA2 = 0;
                berA3 = 0;
                berMaskedA3 = 0;
                berB1 = 0;
                berMaskedB1 = 0;
                berB2 = 0;
                berMaskedB2 = 0;
                berB3 = 0;
                berMaskedB3 = 0;
                berRC = 0;
                berMaskedRC = 0;
                for channelIndex = 1:numChannels
                    % Case A - step 1
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, ...
                        windowLoader.optimizedWindowCaseAStep1, ...
                        receiverRCWindow, ...
                        channels(channelIndex, :), snrValues(snrIndex), ...
                        offset, prefixRemovalLength, ...
                        circularShiftLength, numSubcar, bitsPerSubcar, ...
                        symbolsPerTx, rollOff);
                    berA1 = berA1 + ber;
                    berMaskedA1 = berMaskedA1 + berMasked;
                    % Case A - step 2
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, ...
                        windowLoader.optimizedWindowCaseAStep1, ...
                        windowLoader.optimizedWindowCaseAStep2, ...
                        channels(channelIndex, :), ...
                        snrValues(snrIndex), offset, ...
                        prefixRemovalLength, circularShiftLength, ...
                        numSubcar, bitsPerSubcar, symbolsPerTx, rollOff);
                    berA2 = berA2 + ber;
                    berMaskedA2 = berMaskedA2 + berMasked;
                    % Case A - step 3
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, ...
                        windowLoader.optimizedWindowCaseAStep3, ...
                        windowLoader.optimizedWindowCaseAStep2, ...
                        channels(channelIndex, :), snrValues(snrIndex), ...
                        offset, prefixRemovalLength, ...
                        circularShiftLength, numSubcar, bitsPerSubcar, ...
                        symbolsPerTx, rollOff);
                    berA3 = berA3 + ber;
                    berMaskedA3 = berMaskedA3 + berMasked;
                    % Case B - step 1
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, transmitterRCWindow, ...
                        windowLoader.optimizedWindowCaseBStep1, ...
                        channels(channelIndex, :), ...
                        snrValues(snrIndex), offset, ...
                        prefixRemovalLength, circularShiftLength, ...
                        numSubcar, bitsPerSubcar, symbolsPerTx, rollOff);
                    berB1 = berB1 + ber;
                    berMaskedB1 = berMaskedB1 + berMasked;
                    % Case B - step 2
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, ...
                        windowLoader.optimizedWindowCaseBStep2, ...
                        windowLoader.optimizedWindowCaseBStep1, ...
                        channels(channelIndex, :), snrValues(snrIndex), ...
                        offset, prefixRemovalLength, ...
                        circularShiftLength, numSubcar, bitsPerSubcar, ...
                        symbolsPerTx, rollOff);
                    berB2 = berB2 + ber;
                    berMaskedB2 = berMaskedB2 + berMasked;
                    % Case B - step 3
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, ...
                        windowLoader.optimizedWindowCaseBStep2, ...
                        windowLoader.optimizedWindowCaseBStep3, ...
                        channels(channelIndex, :), ...
                        snrValues(snrIndex), offset, ...
                        prefixRemovalLength, circularShiftLength, ...
                        numSubcar, bitsPerSubcar, symbolsPerTx, rollOff);
                    berB3 = berB3 + ber;
                    berMaskedB3 = berMaskedB3 + berMasked;
                    % RC
                    [berMasked, ber] = run_sim_mc(ensemble, cpLength, ...
                        csLength, tailTx, tailRx, transmitterRCWindow, ...
                        receiverRCWindow, channels(channelIndex, :), ...
                        snrValues(snrIndex), offset, ...
                        prefixRemovalLength, circularShiftLength, ...
                        numSubcar, bitsPerSubcar, symbolsPerTx, rollOff);
                    berRC = berRC + ber;
                    berMaskedRC = berMaskedRC + berMasked;
                end
                berA1SNR(snrIndex) = berA1/numChannels;
                berMaskedA1SNR(snrIndex) = berMaskedA1/numChannels;
                berA2SNR(snrIndex) = berA2/numChannels;
                berMaskedA2SNR(snrIndex) = berMaskedA2/numChannels;
                berA3SNR(snrIndex) = berA3/numChannels;
                berMaskedA3SNR(snrIndex) = berMaskedA3/numChannels;
                berB1SNR(snrIndex) = berB1/numChannels;
                berMaskedB1SNR(snrIndex) = berMaskedB1/numChannels;
                berB2SNR(snrIndex) = berB2/numChannels;
                berMaskedB2SNR(snrIndex) = berMaskedB2/numChannels;
                berB3SNR(snrIndex) = berB3/numChannels;
                berMaskedB3SNR(snrIndex) = berMaskedB3/numChannels;
                berRCSNR(snrIndex) = berRC/numChannels;
                berMaskedRCSNR(snrIndex) = berMaskedRC/numChannels;
            end
            save_opt_multi([resultsPath '/optimized_' fileName], ...
                berA1SNR, berA2SNR, berA3SNR, berB1SNR, berB2SNR, berB3SNR)
            save_masked_opt_multi( ...
                [resultsPath  '/masked_optimized_' fileName], ...
                berMaskedA1SNR, berMaskedA2SNR, berMaskedA3SNR, ...
                berMaskedB1SNR, berMaskedB2SNR, berMaskedB3SNR)
            save_rc([resultsPath '/rc_' fileName], berRCSNR)
            save_masked_rc([resultsPath '/masked_rc_' fileName], ...
                berMaskedRCSNR)
        otherwise
            fprintf('Something is wrong!!\n')
    end
end


function save_opt(fileName, berSNR)
% SAVE_OPT  Save BER results for systems using optimal window.
%   SAVE_OPT(fileName, berSNR)  saves berSNR into fileName.

save(fileName, 'berSNR')
end


function save_masked_opt(fileName, berMaskedSNR)

save(fileName, 'berMaskedSNR')
end


function save_rc(fileName, berRCSNR)

save(fileName, 'berRCSNR')
end


function save_masked_rc(fileName, berMaskedRCSNR)

save(fileName, 'berMaskedRCSNR')
end


function save_opt_multi(fileName, berSNRStep1A, berSNRStep2A, ...
    berSNRStep3A, berSNRStep1B, berSNRStep2B, berSNRStep3B)
% Function to save multiple optimized windows, for WOLA and CPW.
%

save(fileName, 'berSNRStep1A', 'berSNRStep1B', 'berSNRStep2A', ...
    'berSNRStep2B', 'berSNRStep3A', 'berSNRStep3B')
end


function save_masked_opt_multi(fileName, berMaskedSNRStep1A, ...
    berMaskedSNRStep2A, berMaskedSNRStep3A, berMaskedSNRStep1B, ...
    berMaskedSNRStep2B, berMaskedSNRStep3B)
% Function to save multiple optimized windows, for WOLA and CPW.
%

save(fileName, 'berMaskedSNRStep1A', 'berMaskedSNRStep1B', 'berMaskedSNRStep2A', ...
    'berMaskedSNRStep2B', 'berMaskedSNRStep3A', 'berMaskedSNRStep3B')
end


function [berMasked, ber] = run_sim_mc(ensemble, cpLength, csLength, ...
    tailTx, tailRx, windowTx, windowRx, channel, snr, offset, ...
    prefixRemovalLength, circularShiftLength, numSubcar, bitsPerSubcar, ...
    symbolsPerTx, rollOff)
% Function to run simulation in a Monte Carlo process.
%

ber = 0;
berMasked = 0;
for idx = 1:ensemble
    transmittedBits = randi([0 1], numSubcar*bitsPerSubcar/2, symbolsPerTx);
    transmittedSymbols = qammod(transmittedBits, 2^bitsPerSubcar, ...
        'InputType', 'bit', 'UnitAveragePower', true).';
    pilot = transmittedSymbols(1, :);
    [maskedwOFDM, wOFDM] = gen_tx_ofdm(transmittedSymbols, cpLength, ...
        csLength, windowTx, rollOff, symbolsPerTx, numSubcar, offset);
    ber = ber + run_sim(wOFDM, tailTx, cpLength, csLength, channel, ...
        snr, tailRx, prefixRemovalLength, circularShiftLength, windowRx, ...
        pilot, transmittedBits, symbolsPerTx, bitsPerSubcar, offset, numSubcar);
    berMasked = berMasked + run_sim(maskedwOFDM, tailTx, cpLength, ...
        csLength, channel, snr, tailRx, prefixRemovalLength, ...
        circularShiftLength, windowRx, pilot, transmittedBits, ...
        symbolsPerTx, bitsPerSubcar, offset, numSubcar);
end
ber = ber/ensemble;
berMasked = berMasked/ensemble;
end


function [ber] = run_sim(signalTx, tailTx, cpLength, csLength, channel, ...
    snr, tailRx, prefixRemovalLength, circularShiftLength, windowRx, ...
    pilot, transmittedBits, symbolsPerTx, bitsPerSubcar, offset, numSubcar)
% Function to run simulation and estimate BER.
%

signalRx = tx2rx(signalTx, tailTx, cpLength, csLength, channel, snr, ...
    tailRx, prefixRemovalLength, symbolsPerTx, numSubcar);
symbolsRx = wofdm_rx(signalRx, tailRx, numSubcar, prefixRemovalLength, ...
    circularShiftLength, windowRx).';
channelEstimate = symbolsRx(1, offset+1:end-offset)./pilot;
estimatedChannelMatrix = repmat(channelEstimate, symbolsPerTx-1, 1);
estimatedSymbols = symbolsRx(2:end, offset+1:end-offset)./estimatedChannelMatrix;
recoveredBits = qamdemod(estimatedSymbols.', 2^bitsPerSubcar, ...
    'OutputType', 'bit', 'UnitAveragePower', true);
[~, ber] = biterr(transmittedBits(:, 2:end), recoveredBits);
end


function [maskedwOFDMTx, wOFDMTx] = gen_tx_ofdm(transmittedSymbols, ...
    cpLength, csLength, windowTx, rollOff, symbolsPerTx, numSubcar, offset)
% Function to generate transmitted OFDM symbols
%

symbolsInOFDM = [zeros(symbolsPerTx, offset) transmittedSymbols ...
    zeros(symbolsPerTx, offset)];
invertTransformMatrix = dftmtx(numSubcar)'/numSubcar;
OFDMSymbols = invertTransformMatrix*ifftshift(symbolsInOFDM.', 1);
redundancyMatrix = add_redundancy_matrix(numSubcar, cpLength, csLength);
OFDMSymbolsWithRedundancy = redundancyMatrix*OFDMSymbols;
wOFDMTx = (windowTx*OFDMSymbolsWithRedundancy).';
maskedwOFDMTx = dft_rc_filt(wOFDMTx, rollOff);
end


function [filteredSignal] = dft_rc_filt(inputSignal, rollOff)
% DFT_RC_FILT   Filters a signal with an RC window using the DFT.
%

[numSymbols, signalLength] = size(inputSignal);
transformMatrix = dftmtx(2*signalLength - 1);
invertTransformMatrix = dftmtx(2*signalLength - 1)'/(2*signalLength-1);
inputSignalSpectrum = transformMatrix ...
    * ([inputSignal zeros(numSymbols, signalLength - 1)].');
windowRC = gen_raised_cosine(floor(length(inputSignalSpectrum)/2), ...
    rollOff, length(inputSignalSpectrum));
maskedSpectrum = diag(ifftshift(windowRC))*inputSignalSpectrum;
maskedSignal = (invertTransformMatrix*maskedSpectrum).';
filterTail = [zeros(1, signalLength-1); ...
    maskedSignal(1:end, signalLength+1:end)];

filteredSignal = [maskedSignal(:, 1:signalLength); ...
    zeros(1, signalLength)] + [filterTail zeros(numSymbols+1, 1)];
filteredSignal = filteredSignal(1:end-1, :);
end


function [signalRx] = tx2rx(signalTx, tailTx, cpLength, csLength, ...
    channel, snr, tailRx, prefixRemovalLength, symbolsPerTx, numSubcar)
% Function to perfrom the signal transmission and reception.
%

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


function [addRedundancyMatrix] = add_redundancy_matrix(numSubcar, ...
    cpLength, csLength)
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

identityPrefix = eye(cpLength);
zerosPrefix = zeros(cpLength, (numSubcar-cpLength));
identitySuffix = eye(csLength);
zerosSuffix = zeros(csLength, (numSubcar-csLength));
identitySubcarriers = eye(numSubcar);
addRedundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


function [windowRC] = gen_raised_cosine(windowLength, rollOff, ...
    totalLength)
% Function to generate the raised cosine window for the channel mask.
% This considers a centralized spectrum for the OFDM.
%

raisedCosineAxis = (-(rollOff+1)/2+1):1:((rollOff+1)/2-1);
raisedCosine = sin(pi/2*(.5+raisedCosineAxis/rollOff)).^2;

zeroPaddLeft = zeros(1, floor((totalLength-windowLength-2*rollOff)/2));
zeroPaddRight = zeros(1, ceil((totalLength-windowLength-2*rollOff)/2));

windowRC = [zeroPaddLeft raisedCosine ones(1, windowLength) ...
    fliplr(raisedCosine) zeroPaddRight];
end


function receivedSymbols = wofdm_rx(rxSignal, tailRx, numSubcar, ...
    prefixRemovalLength, circularShiftLength, windowRx)
% Function that performs the reception of w-OFDM symbols
%

removeRedundancyMatrix = remove_redundancy_matrix(numSubcar, tailRx, ...
    prefixRemovalLength);
overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx);
circularShiftMatrix = circular_shift_matrix(numSubcar, circularShiftLength);
transformMatrix = dftmtx(numSubcar);
receivedSymbols = transformMatrix*circularShiftMatrix ...
    * overlapAddMatrix*windowRx*removeRedundancyMatrix*rxSignal;
receivedSymbols = ifftshift(receivedSymbols, 1);  % Because of how things were done.
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


function [windowTxRC] = tx_rc_window(numSubcar, cpLength, ...
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


function [windowRx] = rx_rc_window(numSubcar, tailRx)
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


function [csLength, prefixRemovalLength, circularShiftLength] = ...
    calculate_parameters(typeOFDM, cpLength, tailTx, tailRx)
% Function to calculate parameters for each w-OFDM system
%
% - Input:
%   . typeOFDM: w-OFDM system
%   . cyclicPrefix: Number of samples in the cyclic prefix
%   . tailTx: Number of samples in rise and fall tail for Tx
%   . tailRX: Number of samples in rise and fall tail for Rx
%

switch typeOFDM
    case 'wtx'
        csLength = tailTx;
        prefixRemovalLength = cpLength;
        circularShiftLength = 0;
    case 'wrx'
        csLength = tailRx/2;
        prefixRemovalLength = cpLength - tailRx/2;
        circularShiftLength = 0;
    case 'WOLA'
        csLength = tailTx;
        prefixRemovalLength = cpLength - tailRx;
        circularShiftLength = tailRx/2;
    case 'CPW'
        csLength = tailTx + tailRx/2;
        prefixRemovalLength = cpLength - tailRx/2;
        circularShiftLength = 0;
    case 'CPwtx'
        csLength = 0;
        prefixRemovalLength = cpLength - tailTx;
        circularShiftLength = tailTx;
    case 'CPwrx'
        csLength = 0;
        prefixRemovalLength = cpLength - tailRx;
        circularShiftLength = tailRx/2;
end
end


% EoF
