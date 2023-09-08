% main_OOB_figures.m
%
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Jan 20, 2023
%


clc
clear
close all


fprintf('Starting main_OOB_figures.m ... \n\n')


% Definitions
optimizedWindowsFolder = 'optimized_windows';
oobFiguresFolder = 'oob_figures';
settingsFileName = 'settingsData.mat';
dataLoader = load(settingsFileName);
numSubcar = dataLoader.settingsData.generalSettings.numberSubcarriers;
bitsPerSubcar = dataLoader.settingsData.generalSettings.bitsPerSubcarrier;
symbolsPerTx = dataLoader.settingsData.generalSettings.symbolsPerTx;
windowFiles = dir(fullfile(optimizedWindowsFolder));
if ~isdir(oobFiguresFolder)  %#ok
    mkdir(oobFiguresFolder)
end
transmittedBits = randi([0 1], 128*bitsPerSubcar, 1);
transmittedSymbols = qammod(transmittedBits, 2^bitsPerSubcar, ...
    'InputType', 'bit', 'UnitAveragePower', true).';


for fileIndex = 1:length(windowFiles)
    if windowFiles(fileIndex).isdir
        continue
    end
    windowInfo = split(windowFiles(fileIndex).name, '_');
    typeOFDM = windowInfo{3};
    cpLength = str2double(windowInfo{end}(1:end-6));
    switch typeOFDM
        case {'wtx', 'CPwtx'}
            tailTx = dataLoader.settingsData.(typeOFDM).tailTx;
            tailRx = dataLoader.settingsData.(typeOFDM).tailRx;
            csLength = map_cs(typeOFDM, tailTx, tailRx);
            windowLoader = load(fullfile(windowFiles(fileIndex).folder, ...
                windowFiles(fileIndex).name));
            transmitterRCWindow = tx_rc_window(numSubcar, cpLength, ...
                csLength, tailTx);
            [windowedOptSpcetrum, windowedRCSpectrum, rectSpectrum] = ...
                estimate_oob(transmittedSymbols, cpLength, csLength, ...
                transmitterRCWindow, windowLoader.optimizedWindow);
            fig = figure();
            plot(20*log10(abs(rectSpectrum))), hold on
            plot(20*log10(abs(windowedRCSpectrum)))
            plot(20*log10(abs(windowedOptSpcetrum))), hold off, grid on
            legend('CP-OFDM', 'RC Window', 'Optimal Window')
            fileName = strcat('oob_estimate_', typeOFDM, '_CP_', ...
                num2str(cpLength));
            saveas(fig, [oobFiguresFolder '/' fileName '.fig'])
            saveas(fig, [oobFiguresFolder '/' fileName '.eps'], 'epsc')
        case {'CPW', 'WOLA'}
            tailTx = dataLoader.settingsData.(typeOFDM).tailTx;
            tailRx = dataLoader.settingsData.(typeOFDM).tailRx;
            csLength = map_cs(typeOFDM, tailTx, tailRx);
            windowLoader = load(fullfile(windowFiles(fileIndex).folder, ...
                windowFiles(fileIndex).name));
            transmitterRCWindow = tx_rc_window(numSubcar, cpLength, ...
                csLength, tailTx);
            % Case A - step 1 (Tx is optimized once)
            [windowedOptSpectrumCaseAStep1, ...
                windowedRCSpectrumCaseAStep1, rectSpectrumCaseAStep1] = ...
                estimate_oob(transmittedSymbols, cpLength, csLength, ...
                transmitterRCWindow, windowLoader.optimizedWindowCaseAStep1);
            fig = figure();
            plot(20*log10(abs(rectSpectrumCaseAStep1))), hold on
            plot(20*log10(abs(windowedRCSpectrumCaseAStep1)))
            plot(20*log10(abs(windowedOptSpectrumCaseAStep1)))
            hold off, grid on
            legend('CP-OFDM', 'RC Window', 'Optimal Window')
            fileName = strcat('oob_estimate_', typeOFDM, '_CP_', ...
                num2str(cpLength), '_CaseAStep1');
            saveas(fig, [oobFiguresFolder '/' fileName '.fig'])
            saveas(fig, [oobFiguresFolder '/' fileName '.eps'], 'epsc')
            % Case B - step 2 (Tx is optimized once)
            [windowedOptSpectrumCaseBStep2, ...
                windowedRCSpectrumCaseBStep2, rectSpectrumCaseBStep2] = ...
                estimate_oob(transmittedSymbols, cpLength, csLength, ...
                transmitterRCWindow, windowLoader.optimizedWindowCaseBStep2);
            fig = figure();
            plot(20*log10(abs(rectSpectrumCaseBStep2))), hold on
            plot(20*log10(abs(windowedRCSpectrumCaseBStep2)))
            plot(20*log10(abs(windowedOptSpectrumCaseBStep2)))
            hold off, grid on
            legend('CP-OFDM', 'RC Window', 'Optimal Window')
            fileName = strcat('oob_estimate_', typeOFDM, '_CP_', ...
                num2str(cpLength), '_CaseBStep2');
            saveas(fig, [oobFiguresFolder '/' fileName '.fig'])
            saveas(fig, [oobFiguresFolder '/' fileName '.eps'], 'epsc')
            % Case A - step 3 (Tx is optimized twice)
            [windowedOptSpectrumCaseAStep3, ...
                windowedRCSpectrumCaseAStep3, rectSpectrumCaseAStep3] = ...
                estimate_oob(transmittedSymbols, cpLength, csLength, ...
                transmitterRCWindow, windowLoader.optimizedWindowCaseAStep3);
            fig = figure();
            plot(20*log10(abs(rectSpectrumCaseAStep3))), hold on
            plot(20*log10(abs(windowedRCSpectrumCaseAStep3)))
            plot(20*log10(abs(windowedOptSpectrumCaseAStep3)))
            hold off, grid on
            legend('CP-OFDM', 'RC Window', 'Optimal Window')
            fileName = strcat('oob_estimate_', typeOFDM, '_CP_', ...
                num2str(cpLength), '_CaseAStep3');
            saveas(fig, [oobFiguresFolder '/' fileName '.fig'])
            saveas(fig, [oobFiguresFolder '/' fileName '.eps'], 'epsc')
        otherwise
            continue
    end
    close all
end


function [windowedOptSpectrum, windowedRCSpectrum, rectSpectrum] = ...
    estimate_oob(transmittedSymbols, cpLength, csLength, ...
    RCWindowTx, optWindowTx)
% Function to calculate OOB (out of band radiation).
%
% - Input
%   . transmittedSymbols: QAM symbols to be transmitted using OFDM system.
%   . numSubcar: Number of subcarriers in OFDM symbol
%   . cpLength: Number of samples in cyclic prefix
%   . csLength: Number of samples in cyclic suffix
%   . windowTx: Window used in transmitter

offsetLength = (256 - 128)/2;
txOFDMSig = ofdm_tx(transmittedSymbols, cpLength, csLength, offsetLength);
txRCWindowOFDMSig = txOFDMSig .* diag(RCWindowTx);
txOptWindowOFDMSig = txOFDMSig .* diag(optWindowTx);
rectSpectrum = fftshift(fft(txOFDMSig, 2^11));
windowedOptSpectrum = fftshift(fft(txOptWindowOFDMSig, 2^11));
windowedRCSpectrum = fftshift(fft(txRCWindowOFDMSig, 2^11));
end


function txSigOFDM = ofdm_tx(transmittedSymbols, ...
    cpLength, csLength, offsetLength)
% Function to prepare the OFDM symbols to perform the windowing process.
%
% - Input:
%   . transmittedSymbols: Symbols from digital modulation
%   . numSubcar: Number of subcarriers in OFDM symbol
%   . cpLength: Number of samples in cyclic prefix
%   . csLength: Number of samples in cyclic suffix
%

symbolsInOFDM = [zeros(offsetLength, 1); transmittedSymbols.'; ...
    zeros(offsetLength, 1)];
ifftOut = ifft(ifftshift(symbolsInOFDM));
txSigOFDM = [ifftOut(end-cpLength+1:end); ifftOut; ...
    ifftOut(1:csLength)];
end


function windowedOFDM = wofdm_tx(transmittedSymbols, numSubcar, ...
    cpLength, csLength, windowTx)
% Function that prepares the symbols from the digital modulation to be
% transmitted as OFDM symbols.
%
% - Input:
%   . transmittedSymbols: Symbols from digital modulation
%   . numberSubcarriers: Number of subcarriers in OFDM symbol
%   . cyclicPrefix: Number of elements in cyclic prefix
%   . cyclicSuffix: Number of elements in cyclic suffix
%   . windowTx: Transmitter window

invertTransformMatrix = dftmtx(numSubcar)'/numSubcar;
ofdmSymbols = invertTransformMatrix*transmittedSymbols.';
redundancyMatrix = add_redundancy_matrix(numSubcar, ...
    cpLength, csLength);
ofdmSymbolsWithRedundancy = redundancyMatrix*ofdmSymbols;
windowedOFDM = (windowTx*ofdmSymbolsWithRedundancy).';
end


function windowTx = tx_rc_window(numSubcar, cpLength, ...
    csLength, tailTx)
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
onesTransmitted = ones(1, numSubcar+cpLength+csLength ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
windowTx = diag(windowVector);
end


function addRedundancyMatrix = add_redundancy_matrix(numSubcar, ...
    cpLength, csLength)
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

identityPrefix = eye(cpLength);
zerosPrefix = zeros(cpLength, (numSubcar-cpLength));
identitySuffix = eye(csLength);
zerosSuffix = zeros(csLength, (numSubcar-csLength));
identitySubcarriers = eye(numSubcar);
addRedundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


function paddedSignal = pad_spectrum(signalInput, newLength)
    % Function to perform the zero padding. Expects a line-vector.
    padLength = newLength - length(signalInput);
    assert(padLength > 0, ['The signal can not be padded because it is' ...
        ' longer than expected.'])
    
    paddedSignal = [zeros(1, ceil(padLength/2)), signalInput, ...
        zeros(1, floor(padLength/2))];
    
end


function csLength = map_cs(typeOFDM, tailTx, tailRx)
% Function to calculate parameters for each w-OFDM system
%
% - Input:
%   . typeOFDM: w-OFDM system
%   . cyclicPrefix: Number of samples in the cyclic prefix
%   . tailTx: Number of samples in rise and fall tail for Tx
%   . tailRX: Number of samples in rise and fall tail for Rx

switch typeOFDM
    case 'wtx'
        csLength = tailTx;
    case 'WOLA'
        csLength = tailTx;
    case 'CPW'
        csLength = tailTx + tailRx/2;
    case 'CPwtx'
        csLength = 0;
    otherwise
        csLength = 0;
end
end


% EoF
