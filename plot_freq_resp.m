% plot_freq_resp.m

clc
clear
close all


fprintf('Script to compare the frequency response of optimized windows')
fprintf(' and RC windows of the same size for the same system.\n')
fprintf('We have a total of %d windows.\n\n', (6*2+4)*12)

global folderName

folderName = 'optimized_windows_new';
assert(isdir(folderName), 'Something is wrong. Could not find %s', ...
    folderName)  %#ok
windowFiguresFolder = [folderName '/figures'];
if ~isdir(windowFiguresFolder)  %#ok
    mkdir(windowFiguresFolder)
end
resultsFolder = [windowFiguresFolder '/freqeuncy_compare'];
settingsFileName = 'settingsData.mat';
settingsLoader = load(settingsFileName);
numSubcar = settingsLoader.settingsData.generalSettings.numberSubcarriers;
typeOFDMList = {'wtx', 'CPwtx', 'wrx', 'CPwrx', 'CPW', 'WOLA'};


for idx = 1:length(typeOFDMList)
    typeOFDM = typeOFDMList{idx};
    listFiles = filter_by_typeOFDM(typeOFDM);
    tailTx = settingsLoader.settingsData.(typeOFDM).tailTx;
    tailRx = settingsLoader.settingsData.(typeOFDM).tailRx;
    switch typeOFDM
        case {'wtx', 'CPwtx', 'wrx', 'CPwrx'}
            fig0 = figure;
            fig0.Name = ['Windows for ' typeOFDM '-OFDM systems'];
            for plotIndex = 1:length(listFiles)
                fileNameInfo = split(listFiles(plotIndex).name, '_');
                cpLength = fileNameInfo{end}(1:2);
                windowLoader = load([listFiles(plotIndex).folder '/' ...
                    listFiles(plotIndex).name]);
                windowVector = diag(windowLoader.optimizedWindow);
                cpIndex = 1:str2double(cpLength);
                switch typeOFDM
                    case 'wtx'
                        csLength = tailTx;
                        csIndex = length(windowVector)-tailTx:length(windowVector);
                        tailStartIndex = 1:tailTx;
                        tailEndIndex = length(windowVector)-tailTx:length(windowVector);
                        windowRC = transmitter_rc_window(numSubcar, ...
                            str2double(cpLength), csLength, tailTx);
                    case 'CPwtx'
                        csLength = 0;
                        csIndex = [];
                        tailStartIndex = 1:tailTx;
                        tailEndIndex = length(windowVector)-tailTx:length(windowVector);
                        windowRC = transmitter_rc_window(numSubcar, ...
                            str2double(cpLength), csLength, tailTx);
                    case 'wrx'
                        csIndex = length(windowVector)-tailRx/2:length(windowVector);
                        tailStartIndex = 1:tailRx;
                        tailEndIndex = length(windowVector)-tailRx:length(windowVector);
                        windowRC = receiver_rc_window(numSubcar, tailRx);
                    case 'CPwrx'
                        csIndex = [];
                        tailStartIndex = 1:tailRx;
                        tailEndIndex = length(windowVector)-tailRx:length(windowVector);
                        windowRC = receiver_rc_window(numSubcar, tailRx);
                end
                [winOptFreq, freqIndex] = freqz(windowVector, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                % Plot in time
                if plotIndex > ceil(length(listFiles)/2)
                    plotCorrection = plotIndex + ceil(length(listFiles)/2);
                else
                    plotCorrection = plotIndex;
                end
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowVector, 'k'), hold on
                plot(cpIndex, windowVector(cpIndex), 'g')
                plot(csIndex, windowVector(csIndex), 'g')
                plot(tailStartIndex, windowVector(tailStartIndex), 'r')
                plot(tailEndIndex, windowVector(tailEndIndex), 'r')
                hold off
                title([cpLength ' CP - Time domain'])
                xlabel('Window Sample')
                ylabel('Amplitude')
                % Plot in frequency
                subplot(4, ceil(length(listFiles)/2), ...
                    plotCorrection+ceil(length(listFiles)/2))
                plot(freqIndex, 20*log10(abs(winRCFreq))), hold on
                plot(freqIndex, 20*log10(abs(winOptFreq))), hold off
                legend('RC', 'Opt')
                title([cpLength ' CP - Freq domain'])
                xlabel('Normalized Frequency')
                ylabel('Magnitude')
                ylim([-100 100])
            end
    end
end


% --- Functions ---


function filteredList = filter_by_typeOFDM(typeOFDM)
% A function to filter files by system type.
%
global folderName

resultFiles = dir(fullfile(folderName));

dataList = [];
for idx = 1:length(resultFiles)
    if ~endsWith(resultFiles(idx).name, '.mat')
        continue
    end
    fileNameInfo = split(resultFiles(idx).name, '_');
    if isequal(fileNameInfo{3}, typeOFDM)
        dataList = [dataList resultFiles(idx)];  %#ok
    end
end
filteredList = struct(dataList);
end


function windowTxRC = transmitter_rc_window(numSubcar, cpLength, ...
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

raisedCosineAxis = (-(tailTx+1)/2 + 1):1:((tailTx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailTx)).^2;
onesTransmitted = ones(1, numSubcar+cpLength+csLength ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
windowTxRC = windowVector;
end


function windowRxRC = receiver_rc_window(numSubcar, tailRx)
% Function to generate the Raised Cosine receiver window.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail for the receiver
%   window.

raisedCosineAxis = (-(tailRx+1)/2+1):1:((tailRx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailRx)).^2;
onesReceived = ones(1, numSubcar-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
windowRxRC = windowVector;
end

% EoF
