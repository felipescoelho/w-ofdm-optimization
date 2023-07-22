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
resultsFolder = [windowFiguresFolder '/frequency_compare'];
if ~isdir(resultsFolder)  %#ok
    mkdir(resultsFolder)
end
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
            case {'CPW', 'WOLA'}
            fig0 = figure;
            fig0.Name = ['Windows for ' typeOFDM ...
                '-OFDM systems Case A - 1st Step'];
            fig1 = figure;
            fig1.Name = ['Windows for ' typeOFDM ...
                '-OFDM systems Case A - 2nd Step'];
            fig2 = figure;
            fig2.Name = ['Windows for ' typeOFDM ...
                '-OFDM systems Case A - 3rd step'];
            fig3 = figure;
            fig3.Name = ['Windows for ' typeOFDM ...
                '-OFDM systems Case B - 1st Step'];
            fig4 = figure;
            fig4.Name = ['Windows for ' typeOFDM ...
                '-OFDM systems Case B - 2nd Step'];
            fig5 = figure;
            fig5.Name = ['Windows for ' typeOFDM ...
                '-OFDM systems Case B - 3rd Step'];
            cpIndex = 1:str2double(cpLength);
            tailTxStartIndex = 1:tailTx;
            tailRxStartIndex = 1:tailRx;
            for plotIndex = 1:length(listFiles)
                % Correction to separate time from frequency when plotting
                if plotIndex > ceil(length(listFiles)/2)
                    plotCorrection = plotIndex + ceil(length(listFiles)/2);
                else
                    plotCorrection = plotIndex;
                end
                fileNameInfo = split(listFiles(plotIndex).name, '_');
                cpLength = fileNameInfo{end}(1:2);
                windowLoader = load([listFiles(plotIndex).folder '/' ...
                    listFiles(plotIndex).name]);
                
                % Case A Step 1 -- Windowing at transmitter
                windowCaseA1 = diag(windowLoader.optimizedWindowCaseAStep1);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseA1)-tailTx-tailRx/2:length(windowCaseA1);
                else
                    csIndex = length(windowCaseA1)-tailTx:length(windowCaseA1);
                end
                tailTxEndIndex = length(windowCaseA1)-tailTx:length(windowCaseA1);
                windowRC = transmitter_rc_window(numSubcar, ...
                            str2double(cpLength), csLength, tailTx);
                [winOptFreq, freqIndex] = freqz(windowCaseA1, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                
                figure(fig0)
                % Plot in time
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowCaseA1, 'k'), hold on
                plot(cpIndex, windowCaseA1(cpIndex), 'g')
                plot(csIndex, windowCaseA1(csIndex), 'g')
                plot(tailTxStartIndex, windowCaseA1(tailTxStartIndex), 'r')
                plot(tailTxEndIndex, windowCaseA1(tailTxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
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
                % Case A Step 2 -- Windowing at receiver
                windowCaseA2 = diag(windowLoader.optimizedWindowCaseAStep2);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseA2)-tailTx-tailRx/2:length(windowCaseA2);
                else
                    csIndex = length(windowCaseA2)-tailTx:length(windowCaseA2);
                end
                tailRxEndIndex = length(windowCaseA2)-tailRx:length(windowCaseA2);
                windowRC = receiver_rc_window(numSubcar, tailRx);
                [winOptFreq, freqIndex] = freqz(windowCaseA2, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                
                figure(fig1)
                % Plot in time
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowCaseA2, 'k'), hold on
                plot(cpIndex, windowCaseA2(cpIndex), 'g')
                plot(csIndex, windowCaseA2(csIndex), 'g')
                plot(tailRxStartIndex, windowCaseA2(tailRxStartIndex), 'r')
                plot(tailRxEndIndex, windowCaseA2(tailRxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
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
                % Case A Step 3 -- Windowing at transmitter
                windowCaseA3 = diag(windowLoader.optimizedWindowCaseAStep3);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseA3)-tailTx-tailRx/2:length(windowCaseA3);
                else
                    csIndex = length(windowCaseA3)-tailTx:length(windowCaseA3);
                end
                tailTxEndIndex = length(windowCaseA3)-tailTx:length(windowCaseA3);
                windowRC = transmitter_rc_window(numSubcar, ...
                            str2double(cpLength), csLength, tailTx);
                [winOptFreq, freqIndex] = freqz(windowCaseA3, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                
                figure(fig2)
                % Plot in time
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowCaseA3, 'k'), hold on
                plot(cpIndex, windowCaseA3(cpIndex), 'g')
                plot(csIndex, windowCaseA3(csIndex), 'g')
                plot(tailTxStartIndex, windowCaseA3(tailTxStartIndex), 'r')
                plot(tailTxEndIndex, windowCaseA3(tailTxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
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
                % Case B Step 1 -- Windowing at receiver
                windowCaseB1 = diag(windowLoader.optimizedWindowCaseBStep1);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseB1)-tailTx-tailRx/2:length(windowCaseB1);
                else
                    csIndex = length(windowCaseB1)-tailTx:length(windowCaseB1);
                end
                tailRxEndIndex = length(windowCaseB1)-tailRx:length(windowCaseB1);
                windowRC = receiver_rc_window(numSubcar, tailRx);
                [winOptFreq, freqIndex] = freqz(windowCaseB1, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                
                figure(fig3)
                % Plot in time
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowCaseB1, 'k'), hold on
                plot(cpIndex, windowCaseB1(cpIndex), 'g')
                plot(csIndex, windowCaseB1(csIndex), 'g')
                plot(tailRxStartIndex, windowCaseB1(tailRxStartIndex), 'r')
                plot(tailRxEndIndex, windowCaseB1(tailRxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
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
                % Case B Step 2 -- Windowing at transmitter
                windowCaseB2 = diag(windowLoader.optimizedWindowCaseBStep2);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseB2)-tailTx-tailRx/2:length(windowCaseB2);
                else
                    csIndex = length(windowCaseB2)-tailTx:length(windowCaseB2);
                end
                tailTxEndIndex = length(windowCaseB2)-tailTx:length(windowCaseB2);
                windowRC = transmitter_rc_window(numSubcar, ...
                            str2double(cpLength), csLength, tailTx);
                [winOptFreq, freqIndex] = freqz(windowCaseB2, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                
                figure(fig4)
                % Plot in time
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowCaseB2, 'k'), hold on
                plot(cpIndex, windowCaseB2(cpIndex), 'g')
                plot(csIndex, windowCaseB2(csIndex), 'g')
                plot(tailTxStartIndex, windowCaseB2(tailTxStartIndex), 'r')
                plot(tailTxEndIndex, windowCaseB2(tailTxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
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
                % Case B Step 3 -- Windowing at receiver
                windowCaseB3 = diag(windowLoader.optimizedWindowCaseBStep3);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseB3)-tailTx-tailRx/2:length(windowCaseB3);
                else
                    csIndex = length(windowCaseB3)-tailTx:length(windowCaseB3);
                end
                tailRxEndIndex = length(windowCaseB3)-tailRx:length(windowCaseB3);
                windowRC = receiver_rc_window(numSubcar, tailRx);
                [winOptFreq, freqIndex] = freqz(windowCaseB3, 1, 1024);
                [winRCFreq, ~] = freqz(windowRC, 1, 1024);
                
                figure(fig5)
                % Plot in time
                subplot(4, ceil(length(listFiles)/2), plotCorrection)
                plot(windowCaseB3, 'k'), hold on
                plot(cpIndex, windowCaseB3(cpIndex), 'g')
                plot(csIndex, windowCaseB3(csIndex), 'g')
                plot(tailRxStartIndex, windowCaseB3(tailRxStartIndex), 'r')
                plot(tailRxEndIndex, windowCaseB3(tailRxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Samples')
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
            saveas(fig0, [resultsFolder '/optimized_windows_' ...
                typeOFDM '_caseA_setp1.fig'])
            saveas(fig1, [resultsFolder '/optimized_windows_' ...
                typeOFDM '_caseA_step2.fig'])
            saveas(fig2, [resultsFolder '/optimized_windows_' ...
                typeOFDM '_caseA_step3.fig'])
            saveas(fig3, [resultsFolder '/optimized_windows_' ...
                typeOFDM '_caseB_step1.fig'])
            saveas(fig4, [resultsFolder '/optimized_windows_' ...
                typeOFDM '_caseB_step2.fig'])
            saveas(fig5, [resultsFolder '/optimized_windows_' ...
                typeOFDM '_caseB_step3.fig'])
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
