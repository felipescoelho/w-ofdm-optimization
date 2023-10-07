% plot_opt_windows.m

clc
clear
close all


fprintf('Script to plot all windows in time.\n')
fprintf('We have a total of %d windows.\n\n', (6*2 + 4)*12)

global folderName

folderName = 'optimized_windows_new0';
assert(isdir(folderName), 'Something is wrong. Could not find %s', ...
    folderName)  %#ok
windowFiguresFolder = [folderName '/figures'];
if ~isdir(windowFiguresFolder)  %#ok
    mkdir(windowFiguresFolder)
end
settingsFileName = 'settingsData.mat';
settingsLoader = load(settingsFileName);
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
                switch typeOFDM
                    case 'wtx'
                        cpIndex = 1:str2double(cpLength);
                        csIndex = length(windowVector)-tailTx:length(windowVector);
                        tailStartIndex = 1:tailTx;
                        tailEndIndex = length(windowVector)-tailTx:length(windowVector);
                    case 'CPwtx'
                        cpIndex = 1:str2double(cpLength);
                        csIndex = [];
                        tailStartIndex = 1:tailTx;
                        tailEndIndex = length(windowVector)-tailTx:length(windowVector);
                    case 'wrx'
                        cpIndex = 1:str2double(cpLength);
                        csIndex = length(windowVector)-tailRx/2:length(windowVector);
                        tailStartIndex = 1:tailRx;
                        tailEndIndex = length(windowVector)-tailRx:length(windowVector);
                    case 'CPwrx'
                        cpIndex = 1:str2double(cpLength);
                        csIndex = [];
                        tailStartIndex = 1:tailRx;
                        tailEndIndex = length(windowVector)-tailRx:length(windowVector);
                end
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowVector, 'k'), hold on
                plot(cpIndex, windowVector(cpIndex), 'g')
                plot(csIndex, windowVector(csIndex), 'g')
                plot(tailStartIndex, windowVector(tailStartIndex), 'r')
                plot(tailEndIndex, windowVector(tailEndIndex), 'r')
                hold off
                title([cpLength ' CP'])
                xlabel('Window Sample')
                ylabel('Amplitude')
            end
            saveas(fig0, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '.fig'])
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
                fileNameInfo = split(listFiles(plotIndex).name, '_');
                cpLength = fileNameInfo{end}(1:2);
                windowLoader = load([listFiles(plotIndex).folder '/' ...
                    listFiles(plotIndex).name]);
                
                windowCaseA1 = diag(windowLoader.optimizedWindowCaseAStep1);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseA1)-tailTx-tailRx/2:length(windowCaseA1);
                else
                    csIndex = length(windowCaseA1)-tailTx:length(windowCaseA1);
                end
                tailTxEndIndex = length(windowCaseA1)-tailTx:length(windowCaseA1);
                figure(fig0)
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowCaseA1, 'k'), hold on
                plot(cpIndex, windowCaseA1(cpIndex), 'g')
                plot(csIndex, windowCaseA1(csIndex), 'g')
                plot(tailTxStartIndex, windowCaseA1(tailTxStartIndex), 'r')
                plot(tailTxEndIndex, windowCaseA1(tailTxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Sample')
                ylabel('Amplitude')
                
                windowCaseA2 = diag(windowLoader.optimizedWindowCaseAStep2);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseA2)-tailTx-tailRx/2:length(windowCaseA2);
                else
                    csIndex = length(windowCaseA2)-tailTx:length(windowCaseA2);
                end
                tailRxEndIndex = length(windowCaseA2)-tailRx:length(windowCaseA2);
                figure(fig1)
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowCaseA2, 'k'), hold on
                plot(cpIndex, windowCaseA2(cpIndex), 'g')
                plot(csIndex, windowCaseA2(csIndex), 'g')
                plot(tailRxStartIndex, windowCaseA2(tailRxStartIndex), 'r')
                plot(tailRxEndIndex, windowCaseA2(tailRxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Sample')
                ylabel('Amplitude')
                
                windowCaseA3 = diag(windowLoader.optimizedWindowCaseAStep3);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseA3)-tailTx-tailRx/2:length(windowCaseA3);
                else
                    csIndex = length(windowCaseA3)-tailTx:length(windowCaseA3);
                end
                tailTxEndIndex = length(windowCaseA3)-tailTx:length(windowCaseA3);
                figure(fig2)
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowCaseA3, 'k'), hold on
                plot(cpIndex, windowCaseA3(cpIndex), 'g')
                plot(csIndex, windowCaseA3(csIndex), 'g')
                plot(tailTxStartIndex, windowCaseA3(tailTxStartIndex), 'r')
                plot(tailTxEndIndex, windowCaseA3(tailTxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Sample')
                ylabel('Amplitude')
                
                windowCaseB1 = diag(windowLoader.optimizedWindowCaseBStep1);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseB1)-tailTx-tailRx/2:length(windowCaseB1);
                else
                    csIndex = length(windowCaseB1)-tailTx:length(windowCaseB1);
                end
                tailRxEndIndex = length(windowCaseB1)-tailRx:length(windowCaseB1);
                figure(fig3)
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowCaseB1, 'k'), hold on
                plot(cpIndex, windowCaseB1(cpIndex), 'g')
                plot(csIndex, windowCaseB1(csIndex), 'g')
                plot(tailRxStartIndex, windowCaseB1(tailRxStartIndex), 'r')
                plot(tailRxEndIndex, windowCaseB1(tailRxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Sample')
                ylabel('Amplitude')
                
                windowCaseB2 = diag(windowLoader.optimizedWindowCaseBStep2);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseB2)-tailTx-tailRx/2:length(windowCaseB2);
                else
                    csIndex = length(windowCaseB2)-tailTx:length(windowCaseB2);
                end
                tailTxEndIndex = length(windowCaseB2)-tailTx:length(windowCaseB2);
                figure(fig4)
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowCaseB2, 'k'), hold on
                plot(cpIndex, windowCaseB2(cpIndex), 'g')
                plot(csIndex, windowCaseB2(csIndex), 'g')
                plot(tailTxStartIndex, windowCaseB2(tailTxStartIndex), 'r')
                plot(tailTxEndIndex, windowCaseB2(tailTxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Sample')
                ylabel('Amplitude')
                
                windowCaseB3 = diag(windowLoader.optimizedWindowCaseBStep3);
                if isequal(typeOFDM, 'CPW')
                    csIndex = length(windowCaseB3)-tailTx-tailRx/2:length(windowCaseB3);
                else
                    csIndex = length(windowCaseB3)-tailTx:length(windowCaseB3);
                end
                tailRxEndIndex = length(windowCaseB3)-tailRx:length(windowCaseB3);
                figure(fig5)
                subplot(2, ceil(length(listFiles)/2), plotIndex)
                plot(windowCaseB3, 'k'), hold on
                plot(cpIndex, windowCaseB3(cpIndex), 'g')
                plot(csIndex, windowCaseB3(csIndex), 'g')
                plot(tailRxStartIndex, windowCaseB3(tailRxStartIndex), 'r')
                plot(tailRxEndIndex, windowCaseB3(tailRxEndIndex), 'r'), hold off
                title([cpLength ' CP'])
                xlabel('Window Samples')
                ylabel('Amplitude')
                
            end
            saveas(fig0, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '_caseA_setp1.fig'])
            saveas(fig1, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '_caseA_step2.fig'])
            saveas(fig2, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '_caseA_step3.fig'])
            saveas(fig3, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '_caseB_step1.fig'])
            saveas(fig4, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '_caseB_step2.fig'])
            saveas(fig5, [windowFiguresFolder '/optimized_windows_' ...
                typeOFDM '_caseB_step3.fig'])
    end
end



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


% EoF