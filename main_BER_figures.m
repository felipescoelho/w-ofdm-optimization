% main_BER_figures.m
%
% Luiz Felipe Coelho -- luizfelipe.coelho@smt.ufrj.br
% Out 14, 2022
%

clc
clear
close all


fprintf('Starting main_BER_figures.m ... \n\n')
fprintf('This script generates figures from the BER calculation.\n')

% Definitions
global settingsFileName resultsFolder figuresSurface
settingsFileName = 'settingsData.mat';
resultsFolder = 'ber_results';
figuresFolder = [resultsFolder '/figures'];
figuresByCPFolder = [figuresFolder '/by_cyclic_prefix'];
figuresSurface = [figuresFolder '/surface'];
if ~isfolder(figuresFolder)
    mkdir(figuresFolder)
end
if ~isfolder(figuresByCPFolder)
    mkdir(figuresByCPFolder)
end
if ~isfolder(figuresSurface)
    mkdir(figuresSurface)
end
resultFiles = dir(fullfile(resultsFolder));
typeOFDMSet = {'wtx' 'wrx' 'WOLA' 'CPW' 'CPwtx' 'CPwrx'};
settingsLoader = load(settingsFileName);
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrValues = settingsLoader.settingsData.generalSettings.snrValues;

for typeOFDMIndex = 1:length(typeOFDMSet)
    typeOFDM = typeOFDMSet{typeOFDMIndex};
    [rcWindowResults, optWindowResults] = read_results(resultFiles, ...
        typeOFDM);
    plot_ber_surf(rcWindowResults, optWindowResults, typeOFDM)
end


function plot_ber_surf(rcWindowResults, optWindowResults, typeOFDM)

global figuresSurface settingsFileName

% Definitions
width = 8.89;
height = width*2/(1+sqrt(5));
fontSize = 10;
horizontalLeftDistance = .75;
verticalBottomDistance = .85;
plotWidth = width/2 - 2*horizontalLeftDistance;
plotHeight = height - 2*verticalBottomDistance;
settingsLoader = load(settingsFileName);
cyclicPrefixVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrVector = settingsLoader.settingsData.generalSettings.snrValues;
switch typeOFDM
    case {'wtx', 'wrx', 'CPwrx', 'CPwtx'}
        fig = figure;
        fig.Name = strcat(typeOFDM, ' BER Surface');
        fig.Color = 'w';
        fig.Units = 'centimeters';
        fig.Position = [10 10 width height];
        ax1 = axes('units', 'centimeters', 'position', ...
            [horizontalLeftDistance+.5, verticalBottomDistance, ...
            plotWidth, plotHeight]);
        ax2 = axes('units', 'centimeters', 'position', ...
            [width/2 + 1.5*horizontalLeftDistance, verticalBottomDistance, ...
            plotWidth, plotHeight]);

        subplot(ax1)
        surf(cyclicPrefixVector, snrVector, ...
            10*log10(optWindowResults.'), 'LineWidth', .1, ...
            'EdgeAlpha', .91)
        set(ax1, 'Ydir', 'reverse')
        set(ax1, 'Xdir', 'reverse')
        set(ax1, 'FontSize', fontSize)
        set(ax1, 'TickLabelInterpreter', 'latex')
        set(ax1, 'linewidth', 1)
        ylim([-30 50])
        xLabelAx1 = xlabel('CP Length', 'interpreter', 'latex', ...
            'FontSize', fontSize, 'Color', 'k');
        xLabelAx1.Units = 'centimeters';
        xLabelAx1.Position(1) = xLabelAx1.Position(1) - .6;
        xLabelAx1.Position(2) = xLabelAx1.Position(2) - .1;
        xLabelAx1.Rotation = 22.5;
        yLabelAx1 = ylabel('SNR, dB', 'interpreter', 'latex', ...
            'FontSize', fontSize, 'Color', 'k');
        yLabelAx1.Units = 'centimeters';
        yLabelAx1.Position(1) = yLabelAx1.Position(1) + .5;
        yLabelAx1.Position(2) = yLabelAx1.Position(2) - .1;
        yLabelAx1.Rotation = -37.5;
        zlabel('BER, dB', 'interpreter', 'latex', 'FontSize', fontSize, ...
            'Color', 'k')
        title('(a)', 'interpreter', 'latex', 'FontSize', fontSize, ...
            'Color', 'k')
        
        subplot(ax2)
        surf(cyclicPrefixVector, snrVector, ...
            10*log10(rcWindowResults.'), 'LineWidth', .1, 'EdgeAlpha', .91)
        set(ax2, 'Ydir', 'reverse')
        set(ax2, 'Xdir', 'reverse')
        set(ax2, 'FontSize', fontSize)
        set(ax2, 'TickLabelInterpreter', 'latex')
        set(ax2, 'linewidth', 1)
        ylim([-30 50])
%         xLabelAx2 = xlabel('CP Length', 'interpreter', 'latex', ...
%             'FontSize', fontSize, 'Color', 'k');
%         xLabelAx2.Units = 'centimeters';
%         xLabelAx2.Position(1) = xLabelAx2.Position(1) - .6;
%         xLabelAx2.Position(2) = xLabelAx2.Position(2) - .1;
%         xLabelAx2.Rotation = 22.5;
%         yLabelAx2 = ylabel('SNR, dB', 'interpreter', 'latex', ...
%             'FontSize', fontSize, 'Color', 'k');
%         yLabelAx2.Units = 'centimeters';
%         yLabelAx2.Position(1) = yLabelAx2.Position(1) + .5;
%         yLabelAx2.Position(2) = yLabelAx2.Position(2) - .1;
%         yLabelAx2.Rotation = -37.5;
        title('(b)', 'interpreter', 'latex', 'FontSize', fontSize, ...
            'Color', 'k')
        
        fileName = strcat('ber_surface_', typeOFDM);
        saveas(fig, [figuresSurface '/' fileName '.eps'], 'epsc')
        saveas(fig, [figuresSurface '/' fileName '.fig'])
    case {'WOLA', 'CPW'}
        caseSet = {'1stepTx', '1stepRx', '2stepTx', '2stepRx', ...
            '3stepTx', '3stepRx'};
        for setIndex = 1:6
            fig = figure;
            fig.Name = strcat(typeOFDM, ' BER Surface ', ...
                caseSet{setIndex});
            fig.Color = 'w';
            fig.Units = 'centimeters';
            fig.Position = [10 10 width height];
            ax1 = axes('units', 'centimeters', 'position', ...
                [horizontalLeftDistance+.5, verticalBottomDistance, ...
                plotWidth, plotHeight]);
            ax2 = axes('units', 'centimeters', 'position', ...
                [width/2 + horizontalLeftDistance, ...
                verticalBottomDistance, plotWidth, plotHeight]);

            subplot(ax1)
            surf(cyclicPrefixVector, snrVector, ...
                10*log10(optWindowResults(:, :, setIndex).'), ...
                'LineWidth', .1, 'EdgeAlpha', .91)

            subplot(ax2)
            surf(cyclicPrefixVector, snrVector, ...
                10*log10(rcWindowResults.'), 'LineWidth', .1, ...
                'EdgeAlpha', .91)
            
            fileName = strcat('ber_surface_', typeOFDM, '_', ...
                caseSet{setIndex});
            saveas(fig, [figuresSurface '/' fileName '.eps'], 'epsc')
            saveas(fig, [figuresSurface '/' fileName '.fig'])
        end
end
end


function plot_ber_cut_cp(rcWindowResults, optWindowResults, cpLength, ...
    typeOFDM, folderPath)
% Funtion to plot different systems for a single CP length.
%
% - Input
%   . rcWindowResults: Matrix (CP x SNR) with results for RC window
%   . optWindowResults: Matrix (CP x SNR) with results for optimal window
%   . cpLength: Number of samples in cyclic prefix
%   . folderPath: Folder to save figure

global settingsFileName
settingsLoader = load(settingsFileName);
snrValues = settingsLoader.settingsData.generalSettings.snrValues;
width = 8.99;
height = 2*width / (1+sqrt(5));
fontSize = 11;
horizontalLeftDistance = 2.5;
verticalBottomDistance = 3;
plotWidth = width - 1.5*horizontalLeftDistance;
plotHeight = height - 1.5*verticalBottomDistance;

fig = figure;
fig.Name = strcat('BER for ', num2str(cyclicPrefix), ' CP');
fig.Units = 'centimeters';
fig.Color = 'w';
fig.Position = [2 2 width height];
ax1 = axes('units', 'centimeters', 'position', [horizontalLeftDistance, ...
    verticalBottomDistance, plotWidth, plotHeight]);
subplot(ax1)

end


function [rcWindowResults, optWindowResults] = read_results( ...
    resultFiles, typeOFDM)
global settingsFileName resultsFolder
settingsLoader = load(settingsFileName);
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrValues = settingsLoader.settingsData.generalSettings.snrValues;
filesByType = file_by_type(resultFiles, typeOFDM);
[rcFiles, optFiles] = separate_rc_opt(filesByType);
switch typeOFDM
    case {'wtx', 'wrx', 'CPwtx', 'CPwrx'}
        rcWindowResults = zeros(length(cpLengthVector), length(snrValues));
        optWindowResults = zeros(length(cpLengthVector), length(snrValues));
        for cpIndex = 1:length(cpLengthVector)
            rcDataLoader = load([resultsFolder '/' rcFiles{cpIndex}]);
            optDataLoader = load([resultsFolder '/' optFiles{cpIndex}]);
            rcWindowResults(cpIndex, :) = rcDataLoader.berRCSNR;
            optWindowResults(cpIndex, :) = optDataLoader.berSNR;
        end
    case {'WOLA', 'CPW'}
        rcWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues));
        optWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues), 6);
        for cpIndex = 1:length(cpLengthVector)
            rcDataLoader = load([resultsFolder '/' rcFiles{cpIndex}]);
            optDataLoader = load([resultsFolder '/' optFiles{cpIndex}]);
            rcWindowResults(cpIndex, :) = rcDataLoader.berRCSNR;
            optWindowResults(cpIndex, :, 1) = optDataLoader.berSNRStep1A;
            optWindowResults(cpIndex, :, 2) = optDataLoader.berSNRStep1B;
            optWindowResults(cpIndex, :, 3) = optDataLoader.berSNRStep2A;
            optWindowResults(cpIndex, :, 4) = optDataLoader.berSNRStep2B;
            optWindowResults(cpIndex, :, 5) = optDataLoader.berSNRStep3A;
            optWindowResults(cpIndex, :, 6) = optDataLoader.berSNRStep3B;
        end
end
end


function [rcFiles, optFiles] = separate_rc_opt(filesByType)
rcFiles = {};
optFiles = {};
for fileIndex = 1:length(filesByType)
    fileNameInfo = split(filesByType(fileIndex), '_');
    if isequal(fileNameInfo{1}, 'rc')
        rcFiles(end+1) = filesByType(fileIndex);  %#ok
    else
        optFiles(end+1) = filesByType(fileIndex);  %#ok
    end
end
end


function filesByType = file_by_type(resultFiles, typeOFDM)
filesByType = {};
for fileIndex = 1:length(resultFiles)
    if resultFiles(fileIndex).isdir
        continue
    end
    fileNameInfo = split(resultFiles(fileIndex).name, '_');
    fileTypeOFDM = fileNameInfo{3};
    if isequal(fileTypeOFDM, typeOFDM)
        filesByType(end+1) = {resultFiles(fileIndex).name};  %#ok
    end
end
end


function filesByCP = file_by_cp(resultFiles, cyclicPrefixLength)
filesByCP = {};
for fileIndex = 1:length(resultFiles)
    if resultFiles(fileIndex).isdir
        continue
    end
    fileNameInfo = split(resultFiles(fileIndex).name, '_');
    cpLength = fileNameInfo{end}(1:end-6);
    if isequal(cpLength, cyclicPrefixLength)
        filesByCP(end+1) = {resultFiles(fileIndex).name};  %#ok
    end
end
end


function cpList = list_cp(filesByType)
cpList = zeros(length(filesByType), 1);
for fileIndex = 1:length(filesByType)
    fileName = filesByType(fileIndex);
    fileNameInfo = split(fileName, '_');
    cpLength = str2double(fileNameInfo{end}(1:end-6));
    cpList(fileIndex) = cpLength;
end
end


% EoF

