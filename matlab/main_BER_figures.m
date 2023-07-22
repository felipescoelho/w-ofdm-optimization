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
global settingsFileName resultsFolder figuresSurface listSystems ...
    figuresByCPFolder
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
listSystems = {'wtx' 'wrx' 'CPwtx' 'CPwrx' 'WOLA' 'CPW'};
settingsLoader = load(settingsFileName);
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrValues = settingsLoader.settingsData.generalSettings.snrValues;

% Plot all surfaces
for typeOFDMIndex = 1:length(listSystems)
    typeOFDM = listSystems{typeOFDMIndex};
    [rcWindowResults, optWindowResults] = read_results(resultFiles, ...
        typeOFDM);
    plot_ber_surf(rcWindowResults, optWindowResults, typeOFDM)
end
close all
% Plot cut
plot_ber_cut_cp(18)

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
            zlabel('BER, dB', 'interpreter', 'latex', 'FontSize', ...
                fontSize, 'Color', 'k')
            title('(a)', 'interpreter', 'latex', 'FontSize', fontSize, ...
                'Color', 'k')

            subplot(ax2)
            surf(cyclicPrefixVector, snrVector, ...
                10*log10(rcWindowResults.'), 'LineWidth', .1, ...
                'EdgeAlpha', .91)
            set(ax2, 'Ydir', 'reverse')
            set(ax2, 'Xdir', 'reverse')
            set(ax2, 'FontSize', fontSize)
            set(ax2, 'TickLabelInterpreter', 'latex')
            set(ax2, 'linewidth', 1)
            ylim([-30 50])
            title('(b)', 'interpreter', 'latex', 'FontSize', fontSize, ...
            'Color', 'k')
            
            fileName = strcat('ber_surface_', typeOFDM, '_', ...
                caseSet{setIndex});
            saveas(fig, [figuresSurface '/' fileName '.eps'], 'epsc')
            saveas(fig, [figuresSurface '/' fileName '.fig'])
        end
end
end


function plot_ber_cut_cp(cpLength)
% Funtion to plot different systems for a single CP length.
%
% - Input
%   . rcWindowResults: Matrix (CP x SNR) with results for RC window
%   . optWindowResults: Matrix (CP x SNR) with results for optimal window
%   . cpLength: Number of samples in cyclic prefix
%   . folderPath: Folder to save figure

global settingsFileName listSystems figuresByCPFolder resultsFolder

settingsLoader = load(settingsFileName);
snrValues = settingsLoader.settingsData.generalSettings.snrValues;
resultFiles = dir(fullfile(resultsFolder));
selectedResults = zeros(length(listSystems), length(snrValues), 2);

for idx = 1:length(resultFiles)
    if resultFiles(idx).isdir
        continue
    end
    fileNameInfo = split(resultFiles(idx).name, '_');
    typeOFDM = fileNameInfo{3};
    typeIndex = find(strcmp(listSystems, typeOFDM));
    switch typeOFDM
        case {'wtx', 'wrx', 'CPwrx', 'CPwtx'}
            if isequal(fileNameInfo{end}(1:2), num2str(cpLength))
                if isequal(fileNameInfo{1}, 'rc')
                    selectedResults(typeIndex, :, 2) = load(strcat(...
                        resultFiles(idx).folder, '/', ...
                        resultFiles(idx).name)).berRCSNR;
                else
                    selectedResults(typeIndex, :, 1) = load(strcat(...
                        resultFiles(idx).folder, '/', ...
                        resultFiles(idx).name)).berSNR;
                end
            end
        case {'WOLA', 'CPW'}
            if isequal(fileNameInfo{end}(1:2), num2str(cpLength))
                if isequal(fileNameInfo{1}, 'rc')
                    selectedResults(typeIndex, :, 2) = load(strcat(...
                        resultFiles(idx).folder, '/', ...
                        resultFiles(idx).name)).berRCSNR;
                else
                    selectedResults(typeIndex, :, 1) = load(strcat(...
                        resultFiles(idx).folder, '/', ...
                        resultFiles(idx).name)).berSNRStep3B;
                end
            end
        otherwise
            fprintf('No such implementation...\n')
    end
end

% Figure settings
width = 22;
height = 2*width / (1+sqrt(5));
fontSize = 14;
lineWidth = 1.5;
markerSize = 13;
horizontalLeftDistance = 2.5;
verticalBottomDistance = 2;
plotWidth = width - 1.5*horizontalLeftDistance;
plotHeight = height - 1.5*verticalBottomDistance;

fig = figure;
fig.Name = strcat('BER for ', num2str(cpLength), ' CP');
fig.Units = 'centimeters';
fig.Color = 'w';
fig.Position = [10 10 width height];
ax1 = axes('units', 'centimeters', 'position', [horizontalLeftDistance, ...
    verticalBottomDistance, plotWidth, plotHeight]);
subplot(ax1)
indexedLineStyle = {'d-', 's-', 'o-', '*-', '+-', 'x-'};
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 1), ...
        indexedLineStyle{idx}, 'LineWidth', lineWidth, 'MarkerSize', ...
        markerSize)
    if idx == 1
        hold on
    end
end
ax1.ColorOrderIndex = 1;
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 2), ...
        strcat(indexedLineStyle{idx}, '-'), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize)
end
% Create zoom box
topLeft = [40 .5*1e-3];
bottomRight = [50 .5^2*1e-4];
draw_box(topLeft, bottomRight, lineWidth)
hold off, grid on
set(ax1, 'FontSize', fontSize)
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'linewidth', lineWidth)
set(ax1, 'XColor', 'k')
set(ax1, 'YColor', 'k')
xlabel('SNR, dB', 'interpreter', 'latex', 'FontSize', fontSize)
ylabel('BER', 'interpreter', 'latex', 'FontSize', fontSize)
lgd = legend(listSystems);
lgd.Units = 'centimeters';
lgd.Interpreter = 'latex';
lgd.FontSize = fontSize;
lgd.Position(1) = horizontalLeftDistance + .25;
lgd.Position(2) = verticalBottomDistance + .25;
lgd.NumColumns = 1;

% Zoom plot definitions:
lineWidthZoom = 1.5;
fontSizeZoom = 14;
plotWidthZoom = plotWidth/3;
plotHeightZoom = plotHeight/3;
horizontalLeftDistanceZoom = horizontalLeftDistance + plotWidth/2 ...
    - plotWidthZoom/2;
verticalBottomDistanceZoom = verticalBottomDistance + 1;
ax2 = axes('units', 'centimeters', 'position', ...
    [horizontalLeftDistanceZoom, verticalBottomDistanceZoom, ...
    plotWidthZoom, plotHeightZoom]);
subplot(ax2)
indexedLineStyle = {'d-', 's-', 'o-', '*-', '+-', 'x-'};
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 1), ...
        indexedLineStyle{idx}, 'LineWidth', lineWidthZoom, ...
        'MarkerSize', markerSize)
    if idx == 1
        hold on
    end
end
ax2.ColorOrderIndex = 1;
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 2), ...
        strcat(indexedLineStyle{idx}, '-'), 'LineWidth', lineWidthZoom, ...
        'MarkerSize', markerSize)
end
hold off, grid on
xlim([topLeft(1) bottomRight(1)])
ylim([bottomRight(2) topLeft(2)])
set(ax2, 'FontSize', fontSizeZoom)
set(ax2, 'TickLabelInterpreter', 'latex')
set(ax2, 'linewidth', lineWidthZoom)
set(ax2, 'XColor', 'k')
set(ax2, 'YColor', 'k')

fileName = strcat('ber_cut_', num2str(cpLength), '_CP.eps');
filePath = [figuresByCPFolder '/' fileName];
saveas(fig, filePath, 'epsc')
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


function draw_box(topLeft, bottomRight, lineWidth)

top = topLeft(2);
left = topLeft(1);
bottom = bottomRight(2);
right = bottomRight(1);


plot([left right], [top top], 'k', 'linewidth', lineWidth,...
    'Handlevisibility', 'off'), hold on
plot([right right], [top bottom], 'k', 'linewidth', lineWidth,...
    'HandleVisibility', 'off')
plot([left right], [bottom bottom], 'k', 'linewidth', lineWidth,...
    'HandleVisibility', 'off')
plot([left left], [bottom top], 'k', 'linewidth', lineWidth,...
    'HandleVisibility', 'off')
end




% EoF

