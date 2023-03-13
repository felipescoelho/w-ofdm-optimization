% ber_figures_with_channel_mask.m
%
% Luiz Felipe da S. Coelho -- luizfelipe.coelho@smt.ufrj.br
% Mar 10, 2023
%


clc
clear
close all


fprintf('Starting ber_figures_with_channel_mask.m ... \n\n')


% Definitions
global settingsFileName resultsFolder figuresFolder listSystems
settingsFileName = 'settingsData.mat';
resultsFolder = 'ber_results/simulation_with_channel_mask';
assert(isdir(resultsFolder), 'Run main_channel_mask.m.\n')  %#ok
figuresFolder = 'ber_results/figures/simulation_with_channel_mask';
if ~isdir(figuresFolder)  %#ok
    mkdir(figuresFolder)
end
resultsFiles = dir(fullfile(resultsFolder));
listSystems = {'wtx', 'wrx', 'CPwtx', 'CPwrx', 'WOLA', 'CPW'};
% for typeOFDMIndex = 1:length(listSystems)
%     [rcWindowResults, rcMaskedWindowResults, optWindowResults, ...
%         optMaskedWindowResults] = read_by_type(listSystems{typeOFDMIndex});
%     plot_ber_surf(rcWindowResults, rcMaskedWindowResults, ...
%         optWindowResults, optMaskedWindowResults, ...
%         listSystems{typeOFDMIndex})
% end
plot_ber_cut_cp(22)


function plot_ber_cut_cp(cpLength)
% PLOT_BER_CUT_CP   Plot different systems for a single CP length.
%

global settingsFileName resultsFolder listSystems figuresFolder

settingsLoader = load(settingsFileName);
snrValues = settingsLoader.settingsData.generalSettings.snrValues;
resultsFiles = dir(fullfile(resultsFolder));
selectedResults = zeros(length(listSystems), length(snrValues), 4);
for fileIndex = 1:length(resultsFiles)
    if resultsFiles(fileIndex).isdir
        continue
    end
    fileNameInfo = split(resultsFiles(fileIndex).name, '_');
    typeOFDM = fileNameInfo{end-1};
    typeIndex = strcmp(listSystems, typeOFDM);
    if isequal(fileNameInfo{end}(1:2), num2str(cpLength))
        if isequal(fileNameInfo{1}, 'rc')
            selectedResults(typeIndex, :, 2) = load(strcat( ...
                resultsFiles(fileIndex).folder, '/', ...
                resultsFiles(fileIndex).name)).berRCSNR;
        elseif isequal(fileNameInfo{1}, 'optimized')
            switch typeOFDM
                case {'wtx', 'wrx', 'CPwtx', 'CPwrx'}
                    selectedResults(typeIndex, :, 1) = load(strcat( ...
                        resultsFiles(fileIndex).folder, '/', ...
                        resultsFiles(fileIndex).name)).berSNR;
                case {'WOLA', 'CPW'}
                    selectedResults(typeIndex, :, 1) = load(strcat( ...
                        resultsFiles(fileIndex).folder, '/', ...
                        resultsFiles(fileIndex).name)).berSNRStep3B;
            end
        else
            if isequal(fileNameInfo{2}, 'rc')
                selectedResults(typeIndex, :, 4) = load(strcat( ...
                    resultsFiles(fileIndex).folder, '/', ...
                    resultsFiles(fileIndex).name)).berMaskedRCSNR;
            else
                switch typeOFDM
                    case {'wtx', 'wrx', 'CPwtx', 'CPwrx'}
                        selectedResults(typeIndex, :, 3) = load(strcat( ...
                            resultsFiles(fileIndex).folder, '/', ...
                            resultsFiles(fileIndex).name)).berMaskedSNR;
                    case {'WOLA', 'CPW'}
                        selectedResults(typeIndex, :, 3) = load(strcat( ...
                            resultsFiles(fileIndex).folder, '/', ...
                            resultsFiles(fileIndex).name)).berMaskedSNRStep3B;
                end
            end
        end
    end
end
width = 22;
height = 2*width/(1+sqrt(5));
fontSize = 14;
lineWidth = 1.5;
markerSize = 13;
horizontalLeftDistance = 2;
verticalBottomDistance = 2;
plotWidth = width/2 - horizontalLeftDistance;
plotHeight = height - 1.5*verticalBottomDistance;
lineWidthZoom = 1.5;
fontSizeZoom = 14;
plotWidthZoom = plotWidth/2 - 1.5;
plotHeightZoom = plotHeight/2;
horizontalLeftDistanceZoom = horizontalLeftDistance + 1.25;
verticalBottomDistanceZoom = verticalBottomDistance + 1;
fig = figure;
fig.Name = strcat('BER for ', num2str(cpLength), ' CP');
fig.Units = 'centimeters';
fig.Color = 'w';
fig.Position = [10 10 width height];
ax1 = axes('units', 'centimeters', 'position', [horizontalLeftDistance, ...
    verticalBottomDistance, plotWidth, plotHeight]);
ax1Zoom = axes('units', 'centimeters', 'position', ...
    [horizontalLeftDistanceZoom, verticalBottomDistanceZoom, ...
    plotWidthZoom, plotHeightZoom]);
ax2 = axes('units', 'centimeters', 'position', ...
    [1.75*horizontalLeftDistance + plotWidth, verticalBottomDistance, ...
    plotWidth, plotHeight]);
ax2Zoom = axes('units', 'centimeters', 'position', ...
    [horizontalLeftDistanceZoom+plotWidth+.75*horizontalLeftDistance, ...
    verticalBottomDistanceZoom, plotWidthZoom, plotHeightZoom]);
% Plot 1
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
topLeft = [40 2*1e-4];
bottomRight = [50 1e-5];
draw_box(topLeft, bottomRight, lineWidth)
hold off, grid on
xlim([-20 50])
set(ax1, 'FontSize', fontSize)
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'linewidth', lineWidth)
set(ax1, 'XColor', 'k')
set(ax1, 'YColor', 'k')
xlabel('SNR, dB', 'interpreter', 'latex', 'FontSize', fontSize, ...
    'Color', 'k')
ylabel('BER', 'interpreter', 'latex', 'FontSize', fontSize, 'Color', 'k')
title('(a)', 'interpreter', 'latex', 'FontSize', fontSize, 'Color', 'k')
% Zoom plot 1
subplot(ax1Zoom)
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 1), ...
        indexedLineStyle{idx}, 'LineWidth', lineWidthZoom, 'MarkerSize', ...
        markerSize)
    if idx == 1
        hold on
    end
end
ax1Zoom.ColorOrderIndex = 1;
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 2), ...
        strcat(indexedLineStyle{idx}, '-'), 'LineWidth', lineWidthZoom, ...
        'MarkerSize', markerSize)
end
hold off, grid on
xlim([topLeft(1) bottomRight(1)])
ylim([bottomRight(2) topLeft(2)])
set(ax1Zoom, 'FontSize', fontSizeZoom)
set(ax1Zoom, 'TickLabelInterpreter', 'latex')
set(ax1Zoom, 'linewidth', lineWidthZoom)
set(ax1Zoom, 'XColor', 'k')
set(ax1Zoom, 'YColor', 'k')
% Plot 2
subplot(ax2)
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 3), ...
        indexedLineStyle{idx}, 'LineWidth', lineWidth, 'MarkerSize', ...
        markerSize)
    if idx == 1
        hold on
    end
end
ax2.ColorOrderIndex = 1;
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 4), ...
        strcat(indexedLineStyle{idx}, '-'), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize)
end
topLeft = [40 3*1e-4];
bottomRight = [50 2*1e-5];
draw_box(topLeft, bottomRight, lineWidth)
lgd = legend(listSystems);
lgd.Units = 'centimeters';
lgd.Interpreter = 'latex';
lgd.FontSize = fontSize-2;
lgd.NumColumns = 1;
hold off, grid on
xlim([-20 50])
set(ax2, 'FontSize', fontSize)
set(ax2, 'TickLabelInterpreter', 'latex')
set(ax2, 'linewidth', lineWidth)
set(ax2, 'XColor', 'k')
set(ax2, 'YColor', 'k')
xlabel('SNR, dB', 'interpreter', 'latex', 'FontSize', fontSize, ...
    'Color', 'k')
title('(b)', 'interpreter', 'latex', 'FontSize', fontSize, 'Color', 'k')
% Zoom plot 2
subplot(ax2Zoom)
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 3), ...
        indexedLineStyle{idx}, 'LineWidth', lineWidth, 'MarkerSize', ...
        markerSize)
    if idx == 1
        hold on
    end
end
ax2Zoom.ColorOrderIndex = 1;
for idx = 1:length(listSystems)
    semilogy(snrValues, selectedResults(idx, :, 4), ...
        strcat(indexedLineStyle{idx}, '-'), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize)
end
hold off, grid on
xlim([topLeft(1) bottomRight(1)])
ylim([bottomRight(2) topLeft(2)])
set(ax2Zoom, 'FontSize', fontSizeZoom)
set(ax2Zoom, 'TickLabelInterpreter', 'latex')
set(ax2Zoom, 'linewidth', lineWidthZoom)
set(ax2Zoom, 'XColor', 'k')
set(ax2Zoom, 'YColor', 'k')

fileName = strcat('ber_cp_cut_', num2str(cpLength), '_CP');
saveas(fig, [figuresFolder '/' fileName '.eps'], 'epsc')
end


function plot_ber_surf(rcWindowResults, rcMaskedWindowResults, ...
    optWindowResults, optMaskedWindowResults, typeOFDM)
% PLOT_BER_SURF     Plots the BER results as surfaces.
%

switch typeOFDM
    case {'wtx', 'wrx', 'CPwtx', 'CPwrx'}
        fileName = strcat('ber_surface_', typeOFDM);
        plot_double_surf(optWindowResults, rcWindowResults, fileName)
        fileNameMasked = strcat('ber_surface_masked_', typeOFDM);
        plot_double_surf(optMaskedWindowResults, rcMaskedWindowResults, ...
            fileNameMasked)
    case {'CPW', 'WOLA'}
        caseSet = {'1stepTx', '1stepRx', '2stepTx', '2stepRx', ...
            '3stepTx', '3stepRx'};
        for setIndex = 1:length(caseSet)
            fileName = strcat('ber_surface_', typeOFDM, '_', ...
                caseSet{setIndex});
            plot_double_surf(optWindowResults(:, :, setIndex), ...
                rcWindowResults, fileName)
            fileNameMasked = strcat('ber_surface_masked_', typeOFDM, ...
                '_', caseSet{setIndex});
            plot_double_surf(optMaskedWindowResults(:, :, setIndex), ...
                rcMaskedWindowResults, fileNameMasked)
        end
end
end


function plot_double_surf(matrixA , matrixB, fileName)
% PLOT_DOUBLE_SURF  Plots a double graph of surface for BER results
%

global figuresFolder settingsFileName

width = 22;
height = width*2/(1+sqrt(5));
fontSize = 14;
horizontalLeftDistance = 2.5;
verticalBottomDistance = 2;
plotWidth = width/2 - 2*horizontalLeftDistance;
plotHeight = height - 2*verticalBottomDistance;
settingsLoader = load(settingsFileName);
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrVector = settingsLoader.settingsData.generalSettings.snrValues;

fig = figure;
fig.Name = fileName;
fig.Color = 'w';
fig.Units = 'centimeters';
fig.Position = [10 10 width height];
ax1 = axes('units', 'centimeters', 'position', ...
    [horizontalLeftDistance+.5, verticalBottomDistance, plotWidth, ...
    plotHeight]);
ax2 = axes('units', 'centimeters', 'position', ...
    [width/2 + 1.5*horizontalLeftDistance, verticalBottomDistance, ...
    plotWidth, plotHeight]);
subplot(ax1)
surf(cpLengthVector, snrVector, 10*log10(matrixA.'), 'LineWidth', .1, ...
    'EdgeAlpha', .91)
set(ax1, 'Ydir', 'reverse')
set(ax1, 'Xdir', 'reverse')
set(ax1, 'FontSize', fontSize)
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'linewidth', 1)
ylim([-30 50])
xLabelAx1 = xlabel('CP Length, $\mu$', 'interpreter', 'latex', ...
    'FontSize', fontSize, 'Color', 'k');
xLabelAx1.Units = 'centimeters';
xLabelAx1.Position(1) = xLabelAx1.Position(1) - .6;
xLabelAx1.Position(2) = xLabelAx1.Position(2) - .1;
xLabelAx1.Rotation = 22.5;
yLabelAx1 = ylabel('SNR, dB', 'interpreter', 'latex', 'FontSize', ...
    fontSize, 'Color', 'k');
yLabelAx1.Units = 'centimeters';
yLabelAx1.Position(1) = yLabelAx1.Position(1) + .5;
yLabelAx1.Position(2) = yLabelAx1.Position(2) - .1;
yLabelAx1.Rotation = -37.5;
zlabel('BER, dB', 'interpreter', 'latex', 'FontSize', fontSize, ...
    'Color', 'k')
title('(a)', 'interpreter', 'latex', 'FontSize', fontSize, ...
    'Color', 'k')
subplot(ax2)
surf(cpLengthVector, snrVector, 10*log10(matrixB.'), 'LineWidth', .1, ...
    'EdgeAlpha', .91)
set(ax2, 'Ydir', 'reverse')
set(ax2, 'Xdir', 'reverse')
set(ax2, 'FontSize', fontSize)
set(ax2, 'TickLabelInterpreter', 'latex')
set(ax2, 'linewidth', 1)
ylim([-30 50])
title('(b)', 'interpreter', 'latex', 'FontSize', fontSize, 'Color', 'k')

saveas(fig, [figuresFolder '/' fileName '.eps'], 'epsc')
end


function [rcWindowResults, rcMaskedWindowResults, optWindowResults, ...
    optMaskedWindowResults] = read_by_type(typeOFDM)
% READ_BY_TYPE  Reads the resulting BER according to the selected system
% type.
%

global settingsFileName resultsFolder
settingsLoader = load(settingsFileName);
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrValues = settingsLoader.settingsData.generalSettings.snrValues;
resultsFiles = dir(fullfile(resultsFolder));
% Select files containing the same w-OFDM system type
switch typeOFDM
    case {'wtx', 'wrx', 'CPwtx', 'CPwrx'}
        rcWindowResults = zeros(length(cpLengthVector), length(snrValues));
        rcMaskedWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues));
        optWindowResults = zeros(length(cpLengthVector), length(snrValues));
        optMaskedWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues));
        for fileIndex = 1:length(resultsFiles)
            if resultsFiles(fileIndex).isdir
                continue
            end
            fileNameInfo = split(resultsFiles(fileIndex).name, '_');
            fileTypeOFDM = fileNameInfo{end-1};
            cpLength = str2double(fileNameInfo{end}(1:2));
            dataLoader = load([resultsFolder '/' ...
                resultsFiles(fileIndex).name]);
            if isequal(fileTypeOFDM, typeOFDM)
                if isequal(fileNameInfo{1}, 'rc')
                    rcWindowResults( ...
                        cpLengthVector==cpLength, :) = dataLoader.berRCSNR;
                elseif isequal(fileNameInfo{1}, 'optimized')
                    optWindowResults( ...
                        cpLengthVector==cpLength, :) = dataLoader.berSNR;
                else
                    if isequal(fileNameInfo{2}, 'rc')
                        rcMaskedWindowResults( ...
                            cpLengthVector==cpLength, : ...
                            ) = dataLoader.berMaskedRCSNR;
                    elseif isequal(fileNameInfo{2}, 'optimized')
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, : ...
                            ) = dataLoader.berMaskedSNR;
                    end
                end
            end
        end
    case {'WOLA', 'CPW'}
        rcWindowResults = zeros(length(cpLengthVector), length(snrValues));
        rcMaskedWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues));
        optWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues), 6);
        optMaskedWindowResults = zeros(length(cpLengthVector), ...
            length(snrValues), 6);
        for fileIndex = 1:length(resultsFiles)
            if resultsFiles(fileIndex).isdir
                continue
            end
            fileNameInfo = split(resultsFiles(fileIndex).name, '_');
            fileTypeOFDM = fileNameInfo{end-1};
            cpLength = str2double(fileNameInfo{end}(1:2));
            dataLoader = load([resultsFolder '/' ...
                resultsFiles(fileIndex).name]);
            if isequal(fileTypeOFDM, typeOFDM)
                if isequal(fileNameInfo{1}, 'rc')
                    rcWindowResults( ...
                        cpLengthVector==cpLength, :) = dataLoader.berRCSNR;
                elseif isequal(fileNameInfo{1}, 'optimized')
                    optWindowResults( ...
                        cpLengthVector==cpLength, :, 1 ...
                        ) = dataLoader.berSNRStep1A;
                    optWindowResults( ...
                        cpLengthVector==cpLength, :, 2 ...
                        ) = dataLoader.berSNRStep1B;
                    optWindowResults( ...
                        cpLengthVector==cpLength, :, 3 ...
                        ) = dataLoader.berSNRStep2A;
                    optWindowResults( ...
                        cpLengthVector==cpLength, :, 4 ...
                        ) = dataLoader.berSNRStep2B;
                    optWindowResults( ...
                        cpLengthVector==cpLength, :, 5 ...
                        ) = dataLoader.berSNRStep3A;
                    optWindowResults( ...
                        cpLengthVector==cpLength, :, 6 ...
                        ) = dataLoader.berSNRStep3B;
                else
                    if isequal(fileNameInfo{2}, 'rc')
                        rcMaskedWindowResults( ...
                            cpLengthVector==cpLength, : ...
                            ) = dataLoader.berMaskedRCSNR;
                    elseif isequal(fileNameInfo{2}, 'optimized')
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, :, 1 ...
                            ) = dataLoader.berMaskedSNRStep1A;
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, :, 2 ...
                            ) = dataLoader.berMaskedSNRStep1B;
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, :, 3 ...
                            ) = dataLoader.berMaskedSNRStep2A;
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, :, 4 ...
                            ) = dataLoader.berMaskedSNRStep2B;
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, :, 5 ...
                            ) = dataLoader.berMaskedSNRStep3A;
                        optMaskedWindowResults( ...
                            cpLengthVector==cpLength, :, 6 ...
                            ) = dataLoader.berMaskedSNRStep3B;
                    end
                end
            end
        end
end
end


function draw_box(topLeft, bottomRight, lineWidth)
% DRAW_BOX  Function to draw box representing zoom area.
%

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
