% main_interference_figures.m
%
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Oct 20, 2022
%


clc
clear
close all


fprintf('Starting main_interference_figures.m ... \n\n')

global figuresFolder settingsFileName resultsFolder listSystems

% Definitions
listSystems = {'wtx', 'wrx', 'CPwtx', 'CPwrx', 'WOLA', 'CPW'};
settingsFileName = 'settingsData.mat';
resultsFolder = 'interference_results';
figuresFolder = [resultsFolder '/figures'];
if ~isdir(resultsFolder)  %#ok
    message = strcat(resultsFolder, ' does not exist.');
    error(message)
end
if ~isdir(figuresFolder)  %#ok
    mkdir(figuresFolder)
end
resultFiles = dir(fullfile(resultsFolder));
for fileIndex = 1:length(resultFiles)
    if resultFiles(fileIndex).isdir
        continue
    end
    fprintf('Working on file %s.\n', resultFiles(fileIndex).name)
    resultInfo = split(resultFiles(fileIndex).name, '_');
    typeOFDM = resultInfo{end}(1:end-4);
    fprintf('Contains results for %s-OFDM.\n', typeOFDM)
    dataLoader = load([resultsFolder '/' resultFiles(fileIndex).name]);
    rcVector = dataLoader.rcWindowInterference;
    switch typeOFDM
        case {'wtx', 'wrx', 'CPwtx', 'CPwrx'} 
            optVector = dataLoader.optWindowInterference;
            plot_single_case(optVector, rcVector, strcat( ...
                'interference_', typeOFDM))
        case {'WOLA', 'CPW'}
            optVectorTxStep1 = dataLoader.optInterferenceCaseAStep1;
            plot_single_case(optVectorTxStep1, rcVector, strcat( ...
                'interferece_', typeOFDM, '_CaseAStep1'))
            optVectorTxStep2 = dataLoader.optInterferenceCaseAStep2;
            plot_single_case(optVectorTxStep2, rcVector, strcat( ...
                'interference_', typeOFDM, '_CaseAStep2'))
            optVectorTxStep3 = dataLoader.optInterferenceCaseAStep3;
            plot_single_case(optVectorTxStep3, rcVector, strcat( ...
                'interference_', typeOFDM, '_CaseAStep3'))
            optVectorRxStep1 = dataLoader.optInterferenceCaseBStep1;
            plot_single_case(optVectorRxStep1, rcVector, strcat( ...
                'interference_', typeOFDM, '_CaseBStep1'))
            optVectorRxStep2 = dataLoader.optInterferenceCaseBStep2;
            plot_single_case(optVectorRxStep2, rcVector, strcat( ...
                'interference_', typeOFDM, '_CaseBStep2'))
            optVectorRxStep3 = dataLoader.optInterferenceCaseBStep3;
            plot_single_case(optVectorRxStep3, rcVector, strcat( ...
                'interference_', typeOFDM, '_CaseBStep3'))
        otherwise
            fprintf('Skipping %s-OFDM', typeOFDM)
            continue
    end
end
% close all
plot_best_cases()


function plot_single_case(optVector, rcVector, fileName)
% Funtion to plot results for power interference.
%

global figuresFolder settingsFileName

width = 22;
height = width*2/(1+sqrt(5));
fontSize = 14;
lineWidth = 1;
markerSize = 15;
horizontalLeftDistance = 2.5;
verticalBottomDistance = 2;
plotWidth = width - 1.5*horizontalLeftDistance;
plotHeight = height - 1.5*verticalBottomDistance;
settingsLoader = load(settingsFileName);
cpVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;

fig = figure;
fig.Name = fileName;
fig.Color = 'w';
fig.Units = 'centimeters';
fig.Position = [10 10 width height];
ax = axes('units', 'centimeters', 'position', [horizontalLeftDistance, ...
    verticalBottomDistance, plotWidth, plotHeight]);
subplot(ax)
semilogy(cpVector, optVector, '*-', 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize), hold on
semilogy(cpVector, rcVector, 'd-', 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize), hold off, grid on
set(ax, 'FontSize', fontSize)
set(ax, 'TickLabelInterpreter', 'latex')
set(ax, 'linewidth', lineWidth)
set(ax, 'XColor', 'k')
set(ax, 'YColor', 'k')
xlabel('Cyclic Prefix Length, $\mu$', 'Interpreter', 'latex', ...
    'FontSize', fontSize)
ylabel('Interference Power', 'Interpreter', 'latex', 'FontSize', fontSize)
lgd = legend('Optimal', 'Raised Cosine');
lgd.Interpreter = 'latex';
lgd.FontSize = fontSize;


filePath = [figuresFolder '/' fileName];
saveas(fig, filePath, 'epsc')
end


function plot_best_cases()

global settingsFileName listSystems figuresFolder

width = 22;
height = width*2/(1+sqrt(5));
fontSize = 14;
lineWidth = 1.5;
markerSize = 15;
horizontalLeftDistance = 2.5;
verticalBottomDistance = 2;
plotWidth = width - 1.5*horizontalLeftDistance;
plotHeight = height - 1.5*verticalBottomDistance;
settingsLoader = load(settingsFileName);
CPVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
fileName = 'interf_results_best';

resultsConcat = get_best_results();

fig = figure;
fig.Name = fileName;
fig.Color = 'w';
fig.Units = 'centimeters';
fig.Position = [10 10 width height];
ax = axes('units', 'centimeters', 'position', [horizontalLeftDistance, ...
    verticalBottomDistance, plotWidth, plotHeight]);
subplot(ax)
indexedLineStyle = {'d-', 's-', 'o-', '*-', '+-', 'x-'};
for idx = 1:length(listSystems)
    semilogy(CPVector, resultsConcat(idx, :, 1), indexedLineStyle{idx}, ...
        'LineWidth', lineWidth, 'MarkerSize', markerSize)
    if idx == 1
        hold on
    end
end
ax.ColorOrderIndex = 1;
for idx = 1:length(listSystems)
    semilogy(CPVector, resultsConcat(idx, :, 2), ...
        strcat(indexedLineStyle{idx}, '-'), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize)
end
hold off, grid on
set(ax, 'FontSize', fontSize)
set(ax, 'TickLabelInterpreter', 'latex')
set(ax, 'linewidth', lineWidth)
set(ax, 'XColor', 'k')
set(ax, 'YColor', 'k')
xlabel('CP Length, $\mu$', 'interpreter', 'latex', 'FontSize', fontSize)
ylabel('Interference Power', 'interpreter', 'latex', 'FontSize', fontSize)
lgd = legend(listSystems);
lgd.Units = 'centimeters';
lgd.Interpreter = 'latex';
lgd.FontSize = fontSize;
lgd.Position(1) = horizontalLeftDistance + .25;
lgd.Position(2) = verticalBottomDistance + .25;
lgd.NumColumns = 1;
xlim([10 32])

fileName = 'best_performing_interf_power.eps';
filePath = [figuresFolder '/' fileName];
saveas(fig, filePath, 'epsc')
end


function [bestResults] = get_best_results()
% Function to get the best interference results.
%
% This function gets all results for wtx, wrx, CPwtx, and CPwrx, but
% selects the best results for CPW and WOLA systems. In this implementation
% we've observed that the best results are achieved at the third step of
% case B.
%
% Returns
% -------
%   bestResults: 3d array
%       3d array with best results for every system. We have the following
%       indexing: (typeOFDM, results, optimal or RC indexing)
%

global resultsFolder settingsFileName listSystems

fprintf('Extracting best results...\n')
settingsLoader = load(settingsFileName);
CPVectorLength = length(settingsLoader.settingsData.generalSettings.cyclicPrefix);
noSystems = length(listSystems);
bestResults = zeros(noSystems, CPVectorLength, 2);
listFromDirectory = dir(fullfile(resultsFolder));
for fileIndex = 1:length(listFromDirectory)
    if listFromDirectory(fileIndex).isdir
        continue
    end
    fileObject = listFromDirectory(fileIndex);
    resultInfo = split(fileObject.name, '_');
    typeOFDM = resultInfo{end}(1:end-4);
    dataLoader = load([resultsFolder '/' fileObject.name]);
    RCVector = dataLoader.rcWindowInterference;
    idx = find(strcmp(listSystems, typeOFDM));
    switch typeOFDM
        case {'wtx', 'wrx', 'CPwtx', 'CPwrx'} 
            optVector = dataLoader.optWindowInterference;
        case {'WOLA', 'CPW'}
            optVector = dataLoader.optInterferenceCaseBStep3;
        otherwise
            fprintf('Skipping %s-OFDM', typeOFDM)
            continue
    end
    bestResults(idx, :, 1) = optVector;
    bestResults(idx, :, 2) = RCVector;
end

end


% EoF

