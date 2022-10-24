% main_interference_figures.m
%
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Oct 20, 2022
%


clc
clear
close all


fprintf('Starting main_interference_figures.m ... \n\n')

global figuresFolder settingsFileName

% Definitions
settingsFileName = 'settingsData.mat';
resultsFolder = 'interference_results';
figuresFolder = [resultsFolder '/figures'];
if ~isfolder(resultsFolder)
    message = strcat(resultsFolder, ' does not exist.');
    error(message)
end
if ~isfolder(figuresFolder)
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


function plot_single_case(optVector, rcVector, fileName)
% Funtion to plot results for power interference.
%

global figuresFolder settingsFileName

width = 22;
height = width*2/(1+sqrt(5));
fontSize = 14;
horizontalLeftDistance = 2;
verticalBottomDistance = 3;
plotWidth = width - 2*horizontalLeftDistance;
plotHeight = height - 2*verticalBottomDistance;
settingsLoader = load(settingsFileName);
cpVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
fileInfo = split(fileName, '_');
typeOFDM = fileInfo{2};

fig = figure;
fig.Name = strcat(fileInfo{1}, ' ', typeOFDM);
fig.Color = 'w';
fig.Units = 'centimeters';
fig.Position = [10 10 width height];
ax1 = axes('units', 'centimeters', 'position', ...
    [horizontalLeftDistance, verticalBottomDistance, ...
    plotWidth, plotHeight]);

subplot(ax1)
semilogy(cpVector, optVector), hold on
semilogy(cpVector, rcVector), hold off, grid on

filePath = [figuresFolder '/' fileName];

end


% EoF

