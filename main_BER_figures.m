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
resultsFolder = 'ber_results';
figuresFolder = [resultsFolder '/figures'];
figuresByCPFolder = [figuresFolder '/by_cyclic_prefix'];
figuresByTypeFolder = [filguresFolder '/by_type'];
if ~isfolder(figuresFolder)
    mkdir(figuresFolder)
end
if ~isfolder(figureByCPsFolder)
    mkdir(figuresByCPFolder)
end
if ~isfolder(figuresByTypeFolder)
    mkdir(figuresByTypeFolder)
end
resultFiles = dir(fullfile(resultsFolder));
typeOFDMSet = {'wtx' 'wrx' 'WOLA' 'CPW' 'CPwtx' 'CPwrx'};
for typeOFDMIndex = 1:length(typeOFDMSet)
    typeOFDM = typeOFDMSet{typeOFDMIndex};
    filesByType = file_by_type(resultFiles, typeOFDM);
    cpList = list_cp(filesByType);
end


function plot_ber_by_cp(cyclicPrefix, folderPath)
% Funtion to plot different systems for a single CP length.
%
% - Input
%   . cyclicPrefix: Lenght of the cyclic prefix
%   . folderPath: Where the figure will be saved

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


function filesByType = file_by_type(resultFiles, typeOFDM)
filesByType = {};
for fileIndex = 1:length(resultFiles)
    if resultFiles(fileIndex).isdir
        continue
    end
    fileNameInfo = split(resultFiles(fileIndex).name, '_');
    fileTypeOFDM = fileNameInfo{2};
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

