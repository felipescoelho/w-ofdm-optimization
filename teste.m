clc
clear
close all

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
typeOFDMSet = {'wtx' 'wrx' 'CPwtx' 'CPwrx' 'WOLA' 'CPW'};
settingsLoader = load(settingsFileName);
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
snrValues = settingsLoader.settingsData.generalSettings.snrValues;


selectedResults = zeros(length(typeOFDMSet), length(snrValues), 2);
for idx = 1:length(resultFiles)
    if resultFiles(idx).isdir
        continue
    end
    fileNameInfo = split(resultFiles(idx).name, '_');
    typeOFDM = fileNameInfo{3};
    typeIndex = find(strcmp(typeOFDMSet, typeOFDM));
    switch typeOFDM
        case {'wtx', 'wrx', 'CPwrx', 'CPwtx'}
            if isequal(fileNameInfo{end}(1:2), '22')
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
            if isequal(fileNameInfo{end}(1:2), '22')
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
    disp(resultFiles(idx).name)
end