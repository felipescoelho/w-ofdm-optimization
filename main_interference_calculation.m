% main_interference_calculation.m
%   This script will perform the calculation of the ICI + ISI power for the
%   different systems and difference cyclic prefix lengths.
%
% Luiz Felipe Coelho - luizfelipe.coelho@smt.ufrj.br
% Out 17, 2022
%


clc
clear
close all


fprintf('Starting main_interference_calculation.m ... \n\n')

global folderName settingsFileName channelsFilePath
folderName = 'optimized_windows';
if ~isdir(folderName)
    error('Missing %s folder.', folderName)
end
settingsFileName = 'settingsData.mat';
% if ~isfile(settingsFileName)
%     error('Missing %s file.', settingsFileName)
% end
channelsFilePath = './channels/vehA200channel2.mat';
% if ~isfile(channelsFilePath)
%     error('Missing %s file.', channelsFilePath)
% end

folderToSave = 'interference_results';
if ~isdir(folderToSave)
    mkdir(folderToSave)
end
typeOFDMSet = {'wtx', 'wrx', 'WOLA', 'CPW', 'CPwrx', 'CPwtx'};
settingsLoader = load(settingsFileName);
numSubcar = settingsLoader.settingsData.generalSettings.numberSubcarriers;
cpLengthVector = settingsLoader.settingsData.generalSettings.cyclicPrefix;
alpha = .99;
for typeIndex = 1:length(typeOFDMSet)
    typeOFDM = typeOFDMSet{typeIndex};
    tailRx = settingsLoader.settingsData.(typeOFDM).tailRx;
    tailTx = settingsLoader.settingsData.(typeOFDM).tailTx;
    filesByType = file_by_type(typeOFDM);
    rcWindowInterference = zeros(length(cpLengthVector), 1);
    switch typeOFDM
        case {'wtx', 'CPwtx'}
            optWindowInterference = zeros(length(cpLengthVector), 1);
            for cpIndex = 1:length(cpLengthVector)
                cpLength = cpLengthVector(cpIndex);
                [csLength, ~, ~] = calculate_parameters(typeOFDM, ...
                    cpLength, tailTx, tailRx);
                fileName = select_by_cp(filesByType, cpLength);
		windowLoader = load([folderName '/' fileName]);
                windowTx = windowLoader.optimizedWindow;
                windowTxRC = transmitter_rc_window(numSubcar, cpLength, ...
                    csLength, tailTx);
                windowRx = receiver_rc_window(numSubcar, tailRx);
                optWindowInterference(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowTx, windowRx, alpha);
                rcWindowInterference(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowTxRC, windowRx, alpha);
            end
            fileToSave = strcat('interference_', typeOFDM);
            save([folderToSave '/' fileToSave], 'rcWindowInterference', ...
                'optWindowInterference')
        case {'wrx', 'CPwrx'}
            optWindowInterference = zeros(length(cpLengthVector), 1);
            for cpIndex = 1:length(cpLengthVector)
                cpLength = cpLengthVector(cpIndex);
                [csLength, ~, ~] = calculate_parameters(typeOFDM, ...
                    cpLength, tailTx, tailRx);
                fileName = select_by_cp(filesByType, cpLength);
                windowTx = transmitter_rc_window(numSubcar, cpLength, ...
                    csLength, tailTx);
	    	windowLoader = load([folderName '/' fileName]);
                windowRx = windowLoader.optimizedWindow;
                windowRxRC = receiver_rc_window(numSubcar, tailRx);
                optWindowInterference(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowTx, windowRx, alpha);
                rcWindowInterference(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowTx, windowRxRC, alpha);
            end
            fileToSave = strcat('interference_', typeOFDM);
            save([folderToSave '/' fileToSave], 'rcWindowInterference', ...
                'optWindowInterference')
        case {'WOLA', 'CPW'}
            optInterferenceCaseAStep1 = zeros(length(cpLengthVector), 1);
            optInterferenceCaseAStep2 = zeros(length(cpLengthVector), 1);
            optInterferenceCaseAStep3 = zeros(length(cpLengthVector), 1);
            optInterferenceCaseBStep1 = zeros(length(cpLengthVector), 1);
            optInterferenceCaseBStep2 = zeros(length(cpLengthVector), 1);
            optInterferenceCaseBStep3 = zeros(length(cpLengthVector), 1);
            for cpIndex = 1:length(cpLengthVector)
                cpLength = cpLengthVector(cpIndex);
                [csLength, ~, ~] = calculate_parameters(typeOFDM, ...
                    cpLength, tailTx, tailRx);
                fileName = select_by_cp(filesByType, cpLength);
                windowTxRC = transmitter_rc_window(numSubcar, cpLength, ...
                    csLength, tailTx);
                windowRxRC = receiver_rc_window(numSubcar, tailRx);
                windowLoader = load([folderName '/' fileName]);
                windowCaseAStep1 = windowLoader.optimizedWindowCaseAStep1;
                windowCaseAStep2 = windowLoader.optimizedWindowCaseAStep2;
                windowCaseAStep3 = windowLoader.optimizedWindowCaseAStep3;
                windowCaseBStep1 = windowLoader.optimizedWindowCaseBStep1;
                windowCaseBStep2 = windowLoader.optimizedWindowCaseBStep2;
                windowCaseBStep3 = windowLoader.optimizedWindowCaseBStep3;
                optInterferenceCaseAStep1(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowCaseAStep1, windowRxRC, alpha);
                optInterferenceCaseAStep2(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowCaseAStep1, windowCaseAStep2, alpha);
                optInterferenceCaseAStep3(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowCaseAStep3, windowCaseAStep2, alpha);
                optInterferenceCaseBStep1(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowTxRC, windowCaseBStep1, alpha);
                optInterferenceCaseBStep2(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowCaseBStep2, windowCaseBStep1, alpha);
                optInterferenceCaseBStep3(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowCaseBStep2, windowCaseBStep3, alpha);
                rcWindowInterference(cpIndex) = ...
                    calculate_interference(cpLength, typeOFDM, ...
                    windowTxRC, windowRxRC, alpha);
            end
            fileToSave = strcat('interference_', typeOFDM);
            save([folderToSave '/' fileToSave], 'rcWindowInterference', ...
                'optInterferenceCaseAStep1', ...
                'optInterferenceCaseAStep2', ...
                'optInterferenceCaseAStep3', ...
                'optInterferenceCaseBStep1', ...
                'optInterferenceCaseBStep2', ...
                'optInterferenceCaseBStep3')
    end
end


function filesByType = file_by_type(typeOFDM)

global folderName

windowFiles = dir(fullfile(folderName));
filesByType = {};
for fileIndex = 1:length(windowFiles)
    if windowFiles(fileIndex).isdir
        continue
    end
    fileNameInfo = split(windowFiles(fileIndex).name, '_');
    fileTypeOFDM = fileNameInfo{3};
    if isequal(fileTypeOFDM, typeOFDM)
        filesByType(end+1) = {windowFiles(fileIndex).name};  %#ok
    end
end
end


function fileName = select_by_cp(filesByType, cpLength)

for fileIndex = 1:length(filesByType)
    fileInfo = split(filesByType{fileIndex}, '_');
    cpFile = str2double(fileInfo{end}(1:end-6));
    if cpFile == cpLength
        fileName = filesByType{fileIndex};
    end
end
end


function interferencePower = calculate_interference(cpLength, typeOFDM, ...
    windowTx, windowRx, alpha)
% Funtion to calculate interference for a given system
%
% - Input
%   . cpLength: Number of samples in the cyclic prefix
%   . typeOFDM: Type of w-OFDM system

global settingsFileName channelsFilePath

settingsLoader = load(settingsFileName);
channelLoader = load(channelsFilePath);
% System information
numSubcar = settingsLoader.settingsData.generalSettings.numberSubcarriers;
tailTx = settingsLoader.settingsData.(typeOFDM).tailTx;
tailRx = settingsLoader.settingsData.(typeOFDM).tailRx;
[csLength, prefixRemovalLength, circularShiftLength] = calculate_parameters( ...
    typeOFDM, cpLength, tailTx, tailRx);
% Channel inforamtion
meanChannelImpulseResponse = mean(channelLoader.vehA200channel2, 1);
interferenceArray = array_ici_isi(meanChannelImpulseResponse, ...
    numSubcar, cpLength, csLength, tailTx, tailRx, prefixRemovalLength);
intercarrierInterference = interferenceArray(:, :, 1);
intersymbolInterference = sum(interferenceArray(:, :, 2:end), 3);
% Generating system matrices for calculations
transformMatrix = transform_matrix(numSubcar);
invertTransformMatrix = conj(transformMatrix)/numSubcar;
circularShiftMatrix = circular_shift_matrix(numSubcar, circularShiftLength);
overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx);
removeRedundancyMatrix = remove_redundancy_matrix(numSubcar, tailRx, ...
    prefixRemovalLength);
addRedundancyMatrix = add_redundancy_matrix(numSubcar, cpLength, csLength);
matrixA0 = transformMatrix*circularShiftMatrix*overlapAddMatrix ...
    * windowRx*removeRedundancyMatrix*intercarrierInterference ...
    * windowTx*addRedundancyMatrix*invertTransformMatrix;
matrixAm = transformMatrix*circularShiftMatrix*overlapAddMatrix ...
    * windowRx*removeRedundancyMatrix*intersymbolInterference ...
    * windowTx*addRedundancyMatrix*invertTransformMatrix;
matrixA0ICI1 = matrixA0 - diag(diag(matrixA0));
matrixAmISI = diag(diag(matrixAm));
matrixAmICI2 = matrixAm - matrixAmISI;
PICI1 = alpha*diag(matrixA0ICI1*matrixA0ICI1');
PISI = (1-alpha)*diag(matrixAmISI*matrixAmISI');
PICI2 = (1-alpha)*diag(matrixAmICI2*matrixAmICI2');

interferencePower = sum(PICI1 + PISI + PICI2);
end


function overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx)
% Function to generate matrix that operates the overlap and add in the
% samples.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tails for the receiver
%   window.

identityHalfTail = eye(tailRx/2);
identitySubcarriers = eye(numSubcar-tailRx);
zerosHalfTail = zeros(tailRx/2);
zerosHalfTailSubcarriers = zeros(tailRx/2, numSubcar-tailRx);
zerosTailSubcarriers = zeros(numSubcar-tailRx, tailRx);
overlapAddMatrix = [zerosHalfTail identityHalfTail ...
    zerosHalfTailSubcarriers zerosHalfTail identityHalfTail; ...
    zerosTailSubcarriers identitySubcarriers zerosTailSubcarriers; ...
    identityHalfTail zerosHalfTail zerosHalfTailSubcarriers ...
    identityHalfTail zerosHalfTail];
end


function circularShiftMatrix = circular_shift_matrix(numSubcar, ...
    circularShiftLength)
% Function to generate matrix that operates the circular shift.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cicularShift: Samples in the circular shift.
%
% - Output:
%   . circularShiftMatrix: A matrix capable of performing the circular
%   shift operation.

identityCircularShift = eye(circularShiftLength);
identitySubcarriersMinusCircularShift = eye(numSubcar-circularShiftLength);
zerosSubcarriersMinusCircularShift = zeros(circularShiftLength, ...
    numSubcar-circularShiftLength);
circularShiftMatrix = [zerosSubcarriersMinusCircularShift.' ...
    identitySubcarriersMinusCircularShift; identityCircularShift ...
    zerosSubcarriersMinusCircularShift];
end


function removeRedundancyMatrix = remove_redundancy_matrix(numSubcar, ...
    tailRx, prefixRemovalLength)
% Function to generate matrix to remove redundancy.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail at the receiver
%   window.
%   . prefixRemoval: Number of samples to be removed at the reception.

removeRedundancyMatrix = [zeros(numSubcar+tailRx, prefixRemovalLength) ...
    eye(numSubcar+tailRx)];
end


function addRedundancyMatrix = add_redundancy_matrix(numSubcar, ...
    cpLength, csLength)
% Function to generate the matrix that adds redundancy.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cpLength: Number of samples in cyclic prefix.
%   . csLength: Number of samples in cyclic suffix.
%
% - Output:
%   . redundancyMatrix: Matrix that adds the cyclic redundancy to the
%   system.

identityCP = eye(cpLength);
zerosCPVsNumSubcarMinCP = zeros(cpLength, (numSubcar-cpLength));
identityCS = eye(csLength);
zerosCSVsNumSubcarMinCS = zeros(csLength, (numSubcar-csLength));
identityNumSubcar = eye(numSubcar);
addRedundancyMatrix = [zerosCPVsNumSubcarMinCP identityCP; ...
    identityNumSubcar; ...
    identityCS zerosCSVsNumSubcarMinCS];
end


function transmitterWindow = transmitter_rc_window(numSubcar, cpLength, ...
    csLength, tailTx)
% Function to generate the Raised Cosine transmitter window.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cyclicPrefix: Number of samples in cyclic prefix.
%   . cyclicSuffix: Number of samples in cyclic suffix.
%   . tailTx: Number of samples in rise and fall tail for the transmitter
%   window.
%
% - Output:
%   . transmitterWindow: Matrix with diagonal equivalent to the transmitted
%   window.

raisedCosineAxis = -(tailTx+1)/2 + 1:(tailTx+1)/2;
raisedCosine = sin(pi/2*(0.5+raisedCosineAxis(1:end-1)/tailTx)).^2;
onesTransmitted = ones(1, numSubcar+cpLength+csLength ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
transmitterWindow = diag(windowVector);
end


function receiverWindow = receiver_rc_window(numSubcar, tailRx)
% Function to generate the Raised Cosine receiver window.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail for the receiver
%   window.

raisedCosineAxis = -(tailRx+1)/2+1:(tailRx+1)/2 ;
raisedCosine = sin(pi/2*(0.5+raisedCosineAxis(1:end-1)/tailRx)).^2;
onesReceived = ones(1, numSubcar-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
receiverWindow = diag(windowVector);
end


function interferenceArray = array_ici_isi(channelVector, numSubcar, ...
    cpLength, csLength, tailTx, tailRx, prefixRemovalLength)
% Function to calculate the ICI and ISI array.
%
% - Input:
%   . channelVector: Channel's impulse response.
%   . numberSubcarrier: Number of subcarriers in the system.
%   . cyclicPrefix: Number of samples in cyclic prefix.
%   . cyclicSuffix: Number of samples in cyclic suffix.
%   . tailTx: Number of samples in rise and fall tail at transmitter.
%   . tailRx: Number of samples in rise and fall tail at receiver.
%   . prefixRemoval: Number of samples (from cyclic prefix) to be removed
%   at the receiver.
%
% - Output:
%   . interferenceArray : An array with elements corresponding to the
%   channel effects in the OFDM block, where the  

chanOrder = length(channelVector) - 1;
numRx = numSubcar+tailRx+prefixRemovalLength;
numTx = numSubcar+cpLength+csLength;
numAffected = ceil((chanOrder+tailTx)/numRx);
channelNoise = numTx-tailTx;
interferenceArray = zeros(numRx, numTx, numAffected+1);
for symbAffected = 0:numAffected
    for rxSample = 0:numRx-1
        for txSample = 0:numTx-1
            indexer = symbAffected*channelNoise + rxSample - txSample;
            if (0 <= indexer) && (indexer <= chanOrder)
                interferenceArray(rxSample+1, txSample+1, symbAffected+1) ...
                    = channelVector(indexer+1);
            end
        end
    end
end
end


function transformMatrix = transform_matrix(numSubcar)

transformMatrix = zeros(numSubcar);
for row = 1:numSubcar
    for col = 1:numSubcar
        transformMatrix(row, col) = exp(-1j*(row-1)*2*pi*(col-1)/numSubcar);
    end
end
end


function [csLength, prefixRemoval, circularShift] = ...
    calculate_parameters(typeOFDM, cpLength, tailTx, tailRx)
% Function to calculate parameters for each w-OFDM system
%
% - Input:
%   . typeOFDM: w-OFDM system
%   . cyclicPrefix: Number of samples in the cyclic prefix
%   . tailTx: Number of samples in rise and fall tail for Tx
%   . tailRX: Number of samples in rise and fall tail for Rx

switch typeOFDM
    case 'wtx'
        csLength = tailTx;
        prefixRemoval = cpLength;
        circularShift = 0;
    case 'wrx'
        csLength = tailRx/2;
        prefixRemoval = cpLength - tailRx/2;
        circularShift = 0;
    case 'WOLA'
        csLength = tailTx;
        prefixRemoval = cpLength - tailRx;
        circularShift = tailRx/2;
    case 'CPW'
        csLength = tailTx + tailRx/2;
        prefixRemoval = cpLength - tailRx/2;
        circularShift = 0;
    case 'CPwtx'
        csLength = 0;
        prefixRemoval = cpLength - tailTx;
        circularShift = tailTx;
    case 'CPwrx'
        csLength = 0;
        prefixRemoval = cpLength - tailRx;
        circularShift = tailRx/2;
end
end


% EoF

