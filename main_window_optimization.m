% main_window_optimization.m
%   This script will perform the optimization of the windowed-OFDM for all
%   w-OFDM systems wtx, wrx, WOLA, CPW, CPwtx, and CPwrx.
%
% Author: Luiz Felipe Coelho - luizfelipe.coelho@smt.ufrj.br
% Set 30, 2022
%


clc
clear
close all


fprintf('Starting main_window_optimization.m ... \n\n')
fprintf('This script runs the optimization process for the windows\n')
fprintf('in the different w-OFDM systems.\n\n')


global folderName


folderName = 'optimized_windows';
if ~isdir(folderName)  %#ok
    mkdir(folderName)
end
settingsFileName = 'settingsData.mat';
% if isfile(settingsFileName)
%     delete(settingsFileName)
% end
fprintf('Generating settingsData structure array.\n')
settingsData = struct();
settingsData.generalSettings = struct('numberSubcarriers', 256, ...
    'bitsPerSubcarrier', 4, 'cyclicPrefix', 10:2:32, 'symbolsPerTx', 16, ...
    'ensemble', 1, 'snrValues', linspace(-20, 50, 30));
settingsData.wtx = struct('tailTx', 8, 'tailRx', 0);
settingsData.wrx = struct('tailTx', 0, 'tailRx', 10);
settingsData.WOLA = struct('tailTx', 8, 'tailRx', 10);
settingsData.CPW = struct('tailTx', 8, 'tailRx', 10);
settingsData.CPwtx = struct('tailTx', 8, 'tailRx', 0);
settingsData.CPwrx = struct('tailTx', 0, 'tailRx', 10);
fprintf('Saving settingsData structure array in %s\n', settingsFileName)
save(settingsFileName, 'settingsData')
fprintf('Successfully saved settingsData.\n')
fprintf('------------------------------------\n')
fprintf('Starting run ...\n\n')
fields = fieldnames(settingsData);
clear settingsData


alpha = .99;
for fieldIndex = 1:length(fields)
    fieldName = fields{fieldIndex};
    dataLoader = load(settingsFileName);
    numSubcar = dataLoader.settingsData.generalSettings.numberSubcarriers;
    cpLengthVect = dataLoader.settingsData.generalSettings.cyclicPrefix;
    fprintf('Working on %s-OFDM ...\n', fieldName)
    switch fieldName
        case 'wtx'
            tailTx = dataLoader.settingsData.(fieldName).tailTx;
            tailRx = dataLoader.settingsData.(fieldName).tailRx;
            for cpIndex = 1:length(cpLengthVect)
                cpLength = cpLengthVect(cpIndex);
                csLength = tailTx;
                prefixRemovalLength = cpLength;
                circularShiftLength = 0;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', numSubcar)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    numSubcar, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength, alpha)
            end
        case 'wrx'
            tailTx = dataLoader.settingsData.(fieldName).tailTx;
            tailRx = dataLoader.settingsData.(fieldName).tailRx;
            for cpIndex = 1:length(cpLengthVect)
                cpLength = cpLengthVect(cpIndex);
                csLength = tailRx/2;
                prefixRemovalLength = cpLength - tailRx/2;
                circularShiftLength = 0;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', numSubcar)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    numSubcar, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength, alpha)
            end
        case 'WOLA'
            tailTx = dataLoader.settingsData.(fieldName).tailTx;
            tailRx = dataLoader.settingsData.(fieldName).tailRx;
            for cpIndex = 1:length(cpLengthVect)
                cpLength = cpLengthVect(cpIndex);
                csLength = tailTx;
                prefixRemovalLength = cpLength - tailRx;
                circularShiftLength = tailRx/2;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', numSubcar)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    numSubcar, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength, alpha)
            end
        case 'CPW'
            tailTx = dataLoader.settingsData.(fieldName).tailTx;
            tailRx = dataLoader.settingsData.(fieldName).tailRx;
            for cpIndex = 1:length(cpLengthVect)
                cpLength = cpLengthVect(cpIndex);
                csLength = tailTx + tailRx/2;
                prefixRemovalLength = cpLength - tailRx/2;
                circularShiftLength = 0;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', numSubcar)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    numSubcar, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength, alpha)
            end
        case 'CPwtx'
            tailTx = dataLoader.settingsData.(fieldName).tailTx;
            tailRx = dataLoader.settingsData.(fieldName).tailRx;
            for cpIndex = 1:length(cpLengthVect)
                cpLength = cpLengthVect(cpIndex);
                csLength = 0;
                prefixRemovalLength = cpLength - tailTx;
                circularShiftLength = tailTx;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', numSubcar)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    numSubcar, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength, alpha)
            end
        case 'CPwrx'
            tailTx = dataLoader.settingsData.(fieldName).tailTx;
            tailRx = dataLoader.settingsData.(fieldName).tailRx;
            for cpIndex = 1:length(cpLengthVect)
                cpLength = cpLengthVect(cpIndex);
                csLength = 0;
                prefixRemovalLength = cpLength - tailRx;
                circularShiftLength = tailRx/2;
                fprintf('--------------------------------------------\n')
                fprintf('Settings:\n')
                fprintf('N = %u (Subcarriers)\n', numSubcar)
                fprintf('mu = %u (CP length)\n', cpLength)
                fprintf('rho = %u (CS length)\n', csLength)
                fprintf('beta = %u (Tx tail)\n', tailTx)
                fprintf('delta = %u (Rx tail)\n', tailRx)
                fprintf('gamma = %u (Samples removed Rx)\n', ...
                    prefixRemovalLength)
                fprintf('kappa = %u (Circ shift)\n', circularShiftLength)
                optimize_window(fieldName, cpLength, csLength, ...
                    numSubcar, tailTx, tailRx, prefixRemovalLength, ...
                    circularShiftLength, alpha)
            end
        otherwise
            continue
    end
end


function optimize_window(typeOFDM, cpLength, csLength, numSubcar, ...
    tailTx, tailRx, prefixRemovalLength, circularShiftLength, alpha)
% Function to generate the optimized windows.
%
% - Input:
%   . typeOFDM: Type of OFDM system
%   . cpLength: Number of samples in cyclic prefix (mu)
%   . csLength: Number of samples in cyclic suffix (rho)
%   . numSubcar: Number of subcarriers in system (N)
%   . tailTx: Number of samples in rise and fall tail at transmitter (beta)
%   . tailRx: Number of samples in rise and fall at receiver (delta)
%   . prefixRemovalLength: Number of samples (from cyclic prefix) to be
%   removed at the receiver. (gamma)
%   . circularShiftLength: Number of samples in circular shift, for WOLA,
%   CPwtx, and CPwrx. (kappa)

global folderName

folderPath = [cd '/' folderName '/'];
load('./channels/vehA200channel2.mat', 'vehA200channel2')
averageChannel = mean(vehA200channel2, 1);
interferenceArray = array_ici_isi(averageChannel, numSubcar, ...
    cpLength, csLength, tailTx, tailRx, prefixRemovalLength);
intercarrierInterference = interferenceArray(:, :, 1);
intersymbolInterference = sum(interferenceArray(:, :, 2:end), 3);
windowTxRC = transmitter_rc_window(numSubcar, cpLength, csLength, tailTx);
windowRxRC = receiver_rc_window(numSubcar, tailRx);
options = optimoptions(@fmincon);
options.Algorithm = 'sqp';
options.Display = 'iter-detailed';
options.MaxIterations = 1000;
options.MaxFunctionEvaluations = 6000;
options.UseParallel = true;
switch typeOFDM
    case {'wtx', 'CPwtx'}
        objectiveFunction = objective_function_tx(windowRxRC, numSubcar, ...
            tailRx, prefixRemovalLength, circularShiftLength, cpLength, ...
            csLength, intercarrierInterference, intersymbolInterference, ...
            alpha);
        initialValue = diag(windowTxRC);
        lowerBounds = [zeros(tailTx, 1); ...
            ones(numSubcar+cpLength+csLength-2*tailTx, 1); ...
            zeros(tailTx, 1)];
        upperBounds = [ones(tailTx, 1); ...
            ones(numSubcar+cpLength+csLength-2*tailTx, 1); ...
            ones(tailTx, 1)];
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                [windowVector, ~, exitFlag] = fmincon(objectiveFunction, initialValue, ...
                    [], [], [], [], lowerBounds, upperBounds, [], options);
                initialValue = windowVector;
            end
        end
        optimizedWindow = diag(windowVector);
        fprintf('Finished optimization process\n')
        % Save results
        fprintf('Saving results...\n')
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cpLength), 'CP');
        save([folderPath fileName], 'optimizedWindow')
        fprintf('Optimized windows were successfully saved into file.\n')
        fprintf('\n')
    case {'wrx', 'CPwrx'}
        objectiveFunction = objective_function_rx(windowTxRC, ...
            numSubcar, cpLength, csLength, circularShiftLength, tailRx, ...
            prefixRemovalLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = diag(windowRxRC);
        lowerBounds = [zeros(tailRx, 1); ...
            ones(numSubcar-tailRx, 1); zeros(tailRx, 1)];
        upperBounds = [ones(tailRx, 1); ...
            ones(numSubcar-tailRx, 1); ones(tailRx, 1)];
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                [windowVector, ~, exitFlag] = fmincon(objectiveFunction, ...
                    initialValue, [], [], [], [], lowerBounds, ...
                    upperBounds, [], options);
                initialValue = windowVector;
            end
        end
        optimizedWindow = diag(windowVector);
        fprintf('Finished optimization process\n')
        % Save results
        fprintf('Saving results...\n')
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cpLength), 'CP');
        save([folderPath fileName], 'optimizedWindow')
        fprintf('Optimized windows were successfully saved into file.\n')
        fprintf('\n')
    case {'WOLA', 'CPW'}
        lowerBoundsTx = [zeros(tailTx, 1); ...
            ones(numSubcar+cpLength+csLength-2*tailTx, 1); ...
            zeros(tailTx, 1)];
        upperBoundsTx = [ones(tailTx, 1); ...
            ones(numSubcar+cpLength+csLength-2*tailTx, 1); ...
            ones(tailTx, 1)];
        lowerBoundsRx = [zeros(tailRx, 1); ...
            ones(numSubcar-tailRx, 1); zeros(tailRx, 1)];
        upperBoundsRx = [ones(tailRx, 1); ...
            ones(numSubcar-tailRx, 1); ones(tailRx, 1)];
        
        objectiveFunctionCaseAStep1 = objective_function_tx(windowRxRC, ...
            numSubcar, tailRx, prefixRemovalLength, circularShiftLength, ...
            cpLength, csLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = diag(windowTxRC);
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                windowVectorAStep1 = fmincon(objectiveFunctionCaseAStep1, ...
                    initialValue, [], [], [], [], lowerBoundsTx, ...
                    upperBoundsTx, [], options);
                initialValue = windowVectorAStep1;
            end
        end
        optimizedWindowCaseAStep1 = diag(windowVectorAStep1);
        
        objectiveFunctionCaseAStep2 = objective_function_rx(windowTxRC, ...
            numSubcar, cpLength, csLength, circularShiftLength, tailRx, ...
            prefixRemovalLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = diag(windowRxRC);
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                windowVectorAStep2 = fmincon(objectiveFunctionCaseAStep2, ...
                    initialValue, [], [], [], [], lowerBoundsRx, ...
                    upperBoundsRx, [], options);
                initialValue = windowVectorAStep2;
            end
        end
        optimizedWindowCaseAStep2 = diag(windowVectorAStep2);
        
        objectiveFunctionCaseAStep3 = objective_function_tx(windowRxRC, ...
            numSubcar, tailRx, prefixRemovalLength, circularShiftLength, ...
            cpLength, csLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = windowVectorAStep1;
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                windowVectorAStep3 = fmincon(objectiveFunctionCaseAStep3, ...
                    initialValue, [], [], [], [], lowerBoundsTx, ...
                    upperBoundsTx, [], options);
                initalValue = windowVectorAStep3;
            end
        end
        optimizedWindowCaseAStep3 = diag(windowVectorAStep3);
        
        objectiveFunctionCaseBStep1 = objective_function_rx(windowTxRC, ...
            numSubcar, cpLength, csLength, circularShiftLength, tailRx, ...
            prefixRemovalLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = diag(windowRxRC);
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                windowVectorBStep1 = fmincon(objectiveFunctionCaseBStep1, ...
                    initialValue, [], [], [], [], lowerBoundsRx, ...
                    upperBoundsRx, [], options);
                initialValue = windowVectorBStep1;
            end
        end
        optimizedWindowCaseBStep1 = diag(windowVectorBStep1);
        
        objectiveFunctionCaseBStep2 = objective_function_tx(windowRxRC, ...
            numSubcar, tailRx, prefixRemovalLength, circularShiftLength, ...
            cpLength, csLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = diag(windowTxRC);
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                windowVectorBStep2 = fmincon(objectiveFunctionCaseBStep2, ...
                    initialValue, [], [], [], [], lowerBoundsTx, ...
                    upperBoundsTx, [], options);
                initialValue = windowVectorBStep2;
            end
        end
        optimizedWindowCaseBStep2 = diag(windowVectorBStep2);
        
        objectiveFunctionCaseBStep3 = objective_function_rx(windowTxRC, ...
            numSubcar, cpLength, csLength, circularShiftLength, tailRx, ...
            prefixRemovalLength, intercarrierInterference, ...
            intersymbolInterference, alpha);
        initialValue = windowVectorBStep1;
        tryCount = 0;
        exitFlag = 0;
        while tryCount < 5
            tryCount = tryCount + 1;
            if exitFlag == 0
                windowVectorBStep3 = fmincon(objectiveFunctionCaseBStep3, ...
                    initialValue, [], [], [], [], lowerBoundsRx, ...
                    upperBoundsRx, [], options);
                initialValue = windowVectorBStep3;
            end
        end
        optimizedWindowCaseBStep3 = diag(windowVectorBStep3);
        
        % Save results
        fprintf('Saving results...\n')
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cpLength), 'CP');
        save([folderPath fileName], 'optimizedWindowCaseAStep1', ...
            'optimizedWindowCaseAStep2', 'optimizedWindowCaseAStep3', ...
            'optimizedWindowCaseBStep1','optimizedWindowCaseBStep2', ...
            'optimizedWindowCaseBStep3')
        fprintf('Optimized windows were successfully saved into file.')
        fprintf('\n')
end

end


function objectiveFunction = objective_function_tx(windowRx, ...
    numSubcar, tailRx, prefixRemovalLength, circularShiftLength, ...
    cpLength, csLength, intercarrierInterference, ...
    intersymbolInterference, alpha)
% Function to help writing the optimization functions for Tx
%
% - Input:
%   . windowRx: Window used at the receiver
%   . transformMatrix: DFT matrix
%   . circularShiftMatrix: Matrix to perform the circular shift
%   . overlapAdd: Matrix to perform the overlap and add operation
%   . redundancyRemove: Matrix to remove redundancy at reception
%   . redundancyMatrix: Matrix to add redundancy at transmission
%   . invertTransform: IDFT matrix
%   . intercarrierInterference: Matrix with intercarrier interference for
%   the channel
%   . intersymbolInterference: Matrix with intersymbol interference for the
%   channel.
%
% - Output:
%   . optimizationFunction: Function to perform the optimization process
%   using MATLAB's fmincon function.


receiverMatrix = rx_wofdm_matrix(windowRx, numSubcar, tailRx, ...
    prefixRemovalLength, circularShiftLength);
addRedundancyMatrix = add_redundancy_matrix(numSubcar, cpLength, csLength);
invTransformMatrix = transform_matrix(numSubcar)'/numSubcar;

objectiveFunction = @(x) alpha*norm(receiverMatrix ...
    * intercarrierInterference*(diag(x)*addRedundancyMatrix ...
    * invTransformMatrix) - diag(diag(receiverMatrix ...
    * intercarrierInterference*(diag(x)*addRedundancyMatrix ...
    * invTransformMatrix))), 'fro') + (1-alpha) ...
    * norm(receiverMatrix*intersymbolInterference ...
    * (diag(x)*addRedundancyMatrix*invTransformMatrix), 'fro');
end


function objectiveFunction = objective_function_rx(windowTx, numSubcar, ...
    cpLength, csLength, circularShiftLength, tailRx, ...
    prefixRemovalLength, intercarrierInterference, ...
    intersymbolInterference, alpha)
% Funtion to help writing tje optimization functions for Rx
%
% - Input:
%   . windowTx: Window used at the receiver
%   . numSubcar: Number of subcarriers in the system
%   . cpLength: Length of the cyclic prefix
%   . csLength: Length of the cyclic suffix
%   . circularShiftLength: Length of the circular shift
%   . tailRx: Length of the rise and fall tail in reception
%   . prefixRemovalLength: Number of samples removed in reception
%   . intersymbolInterference: Matrix with intersymbol interference for the
%   channel.
%
% - Output:
%   . optimizationFunction: Function to perform the optimization process
%   using MATLAB's fmincon function.


transmitterMatrix = tx_wofdm_matrix(windowTx, numSubcar, cpLength, csLength);
transformMatrix = transform_matrix(numSubcar);
circularShiftMatrix = circular_shift_matrix(numSubcar, circularShiftLength);
overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx);
removeRedundancyMatrix = remove_redundancy_matrix(numSubcar, tailRx, ...
    prefixRemovalLength);

objectiveFunction = @(x) alpha*norm((transformMatrix*circularShiftMatrix ...
    * overlapAddMatrix*diag(x)*removeRedundancyMatrix) ...
    * intercarrierInterference*transmitterMatrix - diag(diag(( ...
    transformMatrix*circularShiftMatrix*overlapAddMatrix*diag(x) ...
    * removeRedundancyMatrix)*intercarrierInterference ...
    * transmitterMatrix)), 'fro') + (1-alpha)*norm((transformMatrix ...
    * circularShiftMatrix*overlapAddMatrix*diag(x) ...
    * removeRedundancyMatrix)*intersymbolInterference ...
    * transmitterMatrix, 'fro');
end


function interfArray = array_ici_isi(channelVector, numSubcar, ...
    cpLength, csLength, tailTx, tailRx, prefixRemovalLength)
% Calculates the array corresponding to the channel's effects on the 
% system.
%
% Parameters
% ----------
% channelVector : Vector (row or column)
%   Channel's impulse response.
% numberSubcarrier : double
%   Number of subcarriers in the system.
% cpLength : double 
%   Number of samples in cyclic prefix.
% csLength : double
%   Number of samples in cyclic suffix.
% tailTx : double
%   Number of samples in rise and fall tail at transmitter.
% tailRx : double
%   Number of samples in rise and fall tail at receiver.
% prefixRemovalLength : double
%   Number of samples (from cyclic prefix) to be removed at the receiver.

chanOrder = length(channelVector) - 1;
numRx = numSubcar+tailRx+prefixRemovalLength;
numTx = numSubcar+cpLength+csLength;
numAffected = ceil((chanOrder+tailTx)/numRx);
chanNoise = numTx-tailTx;
interfArray = zeros(numRx, numTx, numAffected+1);
for symbAffected = 0:numAffected
    for rxSample = 0:numRx-1
        for txSample = 0:numTx-1
            indexer = symbAffected*chanNoise + rxSample - txSample;
            if (0 <= indexer) && (indexer <= chanOrder)
                interfArray(rxSample+1, txSample+1, symbAffected+1) ...
                    = channelVector(indexer+1);
            end
        end
    end
end
end


function transmissionMatrix = tx_wofdm_matrix(windowTx, numSubcar, ...
    cpLength, csLength)
% Function to generate a matrix to operate the transmission process.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cpLength: Number of samples in cyclic prefix.
%   . csLength: Number of samples in cyclic suffix.
%   . tailTx : Number of samples in transmitter window tails.
%
% - Output:
%   . transmitterMatrix: Matrix that operates the transmission process.

transformMatrix = transform_matrix(numSubcar);
invTransformMatrix = transformMatrix'/numSubcar;
addRedundancyMatrix = add_redundancy_matrix(numSubcar, cpLength, csLength);
transmissionMatrix = windowTx*addRedundancyMatrix*invTransformMatrix;
end


function receptionMatrix = rx_wofdm_matrix(windowRx, numSubcar, tailRx, ...
    prefixRemovalLength, circularShiftLength)
% Function to generate a matrix to operate the reception process.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail at the receiver
%   window.
%   . prefixRemovalLength: Number of samples to be removed at the reception.
%   . cicularShiftLength: Number of samples to be shifted ("circularly").
%
% - Output:
%   . receptionMatrix: A matrix capable of performing all the process

removeRedundancyMatrix = remove_redundancy_matrix(numSubcar, tailRx, ...
    prefixRemovalLength);
overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx);
circularShiftMatrix = circular_shift_matrix(numSubcar, circularShiftLength);
transformMatrix = transform_matrix(numSubcar);
receptionMatrix = transformMatrix*circularShiftMatrix*overlapAddMatrix ...
    * windowRx*removeRedundancyMatrix;
end


function circularShiftMatrix = circular_shift_matrix(numSubcar, ...
    circularShiftLength)
% Function to generate matrix that operates the circular shift.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cicularShiftLength: Samples in the circular shift.
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
%   . prefixRemovalLength: Number of samples to be removed at the reception.

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
%   . addRedundancyMatrix: Matrix that adds the cyclic redundancy to the
%   system.

identityPrefix = eye(cpLength);
zerosPrefix = zeros(cpLength, (numSubcar-cpLength));
identitySuffix = eye(csLength);
zerosSuffix = zeros(csLength, (numSubcar-csLength));
identitySubcarriers = eye(numSubcar);
addRedundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


function overlapAddMatrix = overlap_and_add_matrix(numSubcar, tailRx)
% Function to generate matrix that operates the overlap and add in the
% samples.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
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
windowTxRC = diag(windowVector);
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
windowRxRC = diag(windowVector);
end


function transformMatrix = transform_matrix(numSubcar)

transformMatrix = zeros(numSubcar);
for row = 1:numSubcar
    for col = 1:numSubcar
        transformMatrix(row, col) = exp(-1j*(row-1)*2*pi*(col-1)/numSubcar);
    end
end
end


% EoF

