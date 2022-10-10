% main_window_optimization.m
%   This script will perform the optimization of the windowed-OFDM for all
%   size tested w-OFDM systems wtx, wrx, WOLA, CPW, CPwtx, and CPwrx.
%
% I have three variables to define...
%   - mu: the length of my cyclic prefix
%   - beta: the length of my rise and fall tail for Tx window (only wtx,
%   WOLA, CPW, and CPwtx)
%   - delta: the length of my rise and fall tail for Rx window (only wrx,
%   CPW, WOLA, and CPwrx)
%
% All other variables are derived from these.
%   - rho: length of the cyclic suffix
%       . wtx, WOLA: rho = beta
%       . CPW: rho = delta/2 + beta
%       . wrx: rho = delta/2
%   - gamma: number of samples removed at reception
%       . wtx: gamma = mu
%       . wrx, CPW: gamma = mu - delta/2
%       . WOLA, CPwrx: gamma = mu - delta
%       . CPwtx: gamma = mu - delta
%   - kappa: number of samples in circular shift
%       . WOLA, CPwrx: delta/2
%       . CPwtx: beta
%
% I want to analyse based on the length of redundancy.
%   - redundancy = mu + rho
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
% Definitions
%
% From paper:
% wtx: mu=32, rho=8, N=256, beta=8, delta=0, gamma=32, kappa=0
% wrx: mu=32, rho=5, N=256, beta=0, delta=10, gamma=27, kappa=0
% WOLA: mu=32, rho=8, N=256, beta=8, delta=10, gamma=22, kappa=5
% CPW: mu=32, rho=13, N=256, beta=8, delta=10, gamma=27, kappa=0
% CPwtx: mu=32, rho=0, N=256, beta=8, delta=0, gamma=24, kappa=8
% CPwrx: mu=32, rho=0, N=256, beta=0, delta=10, gamma=22, kappa=5
%

folderName = 'optimized_windows';
if ~isfolder(folderName)
    mkdir(folderName)
end
settingsFileName = 'settingsData.mat';
fprintf('Generating settingsData structure array.\n')
settingsData = struct();
settingsData.generalSettings = struct('numberSubcarriers', 256, ...
    'bitsPerSubcarrier', 4, 'cyclicPrefix', 32);
settingsData.wtx = struct('cyclicSuffix', 8, 'tailTx', 8, 'tailRx', 0, ...
    'prefixRemoval', 32, 'circularShift', 0);
settingsData.wrx = struct('cyclicSuffix', 5, 'tailTx', 0, 'tailRx', 10, ...
    'prefixRemoval', 27, 'circularShift', 0);
settingsData.WOLA = struct('cyclicSuffix', 8, 'tailTx', 8, 'tailRx', ...
    10, 'prefixRemoval', 22, 'circularShift', 5);
settingsData.CPW = struct('cyclicSuffix', 13, 'tailTx', 8, 'tailRx', ...
    10, 'prefixRemoval', 27, 'circularShift', 0);
settingsData.CPwtx = struct('cyclicSuffix', 0, 'tailTx', 8, 'tailRx', ...
    0, 'prefixRemoval', 24, 'circularShift', 8);
settingsData.CPwrx = struct('cyclicSuffix', 0, 'tailTx', 0, 'tailRx', ...
    10, 'prefixRemoval', 22, 'circularShift', 5);
fprintf('Saving settingsData structure array in %s\n', settingsFileName)
save(settingsFileName, 'settingsData')
fprintf('Successfully saved settingsData.\n')
fprintf('------------------------------------\n')
fprintf('Starting run ...\n\n')
fields = fieldnames(settingsData);
clear settingsData
for fieldIndex = 1:length(fields)
    fieldName = fields{fieldIndex};
    loadingStructure = load(settingsFileName);
    settingsData = loadingStructure.settingsData;
    switch fieldName
        case {'wtx', 'wrx', 'WOLA', 'CPW', 'CPwtx', 'CPwrx'}
            fprintf('Working on %s-OFDM ...\n', fieldName)
            numSubcar = settingsData.generalSettings.numberSubcarriers;
            cpLength = settingsData.generalSettings.cyclicPrefix;
            csLength = settingsData.(fieldName).cyclicSuffix;
            tailTx = settingsData.(fieldName).tailTx;
            tailRx = settingsData.(fieldName).tailRx;
            prefixRemoval = settingsData.(fieldName).prefixRemoval;
            circularShift = settingsData.(fieldName).circularShift;
            fprintf('--------------------------------------\n')
            fprintf('Settings:\n')
            fprintf('N = %u (Number of subcarriers)\n', numSubcar)
            fprintf('mu = %u (CP length)\n', cpLength)
            fprintf('rho = %u (CS length)\n', csLength)
            fprintf('beta = %u (Tx window tail length)\n', tailTx)
            fprintf('delta = %u (Rx window tail length)\n', tailRx)
            fprintf('gamma = %u (Number of samples removed', prefixRemoval)
            fprintf(' in receiver)\n')
            fprintf('kappa = %u (Number of samples in ', circularShift)
            fprintf('circular shift)\n')
            fprintf('--------------------------------------\n')
            optimize_window(fieldName, cpLength, csLength, numSubcar, ...
                tailTx, tailRx, prefixRemoval, circularShift, folderName)
        otherwise
            continue
    end
end


function optimize_window(typeOFDM, cyclicPrefix, cyclicSuffix, ...
    numberSubcarriers, tailTx, tailRx, prefixRemoval, circularShift, ...
    folderName)
% Function to generate the optimized windows.
%
% - Input:
%   . typeOFDM: Type of OFDM system
%   . cyclicPrefix: Number of samples in cyclic prefix (mu)
%   . cyclicSuffix: Number of samples in cyclic suffix (rho)
%   . numberSubcarriers: Number of subcarriers in system (N)
%   . tailTx: Number of samples in rise and fall tail at transmitter (beta)
%   . tailRx: Number of samples in rise and fall at receiver (delta)
%   . prefixRemoval: Number of samples (from cyclic prefix) to be removed
%   at the receiver. (gamma)
%   . circularShift: Number of samples in circular shift, for WOLA, CPwtx,
%   and CPwrx. (kappa)
%   . folderName: Name of the folder to save optimized windows.

folderPath = [cd '/' folderName '/'];
load('./channels/vehA200channel2.mat', 'vehA200channel2')
averageChannel = mean(vehA200channel2, 1);
interferenceArray = array_ici_isi(averageChannel, numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix, tailTx, tailRx, prefixRemoval);
intercarrierInterference = interferenceArray(:, :, 1);
intersymbolInterference = sum(interferenceArray(:, :, 2:end), 3);
switch typeOFDM
    case {'wtx', 'CPwtx'}
        % Transmitter process
        redundancyMatrix = add_redundancy_matrix(numberSubcarriers, ...
            cyclicPrefix, cyclicSuffix);
        invertTransform = dftmtx(numberSubcarriers)'/numberSubcarriers;
        % Receiver process
        receiverMatrix = rx_wofdm_matrix(numberSubcarriers, tailRx, ...
            prefixRemoval, circularShift);
        % Optimization process
        fprintf('Starting optimization process\n')
        funPreTx = @(x) norm(receiverMatrix*intercarrierInterference ...
            * (diag(x)*redundancyMatrix*invertTransform) ...
            - diag(diag(receiverMatrix*intercarrierInterference ...
            * (diag(x)*redundancyMatrix*invertTransform))), 'fro') ...
            + norm(receiverMatrix*intersymbolInterference ...
            * (diag(x)*redundancyMatrix*invertTransform));
        initialValues = diag(transmitter_rc_window(numberSubcarriers, ...
            cyclicPrefix, cyclicSuffix, tailTx));
        lowerBounds = [zeros(tailTx, 1); ...
            ones(numberSubcarriers+cyclicPrefix-tailTx, 1); ...
            zeros(tailTx, 1)];
        upperBounds = [ones(tailTx, 1); ...
            ones(numberSubcarriers+cyclicPrefix-tailTx, 1); ...
            ones(tailTx, 1)];
        [windowVector, result] = fmincon(funPreTx, initialValues, [], ...
            [], [], [], lowerBounds, upperBounds);
        optimizedTransmitterWindow = diag(windowVector);
        fprintf('Finished optimization process\n')
        fprintf('Result for %s: %.4f W.\n', typeOFDM, result)
        % Save results
        fprintf('Saving results...\n')
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cyclicPrefix), 'CP');
        save([folderPath fileName], 'optimizedTransmitterWindow')
        fprintf('Optimized windows were successfully saved into file.\n')
        fprintf('\n')
    case {'wrx', 'CPwrx'}
        % Transmitter process
        transmitterMatrix = tx_wofdm_matrix(numberSubcarriers, ...
            cyclicPrefix, cyclicSuffix, tailTx);
        % Receiver process
        transformMatrix = dftmtx(numberSubcarriers);
        circularShiftMatrix = circular_shift_matrix(numberSubcarriers, ...
            circularShift);
        overlapAdd = overlap_and_add_matrix(numberSubcarriers, tailRx);
        redundancyRemove = remove_redundancy_matrix(numberSubcarriers, ...
            tailRx, prefixRemoval);
        % Optimization process
        fprintf('Starting optimization process\n')
        funPreRx = @(x) norm((transformMatrix*circularShiftMatrix ...
            * overlapAdd*diag(x)*redundancyRemove) ...
            * intercarrierInterference*transmitterMatrix ...
            - diag(diag((transformMatrix*circularShiftMatrix ...
            * overlapAdd*diag(x)*redundancyRemove) ...
            * intercarrierInterference*transmitterMatrix)), 'fro') ...
            + norm((transformMatrix*circularShiftMatrix*overlapAdd ...
            * diag(x)*redundancyRemove)*intersymbolInterference ...
            * transmitterMatrix, 'fro');
        initialValues = diag(receiver_rc_window(numberSubcarriers, ...
            tailRx));
        lowerBounds = [zeros(tailRx, 1); ...
            ones(numberSubcarriers-tailRx, 1); zeros(tailRx, 1)];
        upperBounds = [ones(tailRx, 1); ...
            ones(numberSubcarriers-tailRx, 1); ones(tailRx, 1)];
        [windowVector, result] = fmincon(funPreRx, initialValues, [], ...
            [], [], [], lowerBounds, upperBounds);
        optimizedReceiverWindow = diag(windowVector);
        fprintf('Finished optimization process\n')
        fprintf('Result for %s: %.4f\n.', typeOFDM, result)
        % Save results
        fprintf('Saving results...\n')
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cyclicPrefix), 'CP');
        save([folderPath fileName], 'optimizedReceiverWindow')
        fprintf('Optimized windows were successfully saved into file.\n')
        fprintf('\n')
    case {'WOLA', 'CPW'}
        % Transmitter process
        redundancyMatrix = add_redundancy_matrix(numberSubcarriers, ...
            cyclicPrefix, cyclicSuffix);
        invertTransform = dftmtx(numberSubcarriers)'/numberSubcarriers;
        % Receiver process
        transformMatrix = dftmtx(numberSubcarriers);
        circularShiftMatrix = circular_shift_matrix(numberSubcarriers, ...
            circularShift);
        overlapAdd = overlap_and_add_matrix(numberSubcarriers, tailRx);
        redundancyRemove = remove_redundancy_matrix(numberSubcarriers, ...
            tailRx, prefixRemoval);
        % Optimization process
        % Initialization
        fprintf('Starting optimization process\n')
        initialValuesTx = diag(transmitter_rc_window(numberSubcarriers, ...
            cyclicPrefix, cyclicSuffix, tailTx));
        initialValuesRx = diag(receiver_rc_window(numberSubcarriers, ...
            tailRx));
        lowerBoundsTx = [zeros(tailTx, 1); ...
            ones(numberSubcarriers+cyclicPrefix-tailTx, 1); ...
            zeros(tailTx, 1)];
        upperBoundsTx = [ones(tailTx, 1); ...
            ones(numberSubcarriers+cyclicPrefix-tailTx, 1); ...
            ones(tailTx, 1)];
        lowerBoundsRx = [zeros(tailRx, 1); ...
            ones(numberSubcarriers-tailRx, 1); zeros(tailRx, 1)];
        upperBoundsRx = [ones(tailRx, 1); ...
            ones(numberSubcarriers-tailRx, 1); ones(tailRx, 1)];
        % First step
        fprintf('Starting step 1:\n')
        % Transmitter
        fprintf('Starting optimization for transmitter...\n')
        funTxStep1 = optimization_function_tx(diag(initialValuesRx), ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        [windowVector, result] = fmincon(funTxStep1, initialValuesTx, ...
            [], [], [], [], lowerBoundsTx, upperBoundsTx);
        optimizedTransmitterStep1 = diag(windowVector);
        fprintf('Finished optimization.\n')
        fprintf('Result for %s: %.4f\n', typeOFDM, result)
        % Receiver
        fprintf('Starting optimization for receiver...\n')
        funRxStep1 = optimization_function_rx(diag(initialValuesTx), ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        [windowVector, result] = fmincon(funRxStep1, initialValuesRx, ...
            [], [], [], [], lowerBoundsRx, upperBoundsRx);
        optimizedReceiverStep1 = diag(windowVector);
        fprintf('Finished optimization.\n')
        fprintf('Result for %s: %.4f\n', typeOFDM, result)
        % Second step
        fprintf('Starting step 2:\n')
        % Transmitter
        fprintf('Starting optimization for transmitter...\n')
        funTxStep2 = optimization_function_tx(optimizedReceiverStep1, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        [windowVector, result] = fmincon(funTxStep2, initialValuesTx, ...
            [], [], [], [], lowerBoundsTx, upperBoundsTx);
        optimizedTransmitterStep2 = diag(windowVector);
        fprintf('Finished optimization.\n')
        fprintf('Result for %s: %.4f\n', typeOFDM, result)
        % Receiver
        fprintf('Starting optimization for receiver...\n')
        funRxStep2 = optimization_function_rx(optimizedTransmitterStep1, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        [windowVector, result] = fmincon(funRxStep2, initialValuesRx, ...
            [], [], [], [], lowerBoundsRx, upperBoundsRx);
        optimizedReceiverStep2 = diag(windowVector);
        fprintf('Finished optimization.\n')
        fprintf('Result for %s: %.4f\n', typeOFDM, result)
        % Third step
        fprintf('Starting step 3:\n')
        % Transmitter
        fprintf('Starting optimization for transmitter...\n')
        funTxStep3 = optimization_function_tx(optimizedReceiverStep2, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        [windowVector, result] = fmincon(funTxStep3, initialValuesTx, ...
            [], [], [], [], lowerBoundsTx, upperBoundsTx);
        optimizedTransmitterStep3 = diag(windowVector);
        fprintf('Finished optiization.\n')
        fprintf('Result for %s: %.4f\n', typeOFDM, result)
        % Receiver
        fprintf('Starting optimization for receiver...\n')
        funRxStep3 = optimization_function_rx(optimizedTransmitterStep2, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        [windowVector, result] = fmincon(funRxStep3, initialValuesRx, ...
            [], [], [], [], lowerBoundsRx, upperBoundsTx);
        optimizedReceiverStep3 = diag(windowVector);
        fprintf('Finished optimization.\n')
        fprintf('Result for %s: %.4f\n', typeOFDM, result)
        % Save results
        fprintf('Saving results...\n')
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cyclicPrefix), 'CP');
        save([folderPath fileName], 'optimizedTransmitterStep1', ...
            'optimizedReceiverStep1', 'optimizedTransmitterStep2', ...
            'optimizedReceiverStep2','optimizedTransmitterStep3', ...
            'optimizedReceiverStep3')
        fprintf('Optimized windows were successfully saved into file.')
        fprintf('\n')
end

end


function optimizationFunction = optimization_function_tx(windowRx, ...
    transformMatrix, circularShiftMatrix, overlapAdd, redundancyRemove, ...
    redundancyMatrix, invertTransform, intercarrierInterference, ...
    intersymbolInterference)
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

optimizationFunction = @(x) norm((transformMatrix*circularShiftMatrix ...
    * overlapAdd*windowRx*redundancyRemove) ... 
    * intercarrierInterference*(diag(x)*redundancyMatrix ...
    * invertTransform) - diag(diag((transformMatrix*circularShiftMatrix ...
    * overlapAdd*windowRx*redundancyRemove)*intercarrierInterference ...
    * (diag(x)*redundancyMatrix*invertTransform))), 'fro') ...
    + norm((transformMatrix*circularShiftMatrix*overlapAdd ...
    * windowRx*redundancyRemove)*intersymbolInterference ...
    * (diag(x)*redundancyMatrix*invertTransform));
end


function optimizationFunction = optimization_function_rx(windowTx, ...
    transformMatrix, circularShiftMatrix, overlapAdd, redundancyRemove, ...
    redundancyMatrix, invertTransform, intercarrierInterference, ...
    intersymbolInterference)
% Funtion to help writing tje optimization functions for Rx
%
% - Input:
%   . windowTx: Window used at the receiver
%   . transformMatrix: DFT matrix
%   . circularShiftMatrix: Matrix to perform the circular shift
%   . overlapAdd: Matrix to perform the overlap and add operation
%   . redundancyRemove: Matrix to remove redundancy at recepetion
%   . redundancyMatrix: Matrix to add redundancy at transmission
%   invertTransform: IDFT matrix
%   . intercarrierInterference: Matrix with intercarrier interference for
%   the channel
%   . intersymbolInterference: Matrix with intersymbol interference for the
%   channel.
%
% - Output:
%   . optimizationFunction: Function to perform the optimization process
%   using MATLAB's fmincon function.

optimizationFunction = @(x) norm((transformMatrix*circularShiftMatrix ...
    * overlapAdd*diag(x)*redundancyRemove)*intercarrierInterference ...
    * (windowTx*redundancyMatrix*invertTransform) - diag(diag(( ...
    transformMatrix*circularShiftMatrix*overlapAdd*diag(x) ...
    * redundancyRemove)*intercarrierInterference*(windowTx ...
    * redundancyMatrix*invertTransform))), 'fro') + norm(( ...
    transformMatrix*circularShiftMatrix*overlapAdd*diag(x) ...
    * redundancyRemove)*intersymbolInterference*(windowTx ...
    * redundancyMatrix*invertTransform), 'fro');
end


function [interferenceArray] = array_ici_isi(channelVector, ...
    numberSubcarriers, cyclicPrefix, cyclicSuffix, tailTx, tailRx, ...
    prefixRemoval)
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

channelOrder = length(channelVector) - 1;
numberReceived = numberSubcarriers+tailRx+prefixRemoval;
numberTransmitted = numberSubcarriers+cyclicPrefix+cyclicSuffix;
channelNoise = numberTransmitted-tailTx;
dataVectorsAffected = ceil((channelOrder+tailTx)/numberReceived) + 1;
interferenceArray = zeros(numberReceived, numberTransmitted, ...
    dataVectorsAffected);
for dataVectorAffected = 0:dataVectorsAffected-1
    for receivedSample = 0:numberReceived-1
        for transmittedSample = 0:numberTransmitted-1
            indexer = dataVectorAffected*channelNoise ...
                + transmittedSample - receivedSample;
            if (0 <= indexer) && (indexer <= channelOrder)
                interferenceArray(receivedSample+1, ...
                    transmittedSample+1, dataVectorAffected+1) ...
                    = channelVector(indexer+1);
            end
        end
    end
end

end


function transmissionMatrix = tx_wofdm_matrix(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix, tailTx)
% Function to generate a matrix to operate the transmission process.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cyclicPrefix: Number of samples in cyclic prefix.
%   . cuclicSuffix: Number of samples in cyclic suffix.
%   . tailTx : Number of samples in transmitter window tails.
%
% - Output:
%   . transmissionMatrix: Matrix that operates the transmission process.
%   . transmitterWindow: Matrix containing the transmission window.
%   . redundancyMatrix: Matrix that operates the redundancy addition.
%   . invertTransform: Matrix that operates the invert discrete Fourier
%   transform.

transformMatrix = dftmtx(numberSubcarriers);
invertTransform = transformMatrix'/numberSubcarriers;
redundancyMatrix = add_redundancy_matrix(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix);
transmitterWindow = transmitter_rc_window(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix, tailTx);
transmissionMatrix = transmitterWindow*redundancyMatrix*invertTransform;
end


function receptionMatrix = rx_wofdm_matrix(numberSubcarriers, tailRx, ...
    prefixRemoval, circularShift)
% Function to generate a matrix to operate the reception process.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail at the receiver
%   window.
%   . prefixRemoval: Number of samples to be removed at the reception.
%   . cicularShift: Number of samples to be shifted ("circularly").
%
% - Output:
%   . receptionMatrix: A matrix capable of performing all the process

redundancyRemove = remove_redundancy_matrix(numberSubcarriers, tailRx, ...
    prefixRemoval);
receiverWindow = receiver_rc_window(numberSubcarriers, tailRx);
overlapAdd = overlap_and_add_matrix(numberSubcarriers, tailRx);
circularShiftMatrix = circular_shift_matrix(numberSubcarriers, ...
    circularShift);
transformMatrix = dftmtx(numberSubcarriers);
receptionMatrix = transformMatrix*circularShiftMatrix*overlapAdd ...
    *receiverWindow*redundancyRemove;
end


function circularShiftMatrix = circular_shift_matrix(numberSubcarriers, ...
    circularShift)
% Function to generate matrix that operates the circular shift.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cicularShift: Samples in the circular shift.
%
% - Output:
%   . circularShiftMatrix: A matrix capable of performing the circular
%   shift operation.

identityCircularShift = eye(circularShift);
identitySubcarriersMinusCircularShift = eye(numberSubcarriers- ...
    circularShift);
zerosSubcarriersMinusCircularShift = zeros(circularShift, ...
    numberSubcarriers-circularShift);
circularShiftMatrix = [zerosSubcarriersMinusCircularShift' ...
    identitySubcarriersMinusCircularShift; identityCircularShift ...
    zerosSubcarriersMinusCircularShift];
end


function redundancyRemove = remove_redundancy_matrix(numberSubcarriers, ...
    tailRx, prefixRemoval)
% Function to generate matrix to remove redundancy.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail at the receiver
%   window.
%   . prefixRemoval: Number of samples to be removed at the reception.

redundancyRemove = [zeros(numberSubcarriers+tailRx, prefixRemoval) ...
    eye(numberSubcarriers+tailRx)];
end


function redundancyMatrix = add_redundancy_matrix(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix)
% Function to generate the matrix that adds redundancy.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . cyclicPrefix: Number of samples in cyclic prefix.
%   . cyclicSuffix: Number of samples in cyclic suffix.
%
% - Output:
%   . redundancyMatrix: Matrix that adds the cyclic redundancy to the
%   system.

identityPrefix = eye(cyclicPrefix);
zerosPrefix = zeros(cyclicPrefix, (numberSubcarriers-cyclicPrefix));
identitySuffix = eye(cyclicSuffix);
zerosSuffix = zeros(cyclicSuffix, (numberSubcarriers-cyclicSuffix));
identitySubcarriers = eye(numberSubcarriers);
redundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


function overlapAdd = overlap_and_add_matrix(numberSubcarriers, tailRx)
% Function to generate matrix that operates the overlap and add in the
% samples.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tails for the receiver
%   window.

identityHalfTail = eye(tailRx/2);
identitySubcarriers = eye(numberSubcarriers-tailRx);
zerosHalfTail = zeros(tailRx/2);
zerosHalfTailSubcarriers = zeros(tailRx/2, numberSubcarriers-tailRx);
zerosTailSubcarriers = zeros(numberSubcarriers-tailRx, tailRx);
overlapAdd = [zerosHalfTail identityHalfTail zerosHalfTailSubcarriers ...
    zerosHalfTail identityHalfTail; zerosTailSubcarriers ...
    identitySubcarriers zerosTailSubcarriers; ...
    identityHalfTail zerosHalfTail zerosHalfTailSubcarriers ...
    identityHalfTail zerosHalfTail];
end


function transmitterWindow = transmitter_rc_window(numberSubcarriers, ...
    cyclicPrefix, cyclicSuffix, tailTx)
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

raisedCosineAxis = (-(tailTx+1)/2 + 1):1:((tailTx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailTx)).^2;
onesTransmitted = ones(1, numberSubcarriers+cyclicPrefix+cyclicSuffix ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
transmitterWindow = diag(windowVector);
end


function receiverWindow = receiver_rc_window(numberSubcarriers, tailRx)
% Function to generate the Raised Cosine receiver window.
%
% - Input:
%   . numberSubcarriers: Number of subcarriers in the system.
%   . tailRx: Number of samples in rise and fall tail for the receiver
%   window.

raisedCosineAxis = (-(tailRx+1)/2+1):1:((tailRx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailRx)).^2;
onesReceived = ones(1, numberSubcarriers-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
receiverWindow = diag(windowVector);
end


% EoF
