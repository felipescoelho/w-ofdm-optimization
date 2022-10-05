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

% Definitions for run
MAX_REDUNDANCY = 35;
MIN_REDUNDANCY = 10;

parfor redundancyLength = MIN_REDUNDANCY:MAX_REDUNDANCY
    S = struct('chan', 'VehA200', 'numberSubcarriers', 256, ...
        'transformType', 'DFT');

    % General Definition:
    S.mu = ;  % CP 
    
    % wtx-OFDM:
    S.TypewOFDM = 'wtx';
    S.delta = 0;  % Samples in tail Rx window
    S.beta = 6;  % Samples in tail Tx window
    S.rho = S.beta;  % Length of CS
    S.gamma = cp_length;
    S.kappa = 0;  % Circular shift at receiver
    gen_opt_win(S)
    
    % wrx-OFDM:
    S.TypewOFDM = 'wrx';
    S.delta = 6;
    S.beta = 0;
    S.rho = S.delta/2;
    S.gamma = cp_length - S.rho;
    S.kappa = 0;
    gen_opt_win(S)
    
    % WOLA-OFDM:
    S.TypewOFDM = 'WOLA';
    S.delta = 6;
    S.beta = 6;
    S.rho = S.beta;
    S.gamma = cp_length - S.delta;
    S.kappa = S.delta/2;
    gen_opt_win(S)
    
    % CPW-OFDM:
    S.TypewOFDM = 'CPW';
    S.delta = 6;
    S.beta = 6;
    S.rho = S.beta + S.delta/2;
    S.gamma = cp_length - S.delta/2;
    S.kappa = 0;
    gen_opt_win(S)
    
    % CPwtx-OFDM:
    S.TypewOFDM = 'CPwtx';
    S.delta = 0;
    S.beta = 6;
    S.rho = 0;
    S.gamma = cp_length - S.beta;
    S.kappa = S.beta;
    Opt_win_gen(S)
    
    % CPwrx-OFDM:
    S.TypewOFDM = 'CPwrx';
    S.delta = 6;
    S.beta = 0;
    S.rho = 0;
    S.gamma = cp_length - S.delta;
    S.kappa = S.delta/2;
    Opt_win_gen(S)
    
end


function [cyclicPrefix, tailTx, tailRx] = select_variables( ...
    redundancyLength, seedValue, tryID)
% This function will select randomly values for the system.
%
% - Input:
%   . redundancyLength: Number of samples in cyclic prefix + cyclic suffix
%   in our system
%   . seedValue: a seed value for consistency, if needed
%
% - Output:
%   . cycliPrefix: Number of samples in cyclic prefix
%   . tailTx: Number of samples in rise and fall tail for the transmitter
%   . tailRx: Number of samples in rise and fall tail for the receiver

folderName = 'random_search';
if ~isdir(folderName)  %#ok
    mkdir(folderName)
end
valid = false;
while ~valid
    rng(seedValue)
    cyclicPrefix = randi(redundancyLength);
    
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

if ~isdir(folderName)  %#ok
    mkdir(folderName)
end
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
        windowVector = fmincon(funPreTx, initialValues, [], [], [], [], ...
            lowerBounds, upperBounds);
        optimizedTransmitterWindow = diag(windowVector);
        % Save results
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            cyclicPrefix, 'CP_', cyclicSuffix, 'CS');
        save([folderPath fileName], 'optimizedTransmitterWindow')
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
        windowVector = fmincon(funPreRx, initialValues, [], [], [], [], ...
            lowerBounds, upperBounds);
        optimizedReceiverWindow = diag(windowVector);
        % Save results
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            cyclicPrefix, 'CP_', cyclicSuffix, 'CS');
        save([folderPath fileName], 'optimizedReceiverWindow')
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
        % Transmitter
        funTxStep1 = optimization_function_tx(initialValuesRx, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        windowVector = fmincon(funTxStep1, initalValuesTx, [], [], [], ...
            [], lowerBoundsTx, upperBoundsTx);
        optimizedTransmitterStep1 = diag(windowVector);
        % Receiver
        funRxStep1 = optimization_function_rx(initialValuesTx, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        windowVector = fmincon(funRxStep1, initialValuesRx, [], [], [], ...
            [], lowerBoundsRx, upperBoundsRx);
        optimizedReceiverStep1 = diag(windowVector);
        % Second step
        % Transmitter
        funTxStep2 = optimization_function_tx(optimizedReceiverStep1, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        windowVector = fmincon(funTxStep2, initialValuesTx, [], [], [], ...
            [], lowerBoundsTx, upperBoundsTx);
        optimizedTransmitterStep2 = diag(windowVector);
        % Receiver
        funRxStep2 = optimization_function_rx(optimizedTransmitterStep1, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        windowVector = fmincon(funRxStep2, initialValuesRx, [], [], [], ...
            [], lowerBoundsRx, upperBoundsRx);
        optimizedReceiverStep2 = diag(windowVector);
        % Third step
        % Transmitter
        funTxStep3 = optimization_function_tx(optimizedReceiverStep2, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        windowVector = fmincon(funTxStep3, initialValuesTx, [], [], [], ...
            [], lowerBoundTx, upperBoundTx);
        optimizedTransmitterStep3 = diag(windowVector);
        % Receiver
        funRxStep3 = optimization_function_rx(optimizedTransmitterStep2, ...
            transformMatrix, circularShiftMatrix, overlapAdd, ...
            redundancyRemove, redundancyMatrix, invertTransform, ...
            intercarrierInterference, intersymbolInterference);
        windowVector = fmincon(funRxStep3, initialValuesRx, [], [], [], ...
            [], lowerBoundRx, upperBoundTx);
        optimizedReceiverStep3 = diag(windowVector);
        % Save results
        fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            cyclicPrefix, 'CP_', cyclicSuffix, 'CS');
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
for dataVectorAffected = 1:dataVectorsAffected
    for receivedSample = 1:numberReceived
        for transmittedSample = 1:numberTransmitted
            indexer = (dataVectorAffected-1)*channelNoise ...
                + transmittedSample - receivedSample;
            if (0 <= indexer) && (indexer <= channelOrder)
                interferenceArray(receivedSample, transmittedSample, ...
                    dataVectorAffected) = channelVector(indexer);
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

raisedCosineAxis = (-(tailRX+1)/2+1):1:((tailRx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailRx)).^2;
onesReceived = ones(1, numberSubcarriers-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
receiverWindow = diag(windowVector);
end


% EoF
