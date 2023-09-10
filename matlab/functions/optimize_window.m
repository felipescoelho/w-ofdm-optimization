function optimize_window(typeOFDM, cpLength, csLength, dftLength, ...
    tailTx, tailRx, prefixRemovalLength, circularShiftLength)
% OPTIMIZE_WINDOW   Optimize window according to specifications.
%   OPTIMIZE_WINDOW(typeOFDM, cpLength, csLength, dftLength, tailTx,
%       tailRx, prefixRemovalLength, circularShiftLength) .
%
%   - Parameters:
%       . typeOFDM : Type of the w-OFDM system.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%       . dftLength : Length of DFT in OFDM system.
%       . tailTx : Length of rise and fall tails at the transmitter.
%       . tailRx : Length of rise and fall tails at the receiver.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%       . circularShiftLength: Length of circular shift.

global savePath channelPath

folderPath = [cd '/' savePath '/'];
load(channelPath, 'vehA200channel2')
avgChannel = mean(vehA200channel2, 1);
channelArray = channel_array(avgChannel, dftLength, cpLength, csLength, ...
    tailTx, tailRx, prefixRemovalLength);
fileName = strcat('optimal_win_', typeOFDM, '_VehA200_', ...
            num2str(cpLength), 'CP');
switch typeOFDM
    case {'wtx', 'CPwtx'}
        windowRx = diag(ones(dftLength+tailRx, 1));
        passLength = dftLength+cpLength+csLength-2*tailTx;
        hessianMatrix = objective_tx(windowRx, dftLength, tailRx, ...
            tailTx, prefixRemovalLength, circularShiftLength, cpLength, ...
            csLength, channelArray);
        [A, b, Aeq, beq] = tx_constraints(tailTx);
        f = zeros(tailTx+1, 1);
        [variableVector, estimatedPower, ~, output] = quadprog( ...
            hessianMatrix, f, A, b, Aeq, beq);
        optimizedWindow = recover_window(variableVector, passLength);
        save([folderPath fileName] , 'optimizedWindow')
        fprintf('Estimated Power: %.4e \n', estimatedPower)
        fprintf('Number of iterations: %u \n', output.iterations)
    case {'wrx', 'CPwrx'}
        windowTx = diag(ones(dftLength+cpLength+csLength, 1));
        passLength = dftLength-tailRx;
        hessianMatrix = objective_rx(windowTx, dftLength, tailRx, ...
            prefixRemovalLength, circularShiftLength, cpLength, ...
            csLength, channelArray);
        [A, b, Aeq, beq] = rx_constraints(tailRx);
        f = zeros(tailRx+1, 1);
        [variableVector, estimatedPower, ~, output] = quadprog( ...
            hessianMatrix, f, A, b, Aeq, beq);
        optimizedWindow = recover_window(variableVector, passLength);
        save([folderPath fileName] , 'optimizedWindow')
        fprintf('Estimated Power: %.4e \n', estimatedPower)
        fprintf('Number of iterations: %u \n', output.iterations)
%     case {'WOLA', 'CPW'}
%         passLengthTx = dftLength+cpLength+csLength-2*tailTx;
%         passLengthRx = dftLength-tailRx;
%         [A, b, Aeq, beq] = both_constraints();
        
end
end


% EoF
