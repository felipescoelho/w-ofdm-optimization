function [costFunction, hessianMatrix] = cost_function_both(dftLength, ...
    tailRx, prefixRemovalLength, circularShiftLength, cpLength, ...
    csLength, channelArray)
% COST_FUNCTION_BOTH    Generates the cost function for fmincon and the
%                       Hessian matrix, to help the solver.
%   [costFunction, hessianMatrix] = COST_FUNCTION_BOTH(dftLength, tailRx,
%       prefixRemovalLength, circularShiftLength, cpLength, csLength,
%       channelArray)   .
%
%   - Parameters:
%       . dftLength : Length of DFT in OFDM system.
%       . tailRx : Length of rise and fall tails at the transmitter.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%       . circularShiftLength : Length of circular shift
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%       . channelArray : Array representing channel effects.
%   - Returns:
%       . costFunction : Cost function as a MATLAB function.
%       . hessianMatrix : Hessian matrix of objective function.

removeRedundancyMatrix = remove_redundancy_matrix(dftLength, tailRx, ...
    prefixRemovalLength);
overlapAddMatrix = overlap_and_add_matrix(dftLength, tailRx);
circularShiftMatrix = circular_shift_matrix(dftLength, circularShiftLength);
transformMatrix = dftmtx(dftLength);
addRedundancyMatrix = add_redundancy_matrix(dftLength, cpLength, csLength);
invTransformMatrix = conj(dftmtx(dftLength))/dftLength;
intercarrierInterference = channelArray(:, :, 1);
intersymbolInterference = sum(channelArray(:, :, 2:end), 3);
B = transformMatrix*circularShiftMatrix*overlapAddMatrix;

end


% EoF
