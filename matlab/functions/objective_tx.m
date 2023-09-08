function [hessianTx] = objective_tx(windowRx, dftLength, tailRx, ...
    tailTx, prefixRemovalLength, circularShiftLength, cpLength, ...
    csLength, channelArray)
% OBJECTIVE_TX  Generates the Hessian matrix of the objective function to
%               optimize the transmitter's window, minimizing the
%               interference power for a given channel array.
%   hessianTx = OBJECTIVE_TX(windowRx, dftLength, tailRx,
%       prefixRemovalLength, circularShiftLength, cpLength, csLength,
%       channelArray) .
%
%   - Parameters:
%       . windowRx : Diagonal matrix with receiver's window coefficients.
%       . dftLength : Length of DFT in OFDM system.
%       . tailRx : Length of rise and fall tails at the receiver.
%       . tailTx : Length of rise and fall tails at the transmitter.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%       . circularShiftLength : Length of circular shift.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%       . channelArray : Array representing channel effects.
%   - Returns:
%       . hessianTx : Hessian matrix of the objective function
%           (transmitter's optimization).

removeRedundancyMatrix = remove_redundancy_matrix(dftLength, tailRx, ...
    prefixRemovalLength);
overlapAddMatrix = overlap_and_add_matrix(dftLength, tailRx);
circularShiftMatrix = circular_shift_matrix(dftLength, circularShiftLength);
transformMatrix = dftmtx(dftLength);
addRedundancyMatrix = add_redundancy_matrix(dftLength, cpLength, csLength);
invTransformMatrix = conj(dftmtx(dftLength))/dftLength;
intercarrierInterference = channelArray(:, :, 1);
intersymbolInterference = sum(channelArray(:, :, 2:end), 3);
B = transformMatrix*circularShiftMatrix*overlapAddMatrix*windowRx ...
    * removeRedundancyMatrix*intercarrierInterference;
C = addRedundancyMatrix*invTransformMatrix;
% H1 is symmetric and we don't compute all terms to reduce computational
% time. If you want to compute all terms, please be advised of the
% numerical error caused by approximations. The matrix might not end up
% exactly symmetric or even real.
H1 = zeros(dftLength+cpLength+csLength);
for i = 1:dftLength+cpLength+csLength
    for j = 1:i
        for m = 1:dftLength
            indexMask = ones(1, dftLength);
            indexMask(m) = 0;
            indexMask = logical(indexMask);
            H1(i, j) = H1(i, j) + ctranspose(B(indexMask, i)) ...
                * conj(C(i, m))*B(indexMask, j)*C(j, m);
        end
    end
end
H1 = real(H1 + H1.' - diag(diag(H1)));
B2 = transformMatrix*circularShiftMatrix*overlapAddMatrix*windowRx ...
    * removeRedundancyMatrix*intersymbolInterference;
H2 = real(diag(diag(ctranspose(B2)*B2*C*ctranspose(C))));
H = H1+H2;
reduceVariableMatrix = reduce_variable_matrix( ...
    dftLength+cpLength+csLength-(2*tailTx), tailTx);
hessianTx = 2*transpose(reduceVariableMatrix)*H*reduceVariableMatrix;
end


% EoF
