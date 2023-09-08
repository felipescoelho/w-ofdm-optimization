function hessianRx = objective_rx(windowTx, dftLength, tailRx, ...
    prefixRemovalLength, circularShiftLength, cpLength, csLength, ...
    channelArray)
% OBJECTIVE_RX  Generates the Hessian matrix of the objective function to
%               optimize the receiver's window, minimizing the interference
%               power for a given channel array.
%   hessianRx = OBJECTIVE_RX(windowTx, dftLength, tailRx,
%       prefixRemovalLength, circularShiftLength, cpLength, csLength,
%       channelArray) .
%
%   - Parameters:
%       . windowTx : Diagonal matrix with transmitter's window
%           coefficients.
%       . dftLength : Length of DFT in OFDM system.
%       . tailRx : Length of rise and fall tails at the receiver.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%       . circularShiftLength : Length of circular shift.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%       . channelArray : Array representing channel effects.
%   - Returns:
%       . hessianRx : Hessian matrix of the objective function (receiver's
%           optimization).

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
C = removeRedundancyMatrix*intercarrierInterference*windowTx ...
    *addRedundancyMatrix*invTransformMatrix;
% H1 is symmetric and we don't compute all terms to reduce computational
% time. If you want to compute all terms, please be advised of the
% numerical error caused by approximations. The matrix might not end up
% exactly symmetric or even real.
H1 = zeros(dftLength+tailRx);
for i = 1:dftLength+tailRx
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
C2 = removeRedundancyMatrix*intersymbolInterference*windowTx ...
    * addRedundancyMatrix*invTransformMatrix;
H2 = real(diag(diag(ctranspose(B)*B*C2*ctranspose(C2))));
H = H1 + H2;
reduceVariableMatrix = reduce_variable_matrix(dftLength, tailRx);
hessianRx = 2*(transpose(reduceVariableMatrix)*H*reduceVariableMatrix);
end


% EoF
