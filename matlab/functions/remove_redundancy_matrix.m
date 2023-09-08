function removeRedundancyMatrix = remove_redundancy_matrix(passLength, ...
    tailRx, prefixRemovalLength)
% REMOVE_REDUNDANCY_MATRIX  Generates the matrix that removes redundancy
%                           from the received symbol.
%   removeRedundancyMatrix = REMOVE_REDUNDANCY_MATRIX(dftLength, tailRx,
%       prefixRemovalLength) .
%
%   - Parameters:
%       . passLength : Length of unchanged portion of information block.
%       . tailRx : Length of rise and fall tails at the receiver.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%   - Returns:
%       . removeRedundancyMatrix : Matrix responsible for removing
%           redundancy from the received symbol.

removeRedundancyMatrix = [zeros(passLength+tailRx, prefixRemovalLength) ...
    eye(passLength+tailRx)];
end


% EoF
