function overlapAddMatrix = overlap_and_add_matrix(dftLength, tailRx)
% OVERLAP_AND_ADD_MATRIX    Generates the matrix that operates the overlap
%                           and add at the receiver.
%   overlapAddMatrix = OVERLAP_AND_ADD_MATRIX(dftLength, tailRx) .
%
%   - Parameters:
%       . dftLength : Length of DFT in OFDM system.
%       . tailRx : Length of rise and fall tails.
%   - Returns:
%       . overlpaAddMatrix : Matrix that operates the overlap and add at
%           the receiver.

identityHalfTail = eye(tailRx/2);
identitySubcarriers = eye(dftLength-tailRx);
zerosHalfTail = zeros(tailRx/2);
zerosHalfTailSubcarriers = zeros(tailRx/2, dftLength-tailRx);
zerosTailSubcarriers = zeros(dftLength-tailRx, tailRx);
overlapAddMatrix = [zerosHalfTail identityHalfTail ...
    zerosHalfTailSubcarriers zerosHalfTail identityHalfTail; ...
    zerosTailSubcarriers identitySubcarriers zerosTailSubcarriers; ...
    identityHalfTail zerosHalfTail zerosHalfTailSubcarriers ...
    identityHalfTail zerosHalfTail];
end


% EoF
