function receiverMatrix = rx_wofdm_matrix(windowRx, dftLength, tailRx, ...
    prefixRemovalLength, circularShiftLength)
% RX_WOFDM_MATRIX   Generates the matrix that concatenates all operations
%                   at the receiver for a given window.
%   receiverMatrix = RC_WOFDM_MATRIX(windowRx, dftLength, tailRx,
%       prefixRemovalLength, circularShiftLength) .
%
%   - Parameters:
%       . windowRx : Diagonal matrix with window coefficients.
%       . dftLength : Length of DFT in OFDM system.
%       . tailRx : Length of rise and fall tails at the receiver.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%       . circularShiftLength : Length of circular shift.
%   - Returns:
%       . receiverMatrix : Matrix that concatenates all operations at the
%           receiver.

removeRedundancyMatrix = remove_redundancy_matrix(dftLength, tailRx, ...
    prefixRemovalLength);
overlapAddMatrix = overlap_and_add_matrix(dftLength, tailRx);
circularShiftMatrix = circular_shift_matrix(dftLength, circularShiftLength);
transformMatrix = dftmtx(dftLength);
receiverMatrix = transformMatrix*circularShiftMatrix*overlapAddMatrix ...
    * windowRx*removeRedundancyMatrix;
end


% EoF
