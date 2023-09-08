function transmitterMatrix = tx_wofdm_matrix(windowTx, dftLength, ...
    cpLength, csLength)
% TX_WOFDM_MATRIX   Generates the matrix that concatenates all operations
%                   at the receiver for a given window.
%   transmitterMatrix = TX_WOFDM_MATRIX(windowTx, dftLength, cpLength,
%       csLength) .
%
%   - Parameters:
%       . windowTx : Diagonal matrix with window coefficients.
%       . dftLength : Length of DFT in OFDM system.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%   - Returns:
%       . transmitterMatrix : Matrix that concatenates all operations at
%           the transmitter.

invTransformMatrix = conj(dftmtx(dftLength))/dftLength;
addRedundancyMatrix = add_redundancy_matrix(dftLength, cpLength, csLength);
transmitterMatrix = windowTx*addRedundancyMatrix*invTransformMatrix;
end


% EoF
