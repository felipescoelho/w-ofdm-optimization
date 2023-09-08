function addRedundancyMatrix = add_redundancy_matrix(dftLength, ...
    cpLength, csLength)
% ADD_REDUNDANCY_MATRIX         Generates the matrix that adds redundancy
%                               to the OFDM symbol.
%   addRedundancyMatrix = ADD_REDUNDANCY_MATRIX(dftLength, cpLength,
%       csLength) .
%
%   - Parameters:
%       . dftLength : Length of DFT in OFDM system.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%   - Returns:
%       . addRedundacyMatrix : Matrix responsible for adding redundancy to
%           the symbol.

identityPrefix = eye(cpLength);
zerosPrefix = zeros(cpLength, (dftLength-cpLength));
identitySuffix = eye(csLength);
zerosSuffix = zeros(csLength, (dftLength-csLength));
identitySubcarriers = eye(dftLength);
addRedundancyMatrix = [zerosPrefix identityPrefix; identitySubcarriers; ...
    identitySuffix zerosSuffix];
end


% EoF
