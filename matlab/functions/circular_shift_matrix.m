function circularShiftMatrix = circular_shift_matrix(dftLength, ...
    circularShiftLength)
% CIRCULAR_SHIFT_MATRIX     Generates the matrix responsible for the
%                           circular shift at the receiver.
%   circuShiftMatrix = CIRCULAR_SHIFT_MATRIX(dftLength,
%       circularShiftLength) .
%
%   - Parameters:
%       . dftLength : Length of DFT in OFDM system.
%       . circularShiftLength : Length of circular shift.
%   - Retruns:
%       . circularShiftMatrix : Matrix responsible for performing the
%           circular shift at the receiver.

identityCircularShift = eye(circularShiftLength);
identitySubcarriersMinusCircularShift = eye(dftLength-circularShiftLength);
zerosSubcarriersMinusCircularShift = zeros(circularShiftLength, ...
    dftLength-circularShiftLength);
circularShiftMatrix = [zerosSubcarriersMinusCircularShift.' ...
    identitySubcarriersMinusCircularShift; identityCircularShift ...
    zerosSubcarriersMinusCircularShift];
end


% EoF
