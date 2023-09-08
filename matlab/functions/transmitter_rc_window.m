function windowTxRC = transmitter_rc_window(dftLength, cpLength, ...
    csLength, tailTx)
% TRANSMITTER_RC_WINDOW     Generates diagonal matrix that operates the RC
%                           windowing at the transmitter.
%   windowTxRC = TRANSMITTER_RC_WINDOW(numSubcar, cpLength, csLength,
%       tailTx) .
%
%   - Parameters:
%       . dftLength : Length of DFT in OFDM system.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%       . tailTx : Length of rise and fall tails.
%   - Returns:
%       . windowTxRC : Diagonal matrix containing window's coefficients.

raisedCosineAxis = (-(tailTx+1)/2 + 1):1:((tailTx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailTx)).^2;
onesTransmitted = ones(1, dftLength+cpLength+csLength ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
windowTxRC = diag(windowVector);
end


% EoF
