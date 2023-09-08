function windowRxRC = receiver_rc_window(dftLength, tailRx)
% RECEIVER_RC_WINDOW    Generates diagonal matrix that operates the RC
%                       windowing at the receiver.
%   windowRxRC = RECEIVER_RC_WINDOW(dftLength, tailRx) .
%
%   - Parameters:
%       . dftLength : Length of DFT in OFDM system.
%       . tailRx : Length of rise and fall tails.
%   - Returns:
%       . windowRxRC : Diagonal matrix containing window's coefficients.

raisedCosineAxis = (-(tailRx+1)/2+1):1:((tailRx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailRx)).^2;
onesReceived = ones(1, dftLength-tailRx);
windowVector = [raisedCosine onesReceived fliplr(raisedCosine)];
windowRxRC = diag(windowVector);
end


% EoF
