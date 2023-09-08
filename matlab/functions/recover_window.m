function window = recover_window(optimizedVariable, passLength)
% RECOVER_WINDOW    Recovers the window matrix from the vector after the
%                   optimization process.
%   window = RECOVER_WINDOW(optimizedVariable) .
%
%   - Parameters:
%       . optimizedVariable : Vector containing resulting values from
%           optimization.
%       . passLength : Length of unchanged portion of information block.
%   - Returns:
%       . window : Diagonal matrix with optimizaed window coefficients.

tailLength = length(optimizedVariable) - 1;
recoverWindowMatrix = [zeros(tailLength, 1) fliplr(eye(tailLength));
    ones(passLength, 1) zeros(passLength, tailLength);
    zeros(tailLength, 1) eye(tailLength)];
windowVector = recoverWindowMatrix * optimizedVariable;
window = diag(windowVector);
end


% EoF
