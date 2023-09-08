function [A, b, Aeq, beq] = tx_constraints(tailTx)
% TX_CONSTRAINTS    Generates the constraint matrices for the transmitter
%                   optimization problem.
%   [A, b, Aeq, beq] = TX_CONSTRAINTS(tailTx) .
%
%   - Parameters:
%       . tailTx : Length of the rise and fall tails for the transmitter.
%   - Returns:
%       . A : Matrix that defines variables for linear inequality
%           constraints.
%       . b : Vector with values for linear ineqaulity constraints.
%       . Aeq : Matrix that defines variables for linear equality
%           constriants.
%       . beq : Vector with values for linear equality constraints.

A = [zeros(tailTx, 1) eye(tailTx)];
b = ones(tailTx, 1);
Aeq = [1 zeros(1, tailTx)];
beq = 1;
end


% EoF
