function [A, b, Aeq, beq] = rx_constraints(tailRx)
% RX_CONSTRAINTS    Generates the constraint matrices for the receiver
%                   optimization problem.
%   [A, b, Aeq, beq] = RX_CONSTRAINTS(tailRx)   The generated constraints
%       are desingned for the problem.
%
%           A*x <= b
%           Aeq*x = beq
%
%   - Parameters:
%       . tailRx : Length of the rise and fall tails for the receiver.
%   - Returns:
%       . A : Matrix that defines variables for linear inequality
%           constraints.
%       . b : Vector with values for linear inequality constraints.
%       . Aeq : Matrix that defines variables for linear equality
%           constraints.
%       . beq : Vector with values for linear equality constraints.

A = [zeros(tailRx, 1) eye(tailRx)];
b = ones(tailRx, 1);
Aeq = [1 zeros(1, tailRx);
    zeros(tailRx, 1) (eye(tailRx) - fliplr(eye(tailRx)))];
beq = [1; zeros(tailRx, 1)];
end


% EoF
