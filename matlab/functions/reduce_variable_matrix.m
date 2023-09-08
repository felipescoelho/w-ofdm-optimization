function reduceVariableMatrix = reduce_variable_matrix(passLength, ...
    tailLength)
% REDUCE_VARIABLE_MATRIX     Generates a matrix to remove redundancy
%                            from the optmization variable.
%   reduceVariableMatrix = REDUCE_VARIABLE_MATRIX(dftLength, tailLength) .
%
%   - Parameters:
%       . passLength : Length of unchanged portion of infomation block.
%       . tailLength : Length of rise and fall tails.
%   - Returns:
%       . reduceVariableMatrix : Matrix to reduce variable length for the
%           optimization process.

reduceVariableMatrix = [zeros(tailLength, 1) fliplr(eye(tailLength)); ...
    ones(passLength, 1) zeros(passLength, tailLength); ...
    zeros(tailLength, 1) eye(tailLength)];
end


% EoF
