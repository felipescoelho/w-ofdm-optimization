
[a, b, c] = select_variables(12, 42, 1);

function [cyclicPrefix, tailTx, tailRx] = select_variables( ...
    redundancyLength, seedValue, tryID)
% This function will select randomly values for the system.
%
% - Input:
%   . redundancyLength: Number of samples in cyclic prefix + cyclic suffix
%   in our system
%   . seedValue: a seed value for consistency, if needed
%
% - Output:
%   . cycliPrefix: Number of samples in cyclic prefix
%   . tailTx: Number of samples in rise and fall tail for the transmitter
%   . tailRx: Number of samples in rise and fall tail for the receiver

folderName = 'random_search';
if ~isdir(folderName)  %#ok
    mkdir(folderName)
end
valid = false;
while ~valid
    rng(seedValue)
    cyclicPrefix = randi(redundancyLength);
    
end
end