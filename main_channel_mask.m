% main_channel_mask.m
%   This script will simulate a communication system using the channel mask
%   for out-of-band randiation (OOB) reduction, and will compare results
%   from systems using the optimized windows and the RC windows.
%
% Author: Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% Mar 2, 2023
%

clc
clear
close all

fprintf('Starting main_channel_mask.m ... \n')

folderName = 'optimized_windows';
assert(isdir(folderName), 'Missing optimized_windows folder:', ...
    'Run main_window_optimization.m.\n')  %#ok
settingsFileName = 'settingsData.mat';
settingsLoader = load(settingsFileName);
typeOFDM = 'wtx';
tailTx = settingsLoader.settingsData.(typeOFDM).tailTx;
switch typeOFDM
    case 'wtx'
        csLength = tailTx;
        prefixRemovealLength = cpLength;
        circularShiftLength = 0;
    case 'wrx'
end



function windowTxRC = transmitter_rc_window(numSubcar, cpLength, ...
    csLength, tailTx)
% Function to generate the Raised Cosine transmitter window.
%
% - Input:
%   . numSubcar: Number of subcarriers in the system.
%   . cpLength: Number of samples in cyclic prefix.
%   . csLength: Number of samples in cyclic suffix.
%   . tailTx: Number of samples in rise and fall tail for the transmitter
%   window.
%
% - Output:
%   . windowTxRC: Matrix with diagonal equivalent to the transmitted raised
%   cosine window.

raisedCosineAxis = (-(tailTx+1)/2 + 1):1:((tailTx+1)/2 - 1);
raisedCosine = sin(pi/2 * (.5 + raisedCosineAxis/tailTx)).^2;
onesTransmitted = ones(1, numSubcar+cpLength+csLength ...
    - 2*tailTx);
windowVector = [raisedCosine onesTransmitted fliplr(raisedCosine)];
windowTxRC = diag(windowVector);
end


% EoF
