function channelArray = channel_array(channelVector, dftLength, ...
    cpLength, csLength, tailTx, tailRx, prefixRemovalLength)
% CHANNEL_ARRAY     Generates the channel array for the w-OFDM problem,
%                   represented as H in the paper.
%   channelArray = CHANNEL_ARRAY(channelVector, dftLength, cpLength,
%       csLength, tailTx, tailRx, prefixRemovalLength) .
%
%   - Parameters:
%       . channelVector : Tapped delay line channel model.
%       . dftLength : Length of DFT in OFDM system.
%       . cpLength : Length of cyclic prefix.
%       . csLength : Length of cyclic suffix.
%       . tailTx : Length of rise and fall tails at transmitter.
%       . tailRx : Length of rise and fall tails at receiver.
%       . prefixRemovalLength : Number of samples to be removed from the
%           beginning of the received symbol.
%   - Returns:
%       . channelArray : Array representing the channel effects.

chanOrder = length(channelVector) - 1;
numRx = dftLength+tailRx+prefixRemovalLength;
numTx = dftLength+cpLength+csLength;
numAffected = ceil((chanOrder+tailTx)/numRx);
chanNoise = numTx-tailTx;
channelArray = zeros(numRx, numTx, numAffected+1);
for symbAffected = 0:numAffected
    for rxSample = 0:numRx-1
        for txSample = 0:numTx-1
            indexer = symbAffected*chanNoise + rxSample - txSample;
            if (0 <= indexer) && (indexer <= chanOrder)
                channelArray(rxSample+1, txSample+1, symbAffected+1) ...
                    = channelVector(indexer+1);
            end
        end
    end
end
end


% EoF
