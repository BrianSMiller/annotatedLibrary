function channel = decodePamguardChannelMap(channelMap)
%decodePamguardChannelMap  Convert PAMGuard's channel bitmask to a
% 1-based channel number.
%
% PAMGuard's channelMap field is a bitmask: bit N set means channel N+1
% (0-based bit, 1-based channel) contributed to the detection. For the
% single-channel deployments this pipeline currently handles, exactly
% one bit is set per detection.
%
%   channel = decodePamguardChannelMap(channelMap)
%
% Errors loudly if any row has zero or more than one bit set, rather
% than silently picking one -- multi-channel handling is not implemented
% in this pipeline yet, so surfacing it as a hard failure here is
% correct until it's needed.
%
% B. Miller, AAD, 2026

if any(channelMap == 0)
    error('decodePamguardChannelMap:noChannel', ...
        '%d row(s) have channelMap == 0 -- no channel bits set.', ...
        sum(channelMap == 0));
end

nBitsSet = arrayfun(@(x) sum(dec2bin(x) == '1'), channelMap);
if any(nBitsSet > 1)
    error('decodePamguardChannelMap:multiChannel', ...
        '%d row(s) have more than one channel bit set in channelMap. Multi-channel detections are not handled by this pipeline yet.', ...
        sum(nBitsSet > 1));
end

channel = log2(double(channelMap)) + 1;
end