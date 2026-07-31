function t = f20pTableToDetection(f20pTable, soundFolder, siteCode, classification, freqBand)
%t = f20pTableToDetection(f20pTable, soundFolder, siteCode, classification, freqBand)
% Convert a detectFin20Hz.m output table (File, DateTime, Kurt, KurtProd,
% SNRF, SNRT, BW) into the canonical annotatedLibrary detection format,
% matching kooguTableToDetection.m / ravenTableToDetection.m: original
% detector-specific columns are kept, and t0/tEnd/duration/fLow/fHigh/
% freq/channel/soundFolder/siteCode/classification are added alongside.
%
% t0/tEnd are the same +/-0.6s window already computed internally for
% SNRT (computeF20pParamsOneFile.m: start_signal/end_signal around the
% TKEO peak), NOT a per-call measured extent -- this detector reports a
% single instant per pulse, not a time-frequency box, so a fixed +/-0.6s
% window is a deliberate simplification. Near file start/end, the actual
% signal window used internally is edge-clipped shorter than 1.2s; this
% converter doesn't reproduce that clipping, so boundary detections will
% have a slightly generous box relative to what SNRT actually saw. Minor,
% and only affects the handful of detections within 0.6s of a file edge.
%
% fLow/fHigh are the detection band (freqBand), not a measured frequency
% extent -- this detector doesn't measure per-call bandwidth the way a
% CNN box detector does. IMPORTANT: pass the same freqBand you actually
% used in detectFin20Hz -- there's no way for this function to check that
% for you, and a mismatch here would silently mislabel every detection's
% frequency bounds.
%
% classification defaults to 'Bp-20', NOT '' -- unlike
% kooguTableToDetection/ravenTableToDetection, where '' falls through to
% a per-row Tags passthrough or 'none'. This detector only ever targets
% one call type, so defaulting to 'none' would silently mislabel every
% row for anyone who omits the 4th argument. Pass '' explicitly if you
% want that fallback behaviour instead (e.g. if a Tags column gets added
% to this table by some future step).

if nargin < 5
    freqBand = [15, 26];
end
if nargin < 4
    classification = 'Bp-20';
end
if nargin < 3
    siteCode = '';
end
if nargin < 2
    soundFolder = '';
end

t = f20pTable;

nDetect = height(t);
if (nDetect == 0)
    return;
end

%% Create XBAT-like fields to facilitate compatibility with XBAT 'events'
if any(strcmp(t.Properties.VariableNames, 'Channel'))
    t.channel = t.Channel;
else
    t.channel = ones(nDetect, 1);
end

t.t0 = datenum(t.DateTime) - 0.6/86400;
t.tEnd = datenum(t.DateTime) + 0.6/86400;
t.duration = (t.tEnd - t.t0) * 86400; % seconds; should be ~1.2 for every row by construction

t.fLow = repmat(freqBand(1), nDetect, 1);
t.fHigh = repmat(freqBand(2), nDetect, 1);
t.freq = [t.fLow, t.fHigh];

t.soundFolder = cellstr(repmat(soundFolder, nDetect, 1));
t.siteCode = cellstr(repmat(siteCode, nDetect, 1));

if isempty(classification)
    if any(strcmp('Tags', t.Properties.VariableNames))
        t.classification = cellstr(string(t.Tags));
    else
        t.classification = repmat({'none'}, nDetect, 1);
    end
else
    t.classification = cellstr(repmat(classification, nDetect, 1));
end

t = sortrows(t, 't0');

end
