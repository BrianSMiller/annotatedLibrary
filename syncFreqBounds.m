function t = syncFreqBounds(t, authority)
% syncFreqBounds  Keep fLow/fHigh and freq mutually consistent on a
% detection table, deriving whichever's stale from whichever's
% authoritative. Part of the shared minimal contract detection tables
% moving between matchbox, bsnr, and annotatedLibrary are expected to
% satisfy (t0, tEnd, fLow, fHigh, freq, soundFolder, channel, duration --
% add whatever else a stage needs, just don't break these). Every caller
% names the direction explicitly rather than this function guessing
% which of fLow/fHigh or freq is more current.
%
%   t = syncFreqBounds(t, 'fLow')  % fLow/fHigh just set -> freq derived
%   t = syncFreqBounds(t, 'freq')  % freq just set       -> fLow/fHigh derived
%
% Inputs
%   t         - table with fLow/fHigh and/or freq columns
%   authority - 'fLow' or 'freq'. Whichever you just computed/set --
%               that side wins, the other side is overwritten to match.
%
% Output
%   t - same table, with fLow, fHigh, and freq all present and mutually
%       consistent.
%
% Errors if the requested authority's own column(s) don't exist -- there
% has to be something to derive FROM. Does not error or warn if the
% derived side already existed and gets overwritten; that's the whole
% point of naming an authority.
%
% B. Miller, AAD, 2026

switch lower(authority)
    case 'flow'
        if ~all(ismember({'fLow','fHigh'}, t.Properties.VariableNames))
            error('syncFreqBounds:missingAuthority', ...
                'authority=''fLow'' requires both fLow and fHigh columns to be present.');
        end
        t.freq = [t.fLow, t.fHigh];

    case 'freq'
        if ~any(strcmp('freq', t.Properties.VariableNames))
            error('syncFreqBounds:missingAuthority', ...
                'authority=''freq'' requires a freq column to be present.');
        end
        t.fLow  = t.freq(:,1);
        t.fHigh = t.freq(:,2);

    otherwise
        error('syncFreqBounds:unknownAuthority', ...
            'authority must be ''fLow'' or ''freq'', got ''%s''.', authority);
end

end
