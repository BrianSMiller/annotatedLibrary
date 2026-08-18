function verify_detection_contract()
% verify_detection_contract  Sanity-check the fLow/fHigh<->freq contract
% changes: syncFreqBounds.m itself, its integration into all four
% matchbox methods, and its integration into bsnr's snrEstimateImpl.
%
% Fully synthetic -- no real deployment data, WAV files, or site metadata
% required. Needs annotatedLibrary, matchbox, and bsnr (with bsnr/tests
% on the path, for createTestFixture) all resolvable.
%
%   verify_detection_contract
%
% B. Miller, AAD, 2026

fprintf('\n=== verify_detection_contract ===\n');

%% (1) syncFreqBounds itself -- both directions, plus error cases

fprintf('--- (1) syncFreqBounds ---\n');

t = table([100;200], [150;250], 'VariableNames', {'fLow','fHigh'});
t = syncFreqBounds(t, 'fLow');
assert(isequal(t.freq, [t.fLow, t.fHigh]), ...
    'authority=fLow should set freq = [fLow, fHigh]');
fprintf('  [PASS] authority=fLow derives freq correctly\n');

t2 = table([100 150; 200 250], 'VariableNames', {'freq'});
t2 = syncFreqBounds(t2, 'freq');
assert(isequal(t2.fLow, t2.freq(:,1)) && isequal(t2.fHigh, t2.freq(:,2)), ...
    'authority=freq should set fLow/fHigh from freq columns');
fprintf('  [PASS] authority=freq derives fLow/fHigh correctly\n');

% authority side missing -> should error, not silently do nothing
tBad = table([100;200], 'VariableNames', {'fLow'});  % fHigh missing
threw = false;
try
    syncFreqBounds(tBad, 'fLow');
catch err
    threw = strcmp(err.identifier, 'syncFreqBounds:missingAuthority');
end
assert(threw, 'Missing fHigh with authority=fLow should throw syncFreqBounds:missingAuthority');
fprintf('  [PASS] missing authority column errors loudly, not silently\n');

% Overwrite behaviour: stale opposite side gets clobbered, not preserved
t3 = table([100;200], [150;250], [999 999; 999 999], ...
    'VariableNames', {'fLow','fHigh','freq'});
t3 = syncFreqBounds(t3, 'fLow');
assert(isequal(t3.freq, [t3.fLow, t3.fHigh]), ...
    'authority=fLow should overwrite a stale freq, not defer to it');
fprintf('  [PASS] authority side overwrites a stale opposite side\n');

%% (2) matchbox: all four methods populate freq from their own fLow/fHigh

fprintf('--- (2) matchbox methods populate freq ---\n');

t0base = datenum(2020,1,1);
d1 = table(t0base + [0;10;20]/86400, t0base + [5;15;25]/86400, ...
    [100;100;100], [200;200;200], ...
    'VariableNames', {'t0','tEnd','fLow','fHigh'});
d2 = table(t0base + [1;11;40]/86400, t0base + [6;16;45]/86400, ...
    [105;105;105], [205;205;205], ...
    'VariableNames', {'t0','tEnd','fLow','fHigh'});

methodChecks = { ...
    'gridded',        {'gridStep', 5/86400}; ...
    'clustered',      {}; ...
    'pairwise',       {}; ...
    'pointProximity', {} ...
};

for m = 1:size(methodChecks,1)
    methodName = methodChecks{m,1};
    extraArgs  = methodChecks{m,2};
    ch = matchbox(d1, d2, 'method', methodName, 'verbose', false, extraArgs{:});
    assert(ismember('freq', ch.Properties.VariableNames), ...
        sprintf('%s: output missing freq column', methodName));
    assert(isequal(ch.freq, [ch.fLow, ch.fHigh]), ...
        sprintf('%s: freq does not match [fLow, fHigh]', methodName));
    fprintf('  [PASS] %s: freq present and matches fLow/fHigh (%d events)\n', ...
        methodName, height(ch));
end

%% (3) bsnr: snrEstimateImpl populates fLow/fHigh from its own freq

fprintf('--- (3) snrEstimateImpl populates fLow/fHigh ---\n');

[annot, cleanupFn] = createTestFixture('durationSec', 5, 'freq', [150 250]);
cleanupObj = onCleanup(cleanupFn); %#ok<NASGU>

annotBatch = repmat(annot, 2, 1);
res = snrEstimate(annotBatch, 'snrType', 'spectrogramSlices', ...
    'showClips', false, 'verbose', false);

assert(all(ismember({'fLow','fHigh'}, res.Properties.VariableNames)), ...
    'snrEstimate output missing fLow/fHigh columns');
assert(isequal(res.freq, [res.fLow, res.fHigh]), ...
    'snrEstimate fLow/fHigh do not match freq');
fprintf('  [PASS] snrEstimate batch output has fLow/fHigh matching freq (%d rows)\n', ...
    height(res));

%% (4) resolveObserverColumns: freq re-derived when fLow/fHigh are resolved

fprintf('--- (4) resolveObserverColumns re-syncs freq ---\n');

cht = table( ...
    [true;false], [true;true], ...
    [100;NaN], [200;NaN], ...
    [110;120], [210;220], ...
    [999 999; 999 999], ...
    'VariableNames', {'detect_observer1','detect_observer2', ...
        'fLow_observer1','fHigh_observer1', ...
        'fLow_observer2','fHigh_observer2', 'freq'});

cht = resolveObserverColumns(cht, {'fLow','fHigh'});

assert(isequal(cht.fLow,  [100;120]), 'fLow not resolved with expected observer priority');
assert(isequal(cht.fHigh, [200;220]), 'fHigh not resolved with expected observer priority');
assert(isequal(cht.freq, [cht.fLow, cht.fHigh]), ...
    'freq should be re-synced to the newly-resolved fLow/fHigh, not left stale');
fprintf('  [PASS] freq re-synced after fLow/fHigh resolution (was stale [999 999])\n');

% Resolving an unrelated column should NOT touch freq at all
chtOther = table( ...
    [true;false], [true;true], ...
    {'a';''}, {'';'b'}, ...
    'VariableNames', {'detect_observer1','detect_observer2', ...
        'soundFolder_observer1','soundFolder_observer2'});
chtOther = resolveObserverColumns(chtOther, {'soundFolder'});
assert(~ismember('freq', chtOther.Properties.VariableNames), ...
    'Resolving an unrelated column should not create/touch freq');
fprintf('  [PASS] resolving unrelated columns leaves freq untouched\n');

%% (5) qcAttachSnr: fLow/fHigh stay in sync via the shared contract

fprintf('--- (5) qcAttachSnr re-syncs fLow/fHigh ---\n');

[annotQc, cleanupFnQc] = createTestFixture('durationSec', 5, 'freq', [150 250]);
cleanupObjQc = onCleanup(cleanupFnQc); %#ok<NASGU>

detTable = struct2table(repmat(annotQc, 2, 1));
detTable.fLow  = detTable.freq(:,1);
detTable.fHigh = detTable.freq(:,2);

detTable = qcAttachSnr(detTable, 'snrType','spectrogramSlices', 'nfft',128);

assert(all(ismember({'snr','RL','NL'}, detTable.Properties.VariableNames)), ...
    'qcAttachSnr output missing snr/RL/NL columns');
assert(isequal(detTable.freq, [detTable.fLow, detTable.fHigh]), ...
    'qcAttachSnr: fLow/fHigh not in sync with freq');
fprintf('  [PASS] qcAttachSnr output has fLow/fHigh matching freq (%d rows)\n', ...
    height(detTable));

%% (6) qcAttachSnr + kalmanRidge: fLow/fHigh stay in sync under GENUINE refinement
%
% (5) only exercises the pass-through case (spectrogramSlices never
% refines bounds). This exercises the case that actually motivated
% syncFreqBounds in the first place -- a method that genuinely moves
% t0/tEnd/freq, not just copies them through unchanged.

fprintf('--- (6) qcAttachSnr + kalmanRidge (genuine refinement) ---\n');

% Same fixture params as bsnr's own test_snrEstimate_methods.m kalmanRidge
% section (known there to produce a trackable, non-dropout call) -- loose
% [150 250] band around a 200 Hz tone gives kalmanRidge something to
% genuinely tighten, unlike (5)'s pass-through-only method.
[annotKr, cleanupKr] = createTestFixture('sampleRate', 2000, 'durationSec', 8, ...
    'toneFreqHz', 200, 'signalRMS', 1.0, 'noiseRMS', 0.1, 'freq', [150 250]);
cleanupObjKr = onCleanup(cleanupKr); %#ok<NASGU>

detTableKr = struct2table(repmat(annotKr, 2, 1));
detTableKr.fLow  = detTableKr.freq(:,1);
detTableKr.fHigh = detTableKr.freq(:,2);

% Call the real path directly (qcAttachSnr -> snrEstimate -> snrEstimateImpl)
% rather than pre-checking exist('viterbiKalmanPresenceTrack','file') first --
% that check is unreliable here because it depends on the CALLING script's
% own location for private-folder visibility, not on whether the real
% pipeline can resolve it. snrEstimateImpl itself lives inside bsnr/, so it
% always sees its own private/ siblings regardless of where this script
% runs from; only an exist() check from outside bsnr was ever location-
% sensitive. Genuine unavailability still shows up as an UndefinedFunction
% error below and is skipped, not silently ignored.
try
    detTableKr = qcAttachSnr(detTableKr, 'snrType','kalmanRidge');
    krAvailable = true;
catch err
    if strcmp(err.identifier, 'MATLAB:UndefinedFunction') || ...
       contains(err.message, 'viterbiKalmanPresenceTrack') || ...
       contains(err.message, 'deriveBoundsFromPresenceTrack')
        fprintf('  [SKIP] kalmanRidge unavailable: %s\n', err.message);
        krAvailable = false;
    else
        rethrow(err);
    end
end

if krAvailable
    nRefined = sum(detTableKr.boundsRefined);
    assert(nRefined > 0, sprintf(['kalmanRidge should genuinely refine at least one row ' ...
        'on a high-SNR trackable tone -- got %d of %d refined. Dropout is a legitimate ' ...
        'per-row outcome, but 0 of %d is suspicious for this fixture; check the fixture ' ...
        'is still trackable before assuming this is a real regression.'], ...
        nRefined, height(detTableKr), height(detTableKr)));
    assert(isequal(detTableKr.freq, [detTableKr.fLow, detTableKr.fHigh]), ...
        'qcAttachSnr: fLow/fHigh not in sync with freq under genuine kalmanRidge refinement');
    fprintf('  [PASS] %d of %d rows genuinely refined, fLow/fHigh still match freq\n', ...
        nRefined, height(detTableKr));
end

fprintf('\n=== verify_detection_contract: all checks passed ===\n');

end