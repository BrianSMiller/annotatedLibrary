function f = abwTonalFreq(t)
% f = abwTonalFreq(t)
% Lookup the frequency, f, of Antarctic blue whale unit A at a given time,
% t.
%
% The lookup function is the long-term linear fit from the paper Gavrilov
% et al 2012. This paper documented very precisely that the tonal frequency
% of unit A of Antarctic blue whale song is decreasing linearly as a
% function of time. This function only fits the long-term linear trend, and
% does not attempt to estimate seasonal variation.

t0 = datenum([2002.195975503062 zeros(1,5)]);
p = [-0.135/365, 27.665908913585564]; % Linear fit from Gavrilov et al 2012
daysSinceStart =  t - t0;
f = polyval(p,daysSinceStart);