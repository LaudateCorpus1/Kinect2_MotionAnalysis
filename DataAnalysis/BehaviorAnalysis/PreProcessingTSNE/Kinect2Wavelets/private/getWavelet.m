function [ amp, W, freqs, scales, t] = getWavelet(x, Fs, fMin, fMax, nbins, omega0)
% [ amp, W, freqs, scales, t] = getWavelet(x, Fs, fMin, fMax, nbins, omega0)
%
% x:        
% vector - wavelet analysis is applied
% multi-column matrix - wavelet analysis is applied to each column
% separately
%
% Fs - sampling frequency
% fMin, fMax, nbins - parameters for spectral analysis - min/max freqency
% (Hz) and number of frequency bins (Morlet is always used so this is
% frequency per se and not just scale)
% omega0 - non-dimensional parameter for Morlet wavelet
%
% output:
% amp    - the PSD
% W      - wavelet coefficients (complex valued)
% freqs  - vector of wavelet frequencies
% scales - wavelet scales, for computing frequencies
% t      - vector of timeline @ Fs
%
% call: wavelet
%
% 08-oct-12 ES
%
% revisions
% 30-oct-14 jdl2: I stripped the routine down to a basic wrapped around the 
%                 subroutines of Berman et al. 2013. I'm using the previous
%                 code just as a template

if isa( x, 'int16' ), x = single( x ); end


sx = size( x );

if min( sx ) == 1
    nChannels = 1;
    nsamples  = max( sx );
    % make into column vector
    x = x( : );
else
    % Assumes column matrix
    nChannels = sx( 2 );
    nsamples  = sx( 1 );
end

% Sampling Rate
dt      = 1/Fs; 
minT    = 1/fMax;
maxT    = 1/fMin;

% Freqs to match processing of Berman et al. 2014
Ts      = minT.*2.^((0:nbins-1).*log(maxT/minT)/(log(2)*(nbins-1)));
freqs   = fliplr(1./Ts);
scales  = (omega0 + sqrt(2+omega0^2))./(4*pi.*freqs);
t       = ( 1 : nsamples )' / Fs;

% actually compute
%[wave,period,scale,coi] = wavelet(Y,dt,pad,dj,s0,J1,mother,param);
amp = zeros( nbins, nsamples, nChannels );
W   = zeros( nbins, nsamples, nChannels );
for i = 1 : nChannels
    [ amp( :, :, i ), W( :, :, i ) ] = waveletMorletConvolution( x( :, i ), freqs, omega0, dt);
end

return
