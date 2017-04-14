% This function compensates chromatic dispersion with OFDE algorithm.
% It is my implementation of Blind look-up filter.
% See article: Tianhua Xu (and others), "Chromatic dispersion compensation
% in coherent transmission system using digital filters"
% 
% inCompIxQx / inCompIyQy - vertical vectors of complex samples
% samplingRateIn - ADC sampling rate
% lambdaRef - reference laser lambda [m]
% span - fiber length [m]
% D - fiber chromatic dispersion coeff [s/m^2]
% S - chromatic dispersion slope [s/m^3]
% ffLength - chosen FFT length, depends of dispersion (see article)
% L - total overlap length, usually ffLength/2
%
% Copyright by Michal Drzymala, 2017
function [outCompIxQx, outCompIyQy, samplingRateOut] = CDE_OFDE(inCompIxQx, inCompIyQy,...
	samplingRateIn, lambdaRef, span, D, S, fftLength, L)
	%inCompIxQx Ix+i*Qx
	%inCompIyQy Iy+i*Qy

	c = 299792458;
	
	% Truncate filter length if necessary
	NSamples = size(inCompIxQx);
	if fftLength > NSamples(1)
		fftLength = NSamples(1);
	end
	
	% Create frequency grid
	fc = c / lambdaRef;
	Tf = fftLength / samplingRateIn; 
	df = 1 / Tf;
	fg = df * (-fftLength/2 : fftLength/2-1);

	% Create filter coeffs from transfer function
	H_D = -1j*D*span*pi*c / fc^2 * fg.^2;
	H_S = 1j*S*span*pi*c^2 / 3 / fc^4 * fg.^3;
	H = exp(H_D + H_S);

	outCompIxQx = OverlapBothTrans(inCompIxQx, H, L);
	outCompIyQy = OverlapBothTrans(inCompIyQy, H, L);

	% outReal4D = Ix | Qx | Iy | Qy
	%outReal4D = ...
	%    [real(outCompIxQx) imag(outCompIxQx) real(outCompIyQy) imag(outCompIyQy)];
	samplingRateOut = samplingRateIn;
end

% This function convolves together signal x(t) with the filter transfer function H(f)
% and truncates output signal y to length(x) with Overlapp-Both method.
% It's modification of Overlap-Save method taking into account overlap from
% both sides.
%
% x - input signal (time domain, real or complex, vector 1x(N) or (N)x1)
% H - filter transfer function (frequency domain, real or complex, vector 1x(N) or (N)x1)
% L - signal chunks length
%
% FFT length = N = length(H) = L + B (should be power of 2)
% where B is the overlap length (B/2 previous samples and B/2 next samples)
%
% Copyright by Michal Drzymala, 2017
function y = OverlapBothTrans(x, H, L)
	y = [];

	% x must be one-dimensional vector
	xSizeVector = size(x);
	if xSizeVector(1) > 1 && xSizeVector(2) > 1;
		display('Error: x must be one dimensional complex vector'); return; end

	% H must be one-dimensional vector of even length
	HsizeVector = size(H);
	if HsizeVector(1) > 1 && HsizeVector(2) > 1;
		display('Error: H must be one dimensional vector'); return; end
	if mod(length(H), 2) ~= 0;
		display('Error: H must be even length'); return; end

	% L must > 0 and shorter than filter length
	if L <= 0;
		display('Error: L must be > 0'); return; end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	if L > length(H);
		display('Error: L must be shorter than filter length'); return; end

	% Length of input signal must be greater or equal filter length
	if length(x) < length(H);
		display('Error: Signal must be longer or equal filter'); return; end

	% Vectorize (make vectors horizontal, .' is complex-safe transposition)
	x = x(:).';
	H = H(:).';

	% Resize x to be divisible by L
	nx = length(x);
	m = mod(nx, L);
	if m ~= 0; x = [x zeros(1, L-m)]; end

	N = length(H);
	B = N-L;
	nxL = length(x);
	y = zeros(1, length(x));
	% Extend x with leading and trailing zeros
	B2 = B/2;
	x = [zeros(1, B2) x zeros(1, B2)];

	for i = 1:L:nxL
		% Select (next) x chunk of length N
		xc = x(i:i+N-1);
		% Calc N-point FFT of the x chunk
		Xc = fftshift(fft(xc, N));
		% Filter signal (by multiplication of the transforms)
		Yc = Xc.*H;
		% Calc y chunk with the IFFT
		yc = ifft(ifftshift(Yc));
		% Skip next L samples (instead of N samples) and add middle L samples
		%   of yc to the output signal (discarding the aliasing of B samples)
		y(i:i+L-1) = y(i:i+L-1) + yc(1+B2:end-B2);
	end

	% Truncate output signal to the length of the input signal
	y = y(1:nx);

	% Reshape y if x was column vector
	if xSizeVector(1) > 1
		y = y(:);
	end
end