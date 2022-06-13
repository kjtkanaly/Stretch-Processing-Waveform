close all;
clear all;

% Simulation Details
includeNoise	= false;

% Radar Parameters
Radar.T               = 64e-6;       % LFM Period
Radar.fs              = 0.8e9;       % First Sample Frequency
Radar.fs2            = 10e6;        % Second Sampling Frequency
Radar.B              = 70e6;         % LFM Bandwidth

% Target Details
Target.A           = [0 3 6 14 16 14]-20; % snr in dB
Target.rngFnt   = 10e3;
Target.rngOff	= [-10 -5 0 2 5 5.4]; % rng offset in meters wrt center of gate

[signalDC, Signal] = StretchProcessingWaveform(Radar, Target, includeNoise, 4);


function [signalDC, Signal] = StretchProcessingWaveform(radar,target,includeNoise, figureNumber)
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Figure Decleration
	figure(figureNumber);
	sb1 = subplot(2,1,1);
	hold(sb1,'on');
	grid on;
	xlabel("Pulse Width (\mus)",'FontSize',14)

	sb2 = subplot(2,1,2);
	hold(sb2,'on');
	grid on;
	xlabel("Range (m)",'FontSize',14)
	ylabel("dB",'FontSize',14)
	xlim([(target.rngFnt-70) (target.rngFnt + 70)])
	ylim([0 60])

	sgtitle("Target Response using a " + radar.B/10^6 + " MHz Bandwidth",'FontSize',14)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Parameter(s) Decleratrion
	c	= physconst('lightspeed');		% (m/s) - Speed of Light

	% Radar Parameters
	T					= radar.T;
	fs					= radar.fs;
	fs2				= radar.fs2;
	pulses			= radar.pulses;
	B					= radar.B;
	D					= fs/fs2;
	chirp_rate	= B/T;         
	tauM			= 0; 

	% Time Parameters
	ts				= 1/fs;
	Ns			= floor(T/ts);
	time			= linspace(-T/2,T/2,Ns);

	% FFT Parameters
	nFFT	= 8192;
	Fs		= 1/(ts * D);
	freq	= (-nFFT/2:nFFT/2-1)*(Fs/nFFT);
	rng	= (c*freq)./(2.*chirp_rate);
	
	% Target Parameters
	A			= target.A;
	rngFnt	= target.rngFnt;
	rngOff	= target.rngOff;
	V			= target.V;
	T0		= 2.*rngFnt./c;
	tauR	= 2.*rngOff./c;
	
	% Initialize Signal Vectors/Matrices
	Signal = zeros(nFFT, pulses);
	timeDC	= downsample(time, D);
	signalDC	= zeros(length(timeDC), pulses);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculations
	for pCnt = 1:pulses

		% Scatter(s) Responses Post-DeMixing, Eq (12.19) - Basic Radar Analysis
		signal = sum(db2mag(A.') .* exp(1j .* pi .* chirp_rate .* (tauM.^2 - (tauR.').^2)) .* exp(1j .* 2 .* pi .* chirp_rate .* (tauR.' - tauM) .* time), 1);
	
		% Check to Add Noise
		if includeNoise == true
			noise = 1/sqrt(2)*complex(randn(1,length(signal)),randn(1,length(signal)));
			signal = signal + noise;
		end
	
		% Decimate the Signal
		signalDC(:,pCnt) = decimate(signal.',D,'fir');
	
		% Apply an FFT to the TIme-Domain Signal
		Signal(:,pCnt) = fftshift(fft(signalDC(:,pCnt),nFFT));
	
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plotting the Resutls
	plot(sb1, timeDC./T, real((signalDC(:,1))),'linewidth',1.5)
	plot(sb2, rng + rngFnt, mag2db(abs((Signal(:,1)))), 'linewidth',1.5)

end