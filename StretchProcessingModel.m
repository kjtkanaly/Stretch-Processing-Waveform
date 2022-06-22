close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot(s) Decleration
figure;
sb1 = subplot(2,1,1);
hold(sb1,'on');
grid on;
xlabel("Pulse Width (t/T)",'FontSize',14)

sb2 = subplot(2,1,2);
hold(sb2,'on');
grid on;
xlabel("Range (m)",'FontSize',14)
ylabel("dB",'FontSize',14)
ylim([0 60])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Details
c	= physconst('lightspeed');		% Speed of Light (m/s)

includeNoise	= false;

% Radar Parameters
Radar.T				= 64e-6;			% LFM Period (s)
Radar.fs				= 0.8e9;			% First Sample Frequency (Hz)
Radar.fs2				= 10e6;			% Second Sampling Frequency (Hz)
Radar.B				= 500e6;			% LFM Bandwidth (Hz)
Radar.TauR			= 0;					% Estimated time delay of Target (s)

% Target Details
Target.A           = [0 3 6 14 16 14]-20;	% Scatterers' SNR (dB)
Target.rngOff	= [-10 -5 0 2 5 5.4];		% Scatterers' Range


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Model
[signal, time] = generateStretchProcessingSignal(Radar, Target, includeNoise);

% Decimate the Signal
D				= Radar.fs/Radar.fs2;
signalDC	= decimate(signal.',D,'fir');
timeDC    = downsample(time, D);

% FFT Parameters
nFFT	= 8192;
Fs		= 1/((1/Radar.fs) * D);
freq	= (-nFFT/2:nFFT/2-1)*(Fs/nFFT);

% Apply an FFT to the TIme-Domain Signal
Signal = fftshift(fft(signalDC,nFFT));

% High Range Res Vector
chirp_rate	= Radar.B/Radar.T;         
rng	= (c*freq)./(2.*chirp_rate);

% Plotting the Resutls
plot(sb1, timeDC./Radar.T, real((signalDC(:,1))),'linewidth',1.5)
plot(sb2, rng, mag2db(abs((Signal(:,1)))), 'linewidth',1.5)
sgtitle("Target Response using a " + Radar.B/10^6 + " MHz Bandwidth and Range Delay of " + (c*Radar.TauR/2) + " m",'FontSize',14)
xlim([-20 20])

%% Functions

function [signal, time] = generateStretchProcessingSignal(Radar,Target,includeNoise)
	
	% Speed of Light (m/s)
	c	= physconst('lightspeed');

	% Radar Parameters
	tauP				= Radar.T;		% Period of the LFM (s)
	fs					= Radar.fs;	% Sampling Rate (Hz)
	B					= Radar.B;	% LFM Bandwidth
	chirpRate		= B/tauP;		% LFM Chirp Slope

	% Time Parameters
	ts				= 1/fs;												% Sample Size (s)
	Ns			= floor(tauP/ts);								% Size of the time vector (n)
	time			= linspace(-tauP/2,tauP/2,Ns);		% Time Vector (s)
	
	% Target Parameters
	A			= Target.A;											% Amplitude of the scatterers
	rngOff	= Target.rngOff;									% Range of the scatterers
	tauR	= 2.*rngOff./c - Radar.TauR;				% Time-Delay of the scatterers

	signal = sum(db2mag(A.') .* (cos(2 .* pi .* chirpRate .* (tauR.') .* time) + 1j.*sin(2 .* pi .* chirpRate .* (tauR.') .* time)), 1);

	% Check to Add Noise
	if includeNoise == true
		noise = 1/sqrt(2)*complex(randn(1,length(signal)),randn(1,length(signal)));
		signal = signal + noise;
	end

end

