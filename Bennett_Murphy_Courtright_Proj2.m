close all
clc
clear

% Load FIR filter in b
load("FIRFilter.mat")

% Uncomment if you want plots
%PLOTS = ones(1, 15);

% Comment out if you don't want plots
PLOTS = zeros(1, 15);

% To get filter response
impulse = @(n) (n==0);

% Values for different channels
band = 1837;
halfBand = band/2;

% Import and set variables
[y, Fs] = audioread('proj_2.wav');
%sound(y, Fs)
y=y.';

% Variables for Time, Samples and Frequency
t = (0:(length(y)-1)) / Fs;
n = 0:length(y)-1;
numPts = length(y);
f = (-numPts/2 : numPts/2-1)*Fs/numPts;

y_fft = fftshift(fft(y)/numPts);

if (PLOTS(1))
    % Plot signal in time and frequency
    figure
    subplot(2,1,1)
    plot(t, y, 'Linewidth', 1.5)
    title('Raw Signal :: Time Representation')
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    grid on
    subplot(2,1,2)
    plot(f, abs(y_fft), 'Linewidth', 1.5)
    title('Raw Signal :: Frequency Representation')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (V/Hz)')
    grid on
end

% Filter Response
FIR = filter(b, 1, impulse(t));
FIR_fft = fftshift(fft(FIR)/numPts);

if (PLOTS(2))
    % Plot impulse response
    figure
    subplot(2,1,1)
    plot(t, FIR, 'Linewidth', 1.5)
    title('Low Pass Filter Response :: Time Representation')
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    xlim([0 0.025])
    grid on
    subplot(2,1,2)
    plot(f, abs(FIR_fft), 'Linewidth', 1.5)
    title('Low Pass Filter Response :: Frequency Representation')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (V/Hz)')
    xlim([-1020 1020])
    grid on
end

outWave = [];

for k=0:5
    
    % Bandpass filter the desired channel
    bandPass = FIR .* cos(2*pi*n*(halfBand+k*band)/Fs);
    bandPass_fft = fftshift(fft(bandPass)/numPts);
    
    %{
    % THIS CODE JUST SHOWS THE FILTER MOVING IN FREQUENCY

    figure
    subplot(2,1,1)
    plot(t, bandPass, 'Linewidth', 1.5)
    title('Bandpass Filter :: Time Representation')
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    xlim([0 0.025])
    grid on
    subplot(2,1,2)
    plot(f, abs(bandPass_fft), 'Linewidth', 1.5)
    title('Bandpass Filter :: Frequency Representation')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (V/Hz)')
    %xlim([-1020 1020])
    grid on
    %}
    
    % Filter everything else out, Get rid of the negative side band 
    x_f = bandPass_fft .* y_fft;
    x_HT = hilbert(ifft(ifftshift(x_f)));
    
    if (PLOTS(3+k))
        % Plot stuff
        figure
        hold on
        plot(f, abs(x_f), 'Linewidth', 1.5)
        plot(f, abs(fftshift(fft(x_HT))), 'Linewidth', 1.5)
        hold off
        title(['Channel ', num2str(k+1), ' :: Frequency Domain'])
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (V/Hz)')
        legend('Filtered Signal', 'Hilbert Transform')
        grid on
    end
    
    
    % Shift to baseband and take real part
    x_n_lp = x_HT.*exp(-j*2*pi*n*(k)*band/Fs);
    x_n = real(x_n_lp);
    
    if (PLOTS(9+k))
        figure
        subplot(2,1,1)
        plot(t, x_n, 'Linewidth', 1.5)
        title(['Channel ', num2str(k+1), ' at Baseband :: Time Domain'])
        xlabel('Time (s)')
        ylabel('Amplitude (V)')
        grid on
        subplot(2,1,2)
        plot(f, abs(fftshift(fft(x_n))), 'Linewidth', 1.5);
        title(['Channel ', num2str(k+1), ' at Baseband :: Frequency Domain'])
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (V/Hz)')
        grid on
    end
    
    outWave = [outWave x_n./max(x_n)];
end

t2 = (0:(length(outWave)-1)) / (Fs);

if (PLOTS(15))
    figure
    plot(t2, outWave)
    title('Final Output :: Time Domain')
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    grid on
end

%soundsc(outWave, Fs)
audiowrite('Bennett_Murphy_Courtright_Proj2.wav', outWave, Fs);