% This script takes in a data sequence ,modulates it using the BPSK
% modulation technique and plays it over a computer speaker, simultaneously
% recording the played sound throught the computer's microphone after which
% the signal is decoded to get back data message
clc
clear

%% Step 1: Generate Data in Binary
data = 'Rockwell was here';
text = char(data);                   %Converts character to ASCII                                        
y = zeros(length(data)*8,1);         %Converts  data from decimal to binary                                  
for n = 1:1:length(data)  
    a=abs(data(n));                                                              
    f = 8*(n-1)+1;                                                          
    y(f:f+7,1)=(de2bi(a,8));                                                   
end
data_binary = y';
  
%% Step 2: Generate Preamble
preamble = randi(10, 16, 1); % Creating the preamble
preamble_binary = de2bi(preamble); %Converts  preamble from decimal to binary
preamble_binary = reshape(preamble_binary, 1, size(preamble_binary,1)*size(preamble_binary,2));

%% Step 3: Define Variables
fs = 44100; % Sample rate
fc = fs/40; % Carrier frequency
fb = fc/50; %Bit rate

Nsb = fs/fb; % Number of samples per bit
T = 10; 

%% Step 4: Apply Direct BPSK to Preamble
t = (1/fs):(1/fs):T; % time for the carrier wave 

Training_Sequence = sin(2*pi*fc*t); %carrier wave

for counter = 1 : size(preamble_binary,2)
    if preamble_binary(counter) == 1 %Modulation of the carrier wave with preamble (Differential Modulation)
        Training_Sequence(counter*Nsb) = Training_Sequence(counter*Nsb)*1 ;
    elseif preamble_binary(counter) == 0
        Training_Sequence(counter*Nsb) = Training_Sequence(counter*Nsb)* -1;
        for count = 1:Nsb
            Training_Sequence(counter*Nsb + count) = Training_Sequence(counter*Nsb + count)* -1;
        end
    end
end

Training_Sequence = Training_Sequence(1:size(preamble_binary,2)*Nsb);

%% Step 5: Apply Differential BSPK to Payload

Data_Sequence = sin(2*pi*fc*t); %Carrier wave

for counter = 1 : size(data_binary,2)
    if data_binary(counter) == 1 %Modulation of the carrier wave with data (BPSK)
        Data_Sequence(counter*Nsb) = Data_Sequence(counter*Nsb) * -1;
        for count = Nsb*counter:size(Data_Sequence,2)
            Data_Sequence(count) = Data_Sequence(count)* -1;
        end
    elseif data_binary(counter) == 0
        for count = Nsb*counter:size(Data_Sequence,2)
            Data_Sequence(count) = Data_Sequence(count)* 1;
        end
    end
end

Data_Sequence = Data_Sequence(1:size(data_binary,2)*Nsb);

tm = t(1 : (size(preamble_binary,2) + size(data_binary,2))*Nsb); %Modulated time function

mod_wave = [Training_Sequence, Data_Sequence]; % Modulated wave carrying preamble and data

figure(1)
plot(tm, mod_wave)
axis([2.035 2.044 -1.05 1.05])
xlabel('Time (Seconds)')
ylabel('Amplitude')
title('Amplitude Vs Time showing a change of phase')

pkFFT = fft(mod_wave);
psd = 10*log10(pkFFT .* conj(pkFFT));
freqList = linspace(-fs/2,fs/2, length(pkFFT));

figure(2)
plot(freqList,fftshift(psd))
axis([fc-fb*4 fc+fb*4 (max(psd)-60) max(psd)+5])
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Unaveraged Spectra of Generated Data')

nAv = 5;
pkFFT = fft(reshape(mod_wave, length(mod_wave)/nAv, nAv));
matpsd = (pkFFT .* conj(pkFFT));
psd = 10*log10(sum(matpsd,2));
freqList = linspace(-fs/2, fs/2, size(pkFFT,1));

figure(3)
plot(freqList,fftshift(sum(psd,2)))
grid on
axis([fc-fb*4 fc+fb*4 (max(psd)-60) max(psd)+5])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Averaged Spectra of Generated Data')

%% Step 6: Create Objects for Playing and Recording

ap = audioplayer(mod_wave, fs);
ar = audiorecorder(fs, 16, 1);

%% Step 7: Play and Record Modulated data

record(ar, T + 1)
pause(1)
play(ap)
pause(T + 1)
stop(ar)
% play(ar)

recorded_data = getaudiodata(ar,'double');
recorded_data = recorded_data';

tr = (1/fs):(1/fs):T + 1;

figure(4)
plot(tr,recorded_data)
axis([0 T+2 -1.05 1.05]);
xlabel('Time (Seconds)')
ylabel('Amplitude')
title('Recorded Sound')

figure(5)
plot(tr,recorded_data)
axis([4.7 4.91 -1.05 1.05]);
xlabel('Time (Seconds)')
ylabel('Amplitude')
title('Recorded Sound')

pkFFT = fft(recorded_data);
psd = 10*log10(pkFFT .* conj(pkFFT));
freqList = linspace(-fs/2,fs/2, length(pkFFT));

figure(6)
plot(freqList,fftshift(psd))
axis([0 7000 30 90])
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Spectra of Recorded Data (Wide)')

nAv = 5;
pkFFT = fft(reshape(recorded_data, length(recorded_data)/nAv, nAv));
matpsd = (pkFFT .* conj(pkFFT));
psd = 10*log10(sum(matpsd,2));
freqList = linspace(-fs/2, fs/2, size(pkFFT,1));

figure(7)
plot(freqList,fftshift(sum(psd,2)))
grid on
axis([fc-fb*4 fc+fb*4 20 90])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Averaged Signal of Recorded Data (Narrow)')

%% Step 8: Mix Data Down to Baseband

pkFFT = fft(recorded_data);
psd = 10*log10(pkFFT .* conj(pkFFT));
freqList = linspace(-fs/2,fs/2, length(pkFFT));

figure(8)
plot(freqList,fftshift(psd))
axis([-4000 4000 30 90])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Averaged Signal of Recorded Data (Double Sideband)')

ts = (1/fs):(1/fs): T + 1;

lo_inphase = cos(2*pi*fc*ts);
lo_quadrature = sin(2*pi*fc*ts); 

baseband_data = lo_inphase .* recorded_data + 1i*lo_quadrature .* recorded_data;

Training_Sequence_Zeros = [Training_Sequence , zeros(1, size(recorded_data,2) - size(Training_Sequence,2))];
baseband_preamble = lo_inphase .* Training_Sequence_Zeros + 1i*lo_quadrature .* Training_Sequence_Zeros;

pkFFT = fft(baseband_data);
psd = 10*log10(pkFFT.* conj(pkFFT));
freqList = linspace(-fs/2,fs/2, length(pkFFT));

figure(9)
plot(freqList,fftshift(psd))
axis([-4000 4000 30 90])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Spectra After Mixing with Complex fc')

%% Step 9: Apply Filtering
 
Coeff = ones(1, round(0.5 * fs/fb))/round(0.5 * fs/fb);

filtered_baseband_data = filter(Coeff, 1, baseband_data);
filtered_baseband_preamble = filter(Coeff, 1, baseband_preamble);

pkFFT = fft(filtered_baseband_data);
psd = 10*log10(pkFFT.* conj(pkFFT));
freqList = linspace(-fs/2,fs/2, length(pkFFT));
 
figure(10)
plot(freqList,fftshift(psd))
axis([-4000 4000 30 90])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/ Arbitary Reference)')
title('Spectra After Filtering')

nAv = 5;
pkFFT = fft(reshape(filtered_baseband_data, length(filtered_baseband_data)/nAv, nAv));
matpsd = (pkFFT .* conj(pkFFT));
psd = 10*log10(sum(matpsd,2));
freqList = linspace(-fs/2, fs/2, size(pkFFT,1));

figure(11)
plot(freqList,fftshift(psd))
axis([-70 70 30 90])
xlabel('Frequency (Hz)')
ylabel(' Power Spectral Density (dB/ Arbitary Reference)')
title('Spectra After Filtering')

figure(12)
plot(ts, abs(filtered_baseband_data))
axis([0 10 0 .1])
xlabel('Time (sec)')
ylabel('Complex Data Magnitude')
title('Magnitude of Complex Filtered Baseband Data')

figure(13)
plot(ts, angle(filtered_baseband_data))
axis([0 T+2 -4 4])
xlabel('Time (sec)')
ylabel('Complex Data Phase Angle (Radians)') 
title('Phase of Complex Filtered Baseband Data')
        
%% Step 10: Correlate Recorded Data with Training Sequence

Correlatoroutput = ifft(conj(fft(filtered_baseband_preamble)) .* fft(filtered_baseband_data));

figure(14)
plot(ts, abs(Correlatoroutput))
axis([0 T+2 0 14000])
xlabel('Time (sec)')
ylabel('Correlator Output Magnitude')
title('Magnitude of Correlator Output Vs Time')

[max_value, start_sample] = max(Correlatoroutput);

%% Step 11: Extract Samples From Recroded Data

end_sample = start_sample + (Nsb* (size(data_binary,2) + size(preamble_binary,2))) - 1;
extracted_packet = filtered_baseband_data(1, start_sample : end_sample);

tn = ts(start_sample : end_sample);

figure(15)
plot(tn, abs(extracted_packet))
%axis([0 T+2 0 .45])
xlabel('Time (sec)')
ylabel('Complex packet Sample Magnitude')
title('Complex packet Sample Magnitude Vs Time')

figure(16)
plot(tn, angle(extracted_packet))
%axis([0 T+2 -4 4])
xlabel('Time (sec)')
ylabel('Complex Packet Data Phase Angle (Radians)') 
title('Phase of Complex Packet Data Vs Time')

%% Step 12: Reduce Total Number of Samples

sample_bits = extracted_packet(1000:Nsb:end);
max_sample = max(sample_bits);
sample_bits_norm = sample_bits./max_sample;

figure(17)
scatter(real(sample_bits_norm), imag(sample_bits_norm))
axis([-1 1 -1 1])
xlabel('Real Value')
ylabel('Imaginary Value (Normalized)') 
title('Imaginary Vs Real Value')

%% Step 13: Plot the Phase

angle_sample_data_bits = angle(sample_bits(size(preamble_binary,2):1:end));%this will include the first preamble bit so that diff() reamoves that
bit_number = 1:1:(size(data_binary,2));
difference = diff(angle_sample_data_bits);

for z = 1 : 1 : size(difference,2)
    if (difference(z) > pi/2 || difference(z) < -pi/2)
        demod_data(z) = 1;
    else%if (difference(z) < 0)
        demod_data(z) = 0;
    end
end
demod_data(1) = [];
demod_data(120) = 0;
figure(18)
plot(bit_number, angle_sample_data_bits(2:1:end))
axis([0 size(data_binary,2) -4 4])
xlabel('Bit Number')
ylabel('Phase of Bit Sample (Radians)') 
title('Phase of Bit Samples Vs Real Value (Data Only)')

%% Step 14: Print the Received Data

out = char(bi2de(reshape(demod_data,8,15).'))







