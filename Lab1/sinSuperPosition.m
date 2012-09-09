function [ ] = sinSuperPosition(a, A, frequency, phase)
%a = gain for each term
% A = amplitude of each sine wave
% frequency of each in Hz
% phase of each wave

% sampling rate determined 
maxF = max(frequency);
sF = 5 * maxF;
sP = 1 / sF;

T = 0: sP : 10; % simulate 0 to 10 seconds

Sum = zeros(length(T), 1);

for i = 1: length(T)
    for j = 1 : length(a)
        value = a(j) * A(j) * sin(frequency(j) * T(i) + phase(j));
        Sum(i) = Sum(i) + value;
    end
end

plot(T,Sum);
xlabel('Time (seconds)');
ylabel('Sum of sinusoids');
title('Sum of sinusoids');

window = length(Sum); % window length
n = pow2(nextpow2(window));
y = fft(Sum, n);
f = (0: n-1) * (sF / n);
power = y.*conj(y) / n;

figure;
plot(f,power)
xlabel('Frequency (Hz)');
ylabel('Power');
title('Periodogram');
end

