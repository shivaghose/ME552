function [] = experimentalBode(f_min, f_max, num_pts)
%EXPERIMENTALBODE Summary of this function goes here
%   Detailed explanation goes here

%df = 0.01;
%length = 10/df;

dt = 0.0001;
tau = 1;
%freq = dt/tau;

T = 0:dt:50;
F = logspace(f_min,f_max,num_pts);
L = length(F);
%F = zeros(length,1);
A = zeros(L,4);
Phase = zeros(L,4);

i = 1;

%Kalman Filtering
p_init = 1;
steadyState_Data_Variance = 0.07; % calculated from the static case
reduceLag_factor = 100;
%r =  1/(steadyState_Data_Variance * reduceLag_factor);
q = 0.001;
KF_args = [p_init, steadyState_Data_Variance, reduceLag_factor, q];

%Low Pass Filter
tau_LPF = 0.05;
delta_t = dt;
LPF_args = [tau_LPF, delta_t];

%Moving Average Filter
MAvgWindowSize = 100;
MAvgF_args = [MAvgWindowSize];

%High Pass Filter
tau_HPF = 0.05;
delta_t = dt;
HPF_args = [tau_HPF, delta_t];


while i <= L
%for f = df:df:10
    f = F(i);
    data = sin(2*pi*f*T);
    %x0 = data(1);
    %Xfiltered = filter(freq,[1 freq-1],X,x0);
    data_KF = lab6Filters(data', 'KF', KF_args);
    [amp,phase] = estimateSineModified(T,data_KF,f);
    A(i,1) = amp;
    Phase(i,1) = phase;
    
    data_LPF = lab6Filters(data', 'LPF', LPF_args);
    [amp,phase] = estimateSineModified(T,data_LPF,f);
    A(i,2) = amp;
    Phase(i,2) = phase;
    
    data_MAvgF = lab6Filters(data', 'MAvgF', MAvgF_args);
    [amp, phase] = estimateSineModified(T,data_MAvgF,f);
    A(i,3) = amp;
    Phase(i,3) = phase;
    
    data_HPF = lab6Filters(data', 'HPF', HPF_args);
    [amp,phase] = estimateSineModified(T,data_HPF,f);
    A(i,4) = amp;
    Phase(i,4) = phase;
    
    %F(i) = f;
    i = i + 1;
end

% AdB = 20*log(A);
% 
% figure
% hold on
% subplot(2,1,1)
% semilogx(F,AdB);
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (dB)')
% title('Magnitude')
% 
% subplot(2,1,2)
% semilogx(F,Phase);
% xlabel('Frequency (Hz)')
% ylabel('Phase (radians)')
% title('Phase')
% hold off

log_amps = 20*log(A);

freqs = F;

phases = rad2deg(Phase);

figure
hold on
subplot(2,1,1)
plot(freqs',log_amps(:,1))
set(gca,'XScale','log')
KF_title = ['Amplitude Response of Kalman Filter, \sigma = ' num2str(steadyState_Data_Variance) ', q = ' num2str(q)];
title(KF_title,'FontWeight','bold')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

subplot(2,1,2)
plot(freqs',phases(:,1))
set(gca,'XScale','log')
title('Phase Response of Kalman Filter','FontWeight','bold')   
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
hold off

%Plot Low Pass Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,2))
set(gca,'XScale','log')
title(['Amplitude Response of Low Pass Filter, \tau = ' num2str(tau_HPF)],'FontWeight','bold')   
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

subplot(2,1,2)
plot(freqs',phases(:,2))
set(gca,'XScale','log')
title('Phase Response of Low Pass Filter','FontWeight','bold') 
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
hold off

%Plot Moving Avg. Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,3))
set(gca,'XScale','log')
title(['Amplitude Response of Moving Avg. Filter, \mu_{size} = ' num2str(MAvgWindowSize)],'FontWeight','bold')  
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

subplot(2,1,2)
plot(freqs',phases(:,3))
set(gca,'XScale','log')
title('Phase Response of Moving Avg. Filter','FontWeight','bold')   
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
hold off


%Plot Moving Avg. Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,4))
set(gca,'XScale','log')
title(['Amplitude Response of High Pass Filter, \tau = ' num2str(tau_HPF)],'FontWeight','bold')    
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

subplot(2,1,2)
plot(freqs',phases(:,4))
set(gca,'XScale','log')
title('Phase Response of High Pass Filter','FontWeight','bold')   
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
hold off

end

