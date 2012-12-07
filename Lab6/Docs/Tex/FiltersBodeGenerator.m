freqs = logspace(-2,2,1000);
amps = zeros(length(freqs),4);
phases = zeros(length(freqs),4);
timeStepSize = 0.001;
for i = 1:length(freqs)
    if(freqs(i)<=0.01)        
        timeEnd = 200;
    elseif(freqs(i)<=0.1)        
        timeEnd = 20;
    elseif(freqs(i)<=1)        
        timeEnd = 2;
    elseif(freqs(i)<=10)        
        timeEnd = 0.2;
    elseif(freqs(i)<=100)        
        timeEnd = 0.2;
    elseif(freqs(i)<=1000)        
        timeEnd = 0.2;    
    end
    [amps(i,:), phases(i,:)] = bodeGen(freqs(i), timeStepSize, timeEnd);    
end

%amps_ratio = amps;
log_amps = 20*log10(amps);

phases = rad2deg(phases);

%Plot Kalman Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,1))
set(gca,'XScale','log')
title('Amplitude Response of Kalman Filter','FontWeight','bold')    

subplot(2,1,2)
plot(freqs',phases(:,1))
set(gca,'XScale','log')
title('Phase Response of Kalman Filter','FontWeight','bold')        
hold off

%Plot Low Pass Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,2))
set(gca,'XScale','log')
title('Amplitude Response of Low Pass Filter','FontWeight','bold')    

subplot(2,1,2)
plot(freqs',phases(:,2))
set(gca,'XScale','log')
title('Phase Response of Low Pass Filter','FontWeight','bold')        
hold off

%Plot Moving Avg. Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,3))
set(gca,'XScale','log')
title('Amplitude Response of Moving Avg. Filter','FontWeight','bold')    

subplot(2,1,2)
plot(freqs',phases(:,3))
set(gca,'XScale','log')
title('Phase Response of Moving Avg. Filter','FontWeight','bold')        
hold off


%Plot Moving Avg. Filter
figure
hold on    
subplot(2,1,1)
plot(freqs',log_amps(:,4))
set(gca,'XScale','log')
title('Amplitude Response of High Pass Filter','FontWeight','bold')    

subplot(2,1,2)
plot(freqs',phases(:,4))
set(gca,'XScale','log')
title('Phase Response of High Pass Filter','FontWeight','bold')        
hold off