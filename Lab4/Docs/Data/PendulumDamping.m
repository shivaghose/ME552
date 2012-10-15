clear all;
close all;
clc;

%Test 1: A330:A48732
%Test 2: A330:A48179
%Test 3: A330:A50300
%Test 4: A250:A42400
%Test 5: A300:A49600
%Test 6: A300:A52400

[time,txt,raw] = xlsread('PendDamp6.xlsx','A300:A52400');
[amp,txt,raw] = xlsread('PendDamp6.xlsx','B300:B52400');

m = 0.0820142;
g = 9.81;
rBar = 0.0857;  % distance from axis of rotation to CG

s=size(time);



k=1;
for i=1:100:s(1)-100
    j=1;
    peakA(k) = amp(i);
    for j=1:100
        if(amp(j+i) >= peakA(k))
            peakA(k) = amp(j+i);
            peakT(k) = time(j+i);            
        end
        
    end
    k = k+1;
end

s1 = size(peakT);
for i=2:s1(2)-1
    if(peakT(i)==0)
        peakT(i) = (peakT(i-1)+peakT(i+1))/2;
    end
end

peakT = peakT.';
peakA = peakA.';

plot(peakT(1:5:s1(2)), peakA(1:5:s1(2)), 'o')

fitType = 'exp1';
myFit = fit(peakT(1:150), peakA(1:150), fitType)

hold on;

plot(peakT, myFit(peakT),'k');

for i=2:s1(2)
    freq(i-1) = (peakT(i)-peakT(i-1))^-1;
end
for i=5:5:150
    delta = log(peakA(i-4)/peakA(i))/4;
    zeta(i/5) = 1/sqrt(1+(2*pi/delta)^2);
end

dampedFreq = mean(freq)*2*pi
Zeta = mean(zeta)
natFreq = dampedFreq/sqrt(1-Zeta^2)

-natFreq * Zeta
plot(peakT, peakA(1)*exp(-natFreq*Zeta*peakT));
J = m*g*rBar/natFreq^2
B = 2*Zeta*natFreq*J

figure
plot(peakT(s1(2)-50:s1(2)),peakA(s1(2)-50:s1(2)), 'o')

j=1;
for i=s1(2)-49:s1(2)
    tf(j) = m*g*rBar/4*(peakA(i-1)-peakA(i));
    j= j+1;
end

Tf_piecewise = mean(tf)
hold on;
p=polyfit(peakT(s1(2)-50:s1(2)),peakA(s1(2)-50:s1(2)),1);

plot(peakT(s1(2)-50:s1(2)), p(1)*peakT(s1(2)-50:s1(2))+p(2),'k');

Tf_slope = -.5*p(1)*pi/natFreq*m*g*rBar