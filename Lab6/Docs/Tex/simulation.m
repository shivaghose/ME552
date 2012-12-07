function [] = simulation()
% simulates system in vertical configuration
SensY = 0.03400293;
g = 9.81;

%meanVY = 1.95268; % = V_Ybias + SensY * g from the actual data
%VYbias =  meanVY - SensY * g; % calculated from the actual data
VYbias = 1.617433958;
meanVY = VYbias + SensY * g;

%sigma2VY = 5.807739 * 10^(-5); % this is the measured number
sigma2VY = 1.24553 * 10^(-5); % number from training set

EA = 0; EA2 = 0;
ET = 0; ET2 = 0;
cnt = 1000000;

for i = 0:cnt
    Vy = normrnd(meanVY,sqrt(sigma2VY));
    Arg = (Vy - VYbias)/(SensY * g);
    
    if Arg > 1
        Arg = 1;
    elseif Arg < -1
        Arg = -1;
    end
   
    Theta = acos(Arg);
    
    ET = ET + Theta;
    ET2 = ET2 + Theta*Theta;
    
    EA = EA + Arg;
    EA2 = EA2 + Arg*Arg;
    
end

meanTheta = ET / cnt

varianceTheta = ET2 / cnt - meanTheta * meanTheta

meanA = EA / cnt

varianceA = EA2 / cnt - meanA * meanA

end

