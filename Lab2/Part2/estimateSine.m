function [amplitude, phase] = estimateSine(filename, FrequencyHz, startT, endT)
% if we know the frequency component that we are supposed to see can
% estimate the magnitude and phase of the resulting signal with respect to
% it,  expects the first column to be time, 2nd the output

VAL = xlsread(filename);
T = VAL(:,1);
Out = VAL(:,2);

startI = 1; % figure out the correct starting position
while (T(startI) < startT) && (startI < length(T))
    startI = startI + 1;
end


endI = length(T);
while (T(endI) > endT) && (endI > 0)
    endI = endI - 1;
end

% truncate
T = T(startI:endI, 1);
Out = Out(startI:endI, 1);

X = ones(length(T),3);
X(:,2) = cos(2*pi*FrequencyHz*T);
X(:,3) = sin(2*pi*FrequencyHz*T);

beta = X\Out;

amplitude = sqrt(beta(2)^2 + beta(3)^2);

phase = atan2(-beta(3), beta(2)); % might not be correct example was for cossine

phase = phase + pi / 2; % going from cossine to sine

% S = amplitude * sin(2*pi*FrequencyHz*T + phase);
% 
% M = mean(Out);
% for i = 1:length(S)
%     S(i) = S(i) + M;
% end
% 
% figure;
% hold on;
% plot(T,Out,'k');
% plot(T,S,'b');
% hold off;


end

