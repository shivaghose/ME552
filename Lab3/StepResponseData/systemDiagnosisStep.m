function [riseTime settlingTime percentOvershoot steadyStateErrorPercent] = systemDiagnosisStep(stepStart,stepEnd,T,X)
% given the time steps, set point, and response determine the rise time,
% settling time, percent overshoot, and steady state error. 
% OUT = [risetime, settling time, percent overshoot, steady state error percent];
% rise time 90% of final value settling time 2% of final value

% percent overshoot
if (stepEnd - stepStart > 0) % if going up
    % look for the maximum value
    percentOvershoot = 100 * (max(X) - stepEnd) / (stepEnd - stepStart);
else % otherwise
    % look for the minimum value
    percentOvershoot = 100 * (min(X) - stepEnd) / (stepEnd - stepStart);
end

% steady state error
% look at the last 5% of values we have 
index = floor(0.95 * length(X));
finalVal = mean(X(index:length(X)));
steadyStateErrorPercent = 100 * abs((finalVal - stepEnd) / (stepEnd - stepStart));
disp(finalVal)

% determine settling time
for i = length(X): -1 : 1
    if abs((X(i) - finalVal) / finalVal) > 0.02
        index = i + 1; % last index that was within 2% of final value
        break;
    end
end
settlingTime = T(index) - T(1); % change in time

% determine rise time
for i = 1 : length(X)
    if abs((X(i) - finalVal) / finalVal) < 0.1
        index = i; % first index to get witin 10% of final value
        break;
    end
end
riseTime = T(index) - T(1); 

end

