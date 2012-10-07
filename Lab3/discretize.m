function [out] = discretize(fun, t, resolution)
% funciton which simulates an encoder applied to the function fun at time t
% fun return two values the angular position at time t, and the angular
% velocity
% note that resolution is in pulses per revolution

% get the angle and angular velocity in radians, counter clockwise is
% positive
[ft, fdott] = fun(t);

% convert resolution to radians
resolution = resolution / (2 * pi);

% number of pulses which at this point is still a decimal
pulses = ft * resolution;

if fdott > 0 % if rotating counter-clockwise, then round down to nearest pulse
    out = floor(pulses) / resolution;
elseif fdott < 0 % if going clockwise, then round up to nearest pulse
    out = ceil(pulses) / resolution;
else % the number of pulses is not well determined
    % hysteriesis tells us what this should be but that is getting a bit
    % tricky
    out = round(pulses) / resolution;
end
end

