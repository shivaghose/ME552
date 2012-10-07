function [Xt] = interpolate(T,X,t)
% linearly interpolates function values if not exactly one of the values we
% have from simulation
% assume that T monotonically increases
index = inf;
for i = 1:length(T)
    if t <= T(i)
        index = i;
        break;
    end
end

if (index == inf) % if we ran off the end without finding it
    % linearly extrapolate - subject to noise, be careful
    dt = T(end) - T(end - 1);
    slope = (X(end) - X(end - 1)) / dt;
    dt = t - T(end);
    Xt = X(end) + slope * dt;
elseif (t == T(index))
    % we got it exactly
    Xt = X(index);
elseif (t < T(1)) % if off the front
    dt = T(2) - T(1);
    slope = (X(2) - X(1)) / dt;
    dt = t - T(1);
    Xt = X(1) + slope * dt; 
else % linearly interpolate
    dt = T(index) - T(index - 1);
    slope = (X(index) - X(index - 1)) / dt;
    dt = t - T(index - 1);
    Xt = X(index - 1) + slope * dt;
end

