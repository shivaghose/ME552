function [tout] = mod2pi(t)
% returns an angle between pi and -pi

if t < 0
    tout = -mod2pi_sub(-t);
else
    tout = mod2pi_sub(t);
end
end

function [tout] = mod2pi_sub(t)
% only good for positive numbers
tout = t - floor(t*0.5/pi + 0.5) * 2*pi;
end
