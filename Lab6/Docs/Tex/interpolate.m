% function [Xt] = interpolate(T,X,t)
% % linearly interpolates function values if not exactly one of the values we
% % have from simulation
% % assume that T monotonically increases
% index = inf;
% for i = 1:length(T)
%     if t <= T(i)
%         index = i;
%         break;
%     end
% end
% 
% if (index == inf) % if we ran off the end without finding it
%     % linearly extrapolate - subject to noise, be careful
%     dt = T(end) - T(end - 1);
%     slope = (X(end) - X(end - 1)) / dt;
%     dt = t - T(end);
%     Xt = X(end) + slope * dt;
% elseif (t == T(index))
%     % we got it exactly
%     Xt = X(index);
% elseif (t < T(1)) % if off the front
%     dt = T(2) - T(1);
%     slope = (X(2) - X(1)) / dt;
%     dt = t - T(1);
%     Xt = X(1) + slope * dt; 
% else % linearly interpolate
%     dt = T(index) - T(index - 1);
%     slope = (X(index) - X(index - 1)) / dt;
%     dt = t - T(index - 1);
%     Xt = X(index - 1) + slope * dt;
% end
function [Xout] = interpolate(Tin,Xin,Tout)
% linear interpolation of the desired Tout values

Xout = zeros(length(Tout),1);

i = 1; % keep track of positon in the original
for o = 1:length(Tout)
    while (i <= length(Tin)) && (Tout(o) > Tin(i))
        i = i + 1;
    end
    if (i == 1) % if we asked for a T before the first T
        dt = Tin(2) - Tin(1);
        slope = (Xin(2) - Xin(1)) / dt;
        dt = Tout(o) - Tin(1);
        Xout(o) = Xin(1) + slope * dt;
    elseif (i > length(Tin))
        % linearly extrapolate - subject to noise, be careful
        dt = Tin(end) - Tin(end - 1);
        slope = (Xin(end) - Xin(end - 1)) / dt;
        dt = Tout(o) - Tin(end);
        Xout(o) = Xin(end) + slope * dt;
    else
        dt = Tin(i) - Tin(i - 1);
        slope = (Xin(i) - Xin(i - 1)) / dt;
        dt = Tout(o) - Tin(i - 1);
        Xout(o) = Xin(i - 1) + slope * dt;
    end
end

end
