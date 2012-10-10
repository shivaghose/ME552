function [Xt] = flatInterpolate(T,X,t)
% not linear interpolation, but rather just returns the previous value
% makes sense for when interpolating the current values

index = inf;
for i = 1:length(T)
    if t < T(i)
        index = i;
        break;
    end
end

if (index == inf)
    Xt = X(end);
elseif (index == 1)
    Xt = 0; % assume initial value of 0
else
    Xt = X(index - 1);
end

end

