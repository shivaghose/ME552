function [T,Uout,Lout] = extractEnvelop(T,X,window)
% extracts lower and upper Envelopes of the given presumably sinusoidal
% function
assert((window > 0),'Window should have been positive');

% get upper and lower
[TUE, UE] = extractUpperEnvelop(T,X,window);
[TLE, LE] = extractLowerEnvelop(T,X,window);

Uout = interpolate(TUE,UE,T);
Lout = interpolate(TLE,LE,T);

end

function [TE, UE] = extractUpperEnvelop(T,X,w)
UEsub = [];
Tsub = [];
k = 1;
if (X(2) - X(1)) > 0
    increasing = 1;
    UEsub(1) = X(1);
    Tsub(1) = T(1);
else 
    increasing = 0;
end
for i = 2:length(X)
    % if going up, keep going
    if (increasing) && (X(i) > UEsub(k)) 
        UEsub(k) = X(i);
        Tsub(k) = T(i);
    elseif (increasing) % as soon as we stop increasing, stop updating that one
        k = k + 1;
        increasing = 0;
    end
    
    if (~increasing) && (isIncreasing(X,i,w)) % if we weren't going up, now we are
        increasing = 1;
        UEsub(k) = X(i);
        Tsub(k) = T(i);
    end
end

% if we were still increasing when we quit, then get rid of the last point,
% since it isn't representative
if increasing
   UE = UEsub(1:(length(UEsub)-1));
   TE = Tsub(1:(length(Tsub)-1));
else
    UE = UEsub;
    TE = Tsub;
end

end

function [B] = isIncreasing(X,index,window)
% is the function increasing at that particular timestep

% look backwards if we can
if index - window < 1
    start = 1;
    finish = 1 + window;
else
    start = index - window;
    finish = index;
end

B = 1;
for i = (start + 1):finish
    if (X(i) - X(i-1)) <= 0
        B = 0;
        break;
    end
end

end

function [B] = isDecreasing(X,index,window)
% look backwards if we can
if index - window < 1
    start = 1;
    finish = 1 + window;
else
    start = index - window;
    finish = index;
end

B = 1;
for i = (start + 1):finish
    if (X(i) - X(i-1)) >= 0
        B = 0;
        break;
    end
end
end

function [TE, LE] = extractLowerEnvelop(T,X,w)
LEsub = [];
Tsub = [];
k = 1;
if (X(2) - X(1)) < 0
    decreasing = 1;
    LEsub(1) = X(1);
    Tsub(1) = T(1);
else 
    decreasing = 0;
end
for i = 2:length(X)
    % if going up, keep going
    if  (decreasing) && (X(i) < LEsub(k))
        LEsub(k) = X(i);
        Tsub(k) = T(i);
    elseif (decreasing) % as soon as we stop increasing, stop updating that one
        k = k + 1;
        decreasing = 0;
    end
    
    if (~decreasing) && (isDecreasing(X,i,w)) % if we weren't going up, now we are
        decreasing = 1;
        LEsub(k) = X(i);
        Tsub(k) = T(i);
    end
end

if decreasing
   LE = LEsub(1:(length(LEsub)-1));
   TE = Tsub(1:(length(Tsub)-1));
else
    LE = LEsub;
    TE = Tsub;
end
    
end

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