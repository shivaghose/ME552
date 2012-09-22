function [F] = estimateForceDumb(C,I,x)
% really dumb force estimate

F = - (C * I * I) / (x * x);

end

