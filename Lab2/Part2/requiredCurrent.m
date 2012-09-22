function [I] = requiredCurrent(A,B,x)
% how much current do we need to hold it up at this distance according to
% this model
global f;
I = sqrt(-f * (B+x)*(B+x) / A);

end

