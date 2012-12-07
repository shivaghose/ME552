function [MeanZeroRate, VarianceZeroRate] = ZeroRate(Encoder, Rate)
% given encoder measurements, figures out where the first sign of motion is
% and determines the zero rate in Volts
% gives ZeroRate in V and variance in V^2
% in this case we are looking for exactly 0

start = 1;

while (start < length(Encoder)) && (abs(Encoder(start)) == 0)
    start = start + 1;
end

%start
%length(Encoder)

% now just average over those elements
SumRate = 0;
SumRate2 = 0;
for i = 1:start
    SumRate = SumRate + Rate(i);
    SumRate2 = SumRate2 + Rate(i)*Rate(i);
end

MeanZeroRate = SumRate / start;

VarianceZeroRate = SumRate2 / start - MeanZeroRate * MeanZeroRate;
end

